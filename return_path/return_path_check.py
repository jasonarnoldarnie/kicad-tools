#!/usr/bin/env python3
"""Return-path checker for KiCad PCBs.

For every selected signal net this estimates how far the return current has to
travel compared with the forward (trace) path, and flags the places where the
two diverge:

  * VOID      - trace runs over a hole/gap in its reference plane; the return
                current detours around it. Detour is measured as the shortest
                path through plane copper between the entry and exit points
                (A* on a raster of the zone fill), minus the trace length over
                the gap.
  * ISLAND    - entry and exit copper are not connected on the reference
                layer at all: no local return path.
  * SPLIT     - the reference plane changes net under the trace (e.g. GND pour
                to +3V3 pour); return must go through a capacitor.
  * TERMINAL  - trace ends (pad/via) over a void longer than --terminal-tol.
  * VIA       - signal changes layers and the reference plane changes. If both
                planes are the same net, the nearest stitching via/PTH pad of
                that net is found (added return ~ 2 x distance). If they are
                different nets, the nearest capacitor bridging them is found.
  * EDGE      - trace is referenced, but the plane edge is within
                w/2 + k*h of its centreline (k = --edge-k, h = dielectric to
                the reference plane), so the return current is crowded.

Per net it reports trace length, estimated return length, their ratio, and an
estimated loop area:
  referenced trace   length x h (return directly under the trace)
  void crossing      area between the trace and the detour path (plan view)
  layer change       distance to stitching via x plane-to-plane spacing

Nets are classified (rf, clock, high-speed, analog, medium, low-speed, static)
from their names; override with a JSON file mapping net-name globs to classes
(default: return_path_classes.json next to the board). Issues on low-speed
nets are capped at WARN and on static nets at INFO (--no-class-severity to
disable).

Reference planes are taken from the board stackup: each copper layer is
referenced to the nearest (by dielectric thickness) copper layer that carries
zone fill. Override with --ref F.Cu=In1.Cu.

Zones are refilled in memory before analysis (the board file is never
written). The script reads the board from disk, so save in pcbnew first.

Run with KiCad's bundled Python (it provides pcbnew and numpy):

  "%LOCALAPPDATA%\\Programs\\KiCad\\10.0\\bin\\python.exe" return_path/return_path_check.py \\
      path/to/board.kicad_pcb --html return_path.html

  ... --netclass MS_50R --netclass DP_90R     only controlled-impedance nets
  ... --net '/Micro/*' --include-power        glob on net names
  ... --json report.json --fail               machine output, exit 1 on FAIL
  ... --serve                                 live report on http://127.0.0.1:8765 with a Rerun button

Coordinates are board (page) coordinates in mm, as shown in pcbnew's status
bar with the default origin.
"""
import argparse
import datetime
import fnmatch
import heapq
import json
import math
import os
import re
import sys
from collections import defaultdict

import numpy as np
import pcbnew

NM = 1e6  # nm per mm

POWER_RE = re.compile(r"^[+-]|VBUS|VDD|VCC|VSYS|V_SYS|BATT|VBAT|VIN\b|VOUT", re.I)

# Net classes, most to least sensitive. First matching rule wins; rules match
# the signal part of the name (last path component, or the pin of Net-(U1-PIN)).
CLASSES = ["rf", "clock", "high-speed", "analog", "medium", "low-speed", "static", "unclassified"]
CLASS_CAP = {"low-speed": "WARN", "static": "INFO"}
CLASS_RULES = [
    ("rf", r"ANT|LNA|RF"),
    ("clock", r"XTAL|OSC|MCLK|CLKOUT|^CLK_?\d*$"),
    ("high-speed", r"USB.*D[+-PN]$|^D[+-]$|SD[0-3]$|SCK|/CLK$|/CMD$|^CMD$|MOSI|MISO|SPI|ETH|MDI|HDMI|LVDS|MIPI"),
    ("medium", r"LRCK|BCLK|SCLK|WS$|I2S|DSDIN|ASDOUT|DIN$|DOUT$|PDM"),
    ("analog", r"MIC|LIN\d|RIN\d|LOUT|ROUT|VMID|VREF|HP_?[LR]|SPK|AIN|AOUT|RING|TIP|SLEEVE|OUT$|CAP\d"),
    ("static", r"(^|_)EN\d*$|LED|GPIO|ADC|NTC|CC\d$|CHG|PGOOD|SW$|RST|RESET|BOOT|ILIM|ISET|ITERM|TMR|^CE$|STAT|INT$|SHIELD"),
    ("low-speed", r"UART|TXD?$|RXD?$|I2C|SDA|SCL$|DTR|RTS|CTS|MTCK|MTDI|MTDO|MTMS|JTAG|SWD|TCK|TMS|TDI|TDO"),
]
SEV_RANK = {"FAIL": 0, "WARN": 1, "INFO": 2}

# The HTML report draws the predicted return current as a glow reaching out to
# w/2 + k*h from the trace centreline: with the current density in the plane
# falling off as 1/(1 + (x/h)^2), about 80 % of it flows within 3h.
# Capped so traces far from their plane don't produce huge glows.
RETURN_BAND_K = 3.0
RETURN_BAND_MAX = 2.0


# --------------------------------------------------------------------------
# Stackup

def _sexpr_block(text, start):
    depth = 0
    for j in range(start, len(text)):
        c = text[j]
        if c == "(":
            depth += 1
        elif c == ")":
            depth -= 1
            if depth == 0:
                return text[start:j + 1]
    return text[start:]


def _parse_sexpr(block):
    tokens = re.findall(r'\(|\)|"(?:[^"\\]|\\.)*"|[^\s()]+', block)
    stack = [[]]
    for tok in tokens:
        if tok == "(":
            stack.append([])
        elif tok == ")":
            node = stack.pop()
            stack[-1].append(node)
        else:
            stack[-1].append(tok.strip('"'))
    return stack[0][0]


def parse_stackup(path):
    """Return [(layer_name, is_copper, thickness_mm)] top to bottom, or None."""
    text = open(path, encoding="utf-8").read()
    i = text.find("(stackup")
    if i < 0:
        return None
    tree = _parse_sexpr(_sexpr_block(text, i))
    out = []
    for node in tree[1:]:
        if not (isinstance(node, list) and node and node[0] == "layer"):
            continue
        name = node[1]
        kind = next((n[1] for n in node if isinstance(n, list) and n[0] == "type"), "")
        thick = sum(float(n[1]) for n in node if isinstance(n, list) and n[0] == "thickness")
        out.append((name, kind == "copper" or name.endswith(".Cu"), thick))
    return out


def copper_z(board, stackup):
    """copper layer name -> depth (mm) below the top copper, dielectrics only."""
    if not stackup:
        copper = [board.GetLayerName(l) for l in board.GetEnabledLayers().CuStack()]
        return {n: 0.2 * i for i, n in enumerate(copper)}
    z, out = 0.0, {}
    for name, is_cu, thick in stackup:
        if is_cu:
            out[name] = z
        elif "dielectric" in name:
            z += thick
    return out


def reference_map(board, stackup, planar_layers, overrides):
    """signal layer name -> (reference layer name, h_mm)."""
    copper = [board.GetLayerName(l) for l in board.GetEnabledLayers().CuStack()]
    refs = {}
    if stackup:
        seq = [s for s in stackup if s[1] or "dielectric" in s[0]]
        cu_idx = [i for i, s in enumerate(seq) if s[1]]
        for i in cu_idx:
            best = None
            for direction in (-1, 1):
                h, j = 0.0, i + direction
                while 0 <= j < len(seq) and not seq[j][1]:
                    h += seq[j][2]
                    j += direction
                if 0 <= j < len(seq) and seq[j][0] in planar_layers:
                    if best is None or h < best[1]:
                        best = (seq[j][0], h)
            if best:
                refs[seq[i][0]] = best
    else:  # no stackup: assume adjacent copper, 0.2 mm
        for i, name in enumerate(copper):
            for j in (i - 1, i + 1):
                if 0 <= j < len(copper) and copper[j] in planar_layers:
                    refs[name] = (copper[j], 0.2)
                    break
    for ov in overrides:
        sig, ref = ov.split("=")
        h = refs.get(sig, (None, 0.2))[1]
        refs[sig] = (ref, h)
    return refs


# --------------------------------------------------------------------------
# Net classification

def signal_part(net):
    m = re.match(r"^Net-\((.+?)-(.+)\)$", net)
    s = m.group(2) if m else net.rsplit("/", 1)[-1]
    s = s.replace("{slash}", "/")
    return re.sub(r"~\{([^}]*)\}", r"\1", s)


def classify(net, overrides):
    for glob, cls in overrides.items():
        if fnmatch.fnmatchcase(net, glob):
            return cls
    part = signal_part(net).upper()
    for cls, rule in CLASS_RULES:
        if re.search(rule, part):
            return cls
    return "unclassified"


def load_class_overrides(path):
    if not path or not os.path.exists(path):
        return {}
    with open(path, encoding="utf-8") as f:
        data = json.load(f)
    data = data.get("nets", data)
    bad = {v for v in data.values() if v not in CLASSES}
    if bad:
        sys.exit(f"{path}: unknown class(es) {sorted(bad)}; use one of {CLASSES}")
    return {k: v for k, v in data.items() if not k.startswith("_")}


# --------------------------------------------------------------------------
# Plane rasters

def poly_contours(polyset):
    out = []
    for i in range(polyset.OutlineCount()):
        chains = [polyset.Outline(i)] + [polyset.Hole(i, h) for h in range(polyset.HoleCount(i))]
        for ch in chains:
            n = ch.PointCount()
            if n >= 3:
                out.append(np.array([(ch.CPoint(k).x, ch.CPoint(k).y) for k in range(n)], float) / NM)
    return out


def unfractured(polyset):
    """Copy of a zone fill with its holes restored as holes.

    KiCad stores fills "fractured": each hole is joined to the outline by a
    zero-width horizontal cut. The cuts don't change the copper (and the raster
    skips horizontal edges), but drawn as outlines they show up as long lines."""
    p = pcbnew.SHAPE_POLY_SET(polyset)
    try:
        p.Unfracture()
    except TypeError:  # older KiCad versions take a POLYGON_MODE argument
        p.Unfracture(pcbnew.SHAPE_POLY_SET.PM_FAST)
    return p


class Plane:
    """Raster of zone fill on one layer; each cell holds the net code (0 = none)."""

    def __init__(self, origin, shape, res):
        self.x0, self.y0 = origin
        self.res = res
        self.grid = np.zeros(shape, dtype=np.int32)

    def add(self, contours, value):
        if not contours:
            return
        ex0, ey0, ex1, ey1 = np.vstack([np.hstack([c, np.roll(c, -1, axis=0)]) for c in contours]).T
        ny, nx = self.grid.shape
        res, x0 = self.res, self.x0
        for r in range(ny):
            y = self.y0 + (r + 0.5) * res
            m = (ey0 <= y) != (ey1 <= y)
            if not m.any():
                continue
            xs = np.sort(ex0[m] + (y - ey0[m]) * (ex1[m] - ex0[m]) / (ey1[m] - ey0[m]))
            for xa, xb in zip(xs[0::2], xs[1::2]):
                c0 = max(int(math.ceil((xa - x0) / res - 0.5)), 0)
                c1 = min(int(math.floor((xb - x0) / res - 0.5)), nx - 1)
                if c1 >= c0:
                    self.grid[r, c0:c1 + 1] = value

    def cell(self, x, y):
        return int((y - self.y0) / self.res), int((x - self.x0) / self.res)

    def centre(self, r, c):
        return (self.x0 + (c + 0.5) * self.res, self.y0 + (r + 0.5) * self.res)

    def label(self, x, y):
        r, c = self.cell(x, y)
        if 0 <= r < self.grid.shape[0] and 0 <= c < self.grid.shape[1]:
            return int(self.grid[r, c])
        return 0

    def label_near(self, x, y, radius):
        """Nearest plane net within radius (skips antipads). Returns net code or 0."""
        r, c = self.cell(x, y)
        k = int(radius / self.res) + 1
        r0, c0 = max(r - k, 0), max(c - k, 0)
        win = self.grid[r0:r + k + 1, c0:c + k + 1]
        hits = np.argwhere(win > 0)
        if not len(hits):
            return 0
        d2 = (hits[:, 0] + r0 - r) ** 2 + (hits[:, 1] + c0 - c) ** 2
        return int(win[tuple(hits[int(np.argmin(d2))])])


# 16-connected moves (max ~2.7 % over true Euclidean) with corner-cut checks.
_MOVES = []
for _dr, _dc in [(1, 0), (-1, 0), (0, 1), (0, -1)]:
    _MOVES.append((_dr, _dc, 1.0, ()))
for _dr, _dc in [(1, 1), (1, -1), (-1, 1), (-1, -1)]:
    _MOVES.append((_dr, _dc, math.sqrt(2), ((_dr, 0), (0, _dc))))
for _dr, _dc in [(2, 1), (2, -1), (-2, 1), (-2, -1)]:
    _MOVES.append((_dr, _dc, math.sqrt(5), ((_dr // 2, 0), (_dr // 2, _dc))))
for _dr, _dc in [(1, 2), (-1, 2), (1, -2), (-1, -2)]:
    _MOVES.append((_dr, _dc, math.sqrt(5), ((0, _dc // 2), (_dr, _dc // 2))))


def _astar(mask, s, g):
    """Return (length_in_cells, [cells]) or None."""
    ny, nx = mask.shape
    gr, gc = g
    best = {s: 0.0}
    parent = {}
    heap = [(math.hypot(s[0] - gr, s[1] - gc), 0.0, s)]
    while heap:
        _, d, cur = heapq.heappop(heap)
        if cur == g:
            path = [cur]
            while cur in parent:
                cur = parent[cur]
                path.append(cur)
            return d, path[::-1]
        if d > best.get(cur, math.inf):
            continue
        r, c = cur
        for dr, dc, cost, via in _MOVES:
            nr, nc = r + dr, c + dc
            if not (0 <= nr < ny and 0 <= nc < nx) or not mask[nr, nc]:
                continue
            if any(not mask[r + a, c + b] for a, b in via):
                continue
            nd = d + cost
            if nd < best.get((nr, nc), math.inf):
                best[(nr, nc)] = nd
                parent[(nr, nc)] = cur
                heapq.heappush(heap, (nd + math.hypot(nr - gr, nc - gc), nd, (nr, nc)))
    return None


def _simplify(cells):
    """Drop cells where the step direction does not change."""
    if len(cells) < 3:
        return cells
    out = [cells[0]]
    for a, b, c in zip(cells, cells[1:], cells[2:]):
        if (b[0] - a[0], b[1] - a[1]) != (c[0] - b[0], c[1] - b[1]):
            out.append(b)
    out.append(cells[-1])
    return out


def geodesic(plane, net, p, q):
    """Shortest path through `net` copper on `plane` from p to q.

    Returns (length_mm, [(x, y), ...]) or None if not connected."""
    res = plane.res
    grid = plane.grid
    ny, nx = grid.shape
    straight = math.dist(p, q)
    for margin in (2.0, 6.0, 20.0, None):
        if margin is None:
            r0, c0, r1, c1 = 0, 0, ny - 1, nx - 1
        else:
            m = margin + straight / 2
            r0, c0 = plane.cell(min(p[0], q[0]) - m, min(p[1], q[1]) - m)
            r1, c1 = plane.cell(max(p[0], q[0]) + m, max(p[1], q[1]) + m)
            r0, c0 = max(r0, 0), max(c0, 0)
            r1, c1 = min(r1, ny - 1), min(c1, nx - 1)
        mask = grid[r0:r1 + 1, c0:c1 + 1] == net
        ends = []
        for x, y in (p, q):
            r, c = plane.cell(x, y)
            r, c = r - r0, c - c0
            if not (0 <= r < mask.shape[0] and 0 <= c < mask.shape[1] and mask[r, c]):
                hits = np.argwhere(mask[max(r - 4, 0):r + 5, max(c - 4, 0):c + 5])
                if not len(hits):
                    return None
                r, c = hits[0][0] + max(r - 4, 0), hits[0][1] + max(c - 4, 0)
            ends.append((int(r), int(c)))
        found = _astar(mask, ends[0], ends[1])
        if found is not None:
            d = found[0] * res
            # A path that leaves the window is at least ~2*margin long, so a
            # shorter result is globally optimal.
            if margin is None or d <= 2 * margin + straight:
                pts = [plane.centre(r + r0, c + c0) for r, c in _simplify(found[1])]
                return d, [p] + pts + [q]
    return None


# --------------------------------------------------------------------------
# Geometry helpers

def mm(v):
    return (v.x / NM, v.y / NM)


def via_width(via):
    try:
        return via.GetWidth(pcbnew.F_Cu) / NM
    except TypeError:
        return via.GetWidth() / NM


def shoelace(pts):
    a = 0.0
    for (x0, y0), (x1, y1) in zip(pts, pts[1:] + pts[:1]):
        a += x0 * y1 - x1 * y0
    return abs(a) / 2


def rnd(pts):
    return [[round(x, 3), round(y, 3)] for x, y in pts]


def rdp(pts, tol):
    """Ramer-Douglas-Peucker: drop polyline points within tol (mm) of the line."""
    if len(pts) < 3:
        return list(pts)
    keep = [False] * len(pts)
    keep[0] = keep[-1] = True
    stack = [(0, len(pts) - 1)]
    while stack:
        i, j = stack.pop()
        (x0, y0), (x1, y1) = pts[i], pts[j]
        dx, dy = x1 - x0, y1 - y0
        length = math.hypot(dx, dy)
        worst, idx = tol, None
        for k in range(i + 1, j):
            px, py = pts[k]
            d = (abs(dy * (px - x0) - dx * (py - y0)) / length if length
                 else math.hypot(px - x0, py - y0))
            if d > worst:
                worst, idx = d, k
        if idx is not None:
            keep[idx] = True
            stack += [(i, idx), (idx, j)]
    return [p for p, k in zip(pts, keep) if k]


def sample_track(t, step):
    """Return (points[(x,y)], tangents[(tx,ty)], ds) along a track or arc."""
    a, b = mm(t.GetStart()), mm(t.GetEnd())
    if t.GetClass() == "PCB_ARC":
        m = mm(t.GetMid())
        ax, ay = a
        bx, by = m
        cx, cy = b
        dd = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
        if abs(dd) > 1e-12:
            ux = ((ax**2 + ay**2) * (by - cy) + (bx**2 + by**2) * (cy - ay) + (cx**2 + cy**2) * (ay - by)) / dd
            uy = ((ax**2 + ay**2) * (cx - bx) + (bx**2 + by**2) * (ax - cx) + (cx**2 + cy**2) * (bx - ax)) / dd
            rad = math.dist((ux, uy), a)
            a0 = math.atan2(ay - uy, ax - ux)
            am = math.atan2(by - uy, bx - ux)
            a1 = math.atan2(cy - uy, cx - ux)
            sweep = (a1 - a0) % (2 * math.pi)
            if (am - a0) % (2 * math.pi) > sweep:  # clockwise
                sweep -= 2 * math.pi
            length = abs(sweep) * rad
            n = max(1, math.ceil(length / step))
            pts, tans = [], []
            s = 1 if sweep > 0 else -1
            for i in range(n + 1):
                ang = a0 + sweep * i / n
                pts.append((ux + rad * math.cos(ang), uy + rad * math.sin(ang)))
                tans.append((-math.sin(ang) * s, math.cos(ang) * s))
            return pts, tans, length / n
    length = math.dist(a, b)
    n = max(1, math.ceil(length / step))
    tx, ty = ((b[0] - a[0]) / length, (b[1] - a[1]) / length) if length else (1.0, 0.0)
    pts = [(a[0] + (b[0] - a[0]) * i / n, a[1] + (b[1] - a[1]) * i / n) for i in range(n + 1)]
    return pts, [(tx, ty)] * (n + 1), (length / n if length else 0.0)


def runs_of(flags):
    """[(i, j)] for maximal runs of True in flags."""
    out, i, n = [], 0, len(flags)
    while i < n:
        if flags[i]:
            j = i
            while j + 1 < n and flags[j + 1]:
                j += 1
            out.append((i, j))
            i = j + 1
        else:
            i += 1
    return out


class DSU:
    def __init__(self):
        self.p = {}

    def find(self, x):
        self.p.setdefault(x, x)
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def union(self, a, b):
        self.p[self.find(a)] = self.find(b)


# --------------------------------------------------------------------------
# Analysis

def analyse(args):
    board = pcbnew.LoadBoard(args.board)
    if not args.no_fill:
        pcbnew.ZONE_FILLER(board).Fill(board.Zones())

    # NetsByName() keys are wxString proxies; force str so lookups by name work
    netinfo = board.GetNetInfo()
    code_to_name = {ni.GetNetCode(): str(n) for n, ni in netinfo.NetsByName().items()}
    name_to_class = {str(n): str(ni.GetNetClassName()).split(",") for n, ni in netinfo.NetsByName().items()}

    class_file = args.classes or os.path.join(os.path.dirname(os.path.abspath(args.board)),
                                              "return_path_classes.json")
    class_overrides = load_class_overrides(class_file)

    # Board raster frame
    bb = board.GetBoardEdgesBoundingBox()
    if bb.GetWidth() == 0:
        bb = board.ComputeBoundingBox(False)
    pad = 1.0
    x0, y0 = bb.GetX() / NM - pad, bb.GetY() / NM - pad
    nx = int((bb.GetWidth() / NM + 2 * pad) / args.res) + 1
    ny = int((bb.GetHeight() / NM + 2 * pad) / args.res) + 1

    geo = dict(bbox=[round(bb.GetX() / NM, 3), round(bb.GetY() / NM, 3),
                     round(bb.GetWidth() / NM, 3), round(bb.GetHeight() / NM, 3)],
               outline=[], planes=defaultdict(list), tracks=[], vias=[], pads=[], refs=[])
    outline = pcbnew.SHAPE_POLY_SET()
    try:
        ok = board.GetBoardPolygonOutlines(outline, True)
    except TypeError:  # older API without aInferOutlineIfNecessary
        ok = board.GetBoardPolygonOutlines(outline)
    if ok:
        geo["outline"] = [rnd(c) for c in poly_contours(outline)]

    planes = {}
    zone_nets = set()
    for z in board.Zones():
        if z.GetIsRuleArea() or not z.IsFilled():
            continue
        for lid in z.GetLayerSet().Seq():
            if not pcbnew.IsCopperLayer(lid):
                continue
            fill = z.GetFilledPolysList(lid)
            contours = poly_contours(fill)
            if not contours:
                continue
            lname = board.GetLayerName(lid)
            if lname not in planes:
                planes[lname] = Plane((x0, y0), (ny, nx), args.res)
            planes[lname].add(contours, z.GetNetCode())
            zone_nets.add(z.GetNetname())
            geo["planes"][lname].append(dict(net=z.GetNetname(),
                                             c=[rnd(c) for c in poly_contours(unfractured(fill))]))

    stackup = parse_stackup(args.board)
    refs = reference_map(board, stackup, set(planes), args.ref)
    zpos = copper_z(board, stackup)

    # Net selection
    def selected(name):
        if not name or name in zone_nets or name.upper().startswith("GND"):
            return False
        if args.net and not any(fnmatch.fnmatchcase(name, g) for g in args.net):
            return False
        if args.netclass and not set(args.netclass) & set(name_to_class.get(name, [])):
            return False
        if not args.include_power and not args.net and not args.netclass:
            if POWER_RE.search(name.rsplit("/", 1)[-1]):
                return False
        return True

    tracks_by_net = defaultdict(list)
    vias_by_net = defaultdict(list)
    for t in board.GetTracks():
        name = t.GetNetname()
        if t.GetClass() == "PCB_VIA":
            vias_by_net[name].append(t)
            x, y = mm(t.GetPosition())
            geo["vias"].append([round(x, 3), round(y, 3), round(via_width(t), 3), name])
        else:
            tracks_by_net[name].append(t)
            pts, _, _ = sample_track(t, 0.25) if t.GetClass() == "PCB_ARC" else ([mm(t.GetStart()), mm(t.GetEnd())], 0, 0)
            geo["tracks"].append(dict(n=name, l=board.GetLayerName(t.GetLayer()),
                                      w=round(t.GetWidth() / NM, 3), p=rnd(pts)))

    pads_by_net = defaultdict(list)
    caps = []
    for fp in board.GetFootprints():
        ref = fp.GetReference()
        fp_pads = list(fp.Pads())
        side = "B" if fp.IsFlipped() else "F"
        fx, fy = mm(fp.GetPosition())
        geo["refs"].append([ref, round(fx, 2), round(fy, 2), side])
        for p in fp_pads:
            pads_by_net[p.GetNetname()].append(p)
            pb = p.GetBoundingBox()
            pside = "*" if p.GetAttribute() == pcbnew.PAD_ATTRIB_PTH else side
            geo["pads"].append([round(pb.GetX() / NM, 3), round(pb.GetY() / NM, 3),
                                round(pb.GetWidth() / NM, 3), round(pb.GetHeight() / NM, 3),
                                pside, p.GetNetname()])
        if re.match(r"^C\d", ref):
            caps.append((ref, [(p.GetNetname(), mm(p.GetPosition())) for p in fp_pads]))

    issues = []
    overlays = dict(voids=[], vias=[], edges=[], returns=[])
    summary = {}

    for net, tracks in sorted(tracks_by_net.items()):
        if not selected(net):
            continue
        cls = classify(net, class_overrides)
        s = dict(net=net, cls=cls,
                 netclass=",".join(c for c in name_to_class.get(net, []) if c != "Default") or "Default",
                 length=0.0, unref=0.0, edge=0.0, voids=0, void_extra=0.0,
                 transitions=0, via_extra=0.0, worst_stitch=None, broken=False,
                 area_ref=0.0, area_void=0.0, area_via=0.0,
                 sig_by_layer=defaultdict(float), under_by_plane=defaultdict(float))
        summary[net] = s

        def issue(kind, severity, layer, pos, value, msg):
            raw = severity
            if not args.no_class_severity and cls in CLASS_CAP:
                cap = CLASS_CAP[cls]
                if SEV_RANK[severity] < SEV_RANK[cap]:
                    severity = cap
            issues.append(dict(kind=kind, severity=severity, raw=raw, net=net, layer=layer,
                               x=None if pos is None else round(pos[0], 3),
                               y=None if pos is None else round(pos[1], 3),
                               value=None if value is None else round(value, 3), msg=msg))

        # ---- sample tracks against their reference plane
        runs = []
        edge_by_layer = defaultdict(float)
        for t in tracks:
            lname = board.GetLayerName(t.GetLayer())
            pts, tans, ds = sample_track(t, args.step)
            seg_len = ds * (len(pts) - 1)
            s["length"] += seg_len
            s["sig_by_layer"][lname] += seg_len
            ref = refs.get(lname)
            if not ref or ref[0] not in planes:
                s["unref"] += seg_len
                issue("NO_PLANE", "WARN", lname, pts[len(pts) // 2], seg_len,
                      f"no reference plane for {lname}")
                continue
            plane, h = planes[ref[0]], ref[1]
            labels = [plane.label(x, y) for x, y in pts]
            w = t.GetWidth() / NM
            off = w / 2 + args.edge_k * h
            hw = min(w / 2 + RETURN_BAND_K * h, RETURN_BAND_MAX)
            n = len(pts)

            # predicted return current: directly under each referenced stretch of trace
            for i, j in runs_of([lab != 0 for lab in labels]):
                if j > i:
                    s["under_by_plane"][ref[0]] += (j - i) * ds
                    overlays["returns"].append(dict(net=net, plane=ref[0], hw=round(hw, 3),
                                                    p=rnd(rdp(pts[i:j + 1], 0.005))))

            for i, j in runs_of([lab == 0 for lab in labels]):
                bounds, ends = [], []
                if i > 0:
                    bounds.append((pts[i - 1], labels[i - 1]))
                else:
                    ends.append((round(pts[0][0], 2), round(pts[0][1], 2), lname))
                if j < n - 1:
                    bounds.append((pts[j + 1], labels[j + 1]))
                else:
                    ends.append((round(pts[-1][0], 2), round(pts[-1][1], 2), lname))
                runs.append(dict(len=(j - i + 1) * ds, bounds=bounds, ends=ends, layer=lname,
                                 ref=ref[0], h=h, hw=hw, w=w, mid=pts[(i + j) // 2],
                                 pts=pts[max(i - 1, 0):j + 2]))
                s["unref"] += (j - i + 1) * ds

            edge_flags = [False] * n
            for i in range(n):
                if not labels[i]:
                    continue
                s["area_ref"] += ds * h
                if i > 0 and labels[i - 1] and labels[i - 1] != labels[i]:
                    issue("SPLIT", "FAIL", lname, pts[i], None,
                          f"reference on {ref[0]} changes {code_to_name.get(labels[i-1])} -> "
                          f"{code_to_name.get(labels[i])} under trace")
                if args.edge_k > 0:
                    x, y = pts[i]
                    tx, ty = tans[i]
                    if (plane.label(x - ty * off, y + tx * off) != labels[i] or
                            plane.label(x + ty * off, y - tx * off) != labels[i]):
                        edge_by_layer[lname] += ds
                        edge_flags[i] = True
            for i, j in runs_of(edge_flags):
                if j > i:
                    overlays["edges"].append(dict(net=net, layer=lname, p=rnd(pts[i:j + 1])))

        for lname, e in edge_by_layer.items():
            s["edge"] += e
            if e >= args.edge_report:
                issue("EDGE", "WARN", lname, None, e,
                      f"{e:.2f} mm within w/2+{args.edge_k:g}h of plane edge on {refs[lname][0]}")

        # ---- merge void runs that continue across segment joints
        dsu = DSU()
        by_end = defaultdict(list)
        for k, r in enumerate(runs):
            dsu.find(k)
            for e in r["ends"]:
                by_end[e].append(k)
        for ks in by_end.values():
            for k in ks[1:]:
                dsu.union(ks[0], k)
        groups = defaultdict(list)
        for k in range(len(runs)):
            groups[dsu.find(k)].append(runs[k])

        for g in groups.values():
            total = sum(r["len"] for r in g)
            bounds = [b for r in g for b in r["bounds"]]
            layer, refl, mid, h = g[0]["layer"], g[0]["ref"], g[0]["mid"], g[0]["h"]
            trace_pts = [rnd(r["pts"]) for r in g]
            ov = dict(net=net, layer=layer, ref=refl, hw=round(g[0]["hw"], 3), mid=rnd([mid])[0],
                      w=round(max(r["w"] for r in g), 3), trace=trace_pts, detours=[], extra=0.0,
                      kind="VOID", minor=False)
            if len(bounds) < 2:
                s["area_void"] += total * h  # lower bound: return path unknown
                if total > args.terminal_tol:
                    ov["kind"] = "TERMINAL"
                    overlays["voids"].append(ov)
                    issue("TERMINAL", "WARN", layer, mid, total,
                          f"{total:.2f} mm of trace at a pad/via end has no plane on {refl}")
                continue
            (p0, l0) = bounds[0]
            worst = 0.0
            for (p1, l1) in bounds[1:]:
                if l1 != l0:
                    ov["kind"] = "SPLIT"
                    issue("SPLIT", "FAIL", layer, mid, total,
                          f"trace crosses from {code_to_name.get(l0)} to {code_to_name.get(l1)} "
                          f"copper on {refl} over a {total:.2f} mm gap")
                    s["broken"] = True
                    continue
                found = geodesic(planes[refl], l0, p0, p1)
                if found is None:
                    ov["kind"] = "ISLAND"
                    issue("ISLAND", "FAIL", layer, mid, total,
                          f"{code_to_name.get(l0)} copper on {refl} is not connected across this "
                          f"{total:.2f} mm gap - no return path")
                    s["broken"] = True
                    continue
                geo_len, path = found
                extra = max(geo_len - (total + args.step), 0.0)
                worst = max(worst, extra)
                ov["detours"].append(rnd(path))
                # loop = region between the trace over the gap and the detour
                if len(g) == 1:
                    s["area_void"] += shoelace(list(g[0]["pts"]) + path[::-1])
                else:
                    s["area_void"] += shoelace(path)
            s["voids"] += 1
            s["void_extra"] += worst
            ov["extra"] = round(worst, 3)
            # Minor detours (round an antipad, say) are only drawn as part of a
            # selected net's return path, so the board view stays readable.
            ov["minor"] = ov["kind"] == "VOID" and worst < args.report_min
            overlays["voids"].append(ov)
            if worst >= args.report_min:
                sev = "FAIL" if worst > args.max_detour else "WARN"
                issue("VOID", sev, layer, mid, worst,
                      f"{total:.2f} mm over a void in {refl}: return detours +{worst:.2f} mm")

        # ---- layer transitions (vias and PTH pads of this net)
        transitions = [(v, mm(v.GetPosition()), via_width(v) / 2) for v in vias_by_net.get(net, [])]
        transitions += [(p, mm(p.GetPosition()), max(p.GetSize().x, p.GetSize().y) / NM / 2)
                        for p in pads_by_net.get(net, []) if p.GetAttribute() == pcbnew.PAD_ATTRIB_PTH]
        for obj, pos, rad in transitions:
            layers = set()
            for t in tracks:
                for e in (mm(t.GetStart()), mm(t.GetEnd())):
                    if math.dist(e, pos) <= rad + 0.01:
                        layers.add(board.GetLayerName(t.GetLayer()))
            for p in pads_by_net.get(net, []):
                if p.GetAttribute() == pcbnew.PAD_ATTRIB_SMD and p.HitTest(obj.GetPosition()):
                    layers.update(board.GetLayerName(l) for l in p.GetLayerSet().CuStack())
            ref_layers = sorted({refs[l][0] for l in layers if l in refs and refs[l][0] in planes},
                                key=lambda n: board.GetLayerID(n))
            if len(ref_layers) < 2:
                continue
            s["transitions"] += 1
            ra, rb = ref_layers[0], ref_layers[-1]
            dz = abs(zpos.get(rb, 0) - zpos.get(ra, 0))
            na = planes[ra].label_near(*pos, 1.5)
            nb = planes[rb].label_near(*pos, 1.5)
            where = f"{'/'.join(sorted(layers))} via"
            ov = dict(net=net, pos=rnd([pos])[0], layers=sorted(layers), planes=[ra, rb],
                      to=None, d=None, kind="stitch")
            overlays["vias"].append(ov)
            if not na or not nb:
                ov["sev"] = "FAIL"
                issue("VIA", "FAIL", where, pos, None,
                      f"no plane copper near via on {ra if not na else rb}")
                s["broken"] = True
                continue
            if na == nb:
                gnd = code_to_name[na]
                cands = [mm(v.GetPosition()) for v in vias_by_net.get(gnd, [])
                         if v.IsOnLayer(board.GetLayerID(ra)) and v.IsOnLayer(board.GetLayerID(rb))]
                cands += [mm(p.GetPosition()) for p in pads_by_net.get(gnd, [])
                          if p.GetAttribute() == pcbnew.PAD_ATTRIB_PTH]
                if not cands:
                    ov["sev"] = "FAIL"
                    issue("VIA", "FAIL", where, pos, None, f"no {gnd} via connects {ra} and {rb}")
                    s["broken"] = True
                    continue
                to = min(cands, key=lambda c: math.dist(pos, c))
                d = math.dist(pos, to)
                s["via_extra"] += 2 * d
                s["area_via"] += d * dz
                s["worst_stitch"] = max(s["worst_stitch"] or 0, d)
                sev = "FAIL" if d > args.max_stitch else "INFO"
                ov.update(to=rnd([to])[0], d=round(d, 3), sev=sev, via_net=gnd)
                issue("VIA", sev, where, pos, d,
                      f"{ra}->{rb} ({gnd}): nearest stitching via {d:.2f} mm, return +{2*d:.2f} mm")
            else:
                a_name, b_name = code_to_name[na], code_to_name[nb]
                best = None
                for cref, cpads in caps:
                    pa = [q for n_, q in cpads if n_ == a_name]
                    pb = [q for n_, q in cpads if n_ == b_name]
                    if pa and pb:
                        for x in pa:
                            for y in pb:
                                loop = math.dist(pos, x) + math.dist(x, y) + math.dist(y, pos)
                                if best is None or loop < best[0]:
                                    best = (loop, cref, x, y)
                if best is None:
                    ov["sev"] = "FAIL"
                    issue("VIA", "FAIL", where, pos, None,
                          f"{ra}({a_name})->{rb}({b_name}): no capacitor bridges the planes")
                    s["broken"] = True
                    continue
                loop, cref, cpos, cpos_b = best
                s["via_extra"] += loop
                s["area_via"] += loop / 2 * dz
                s["worst_stitch"] = max(s["worst_stitch"] or 0, loop / 2)
                sev = "FAIL" if loop / 2 > args.max_cap else "WARN"
                ov.update(to=rnd([cpos])[0], pads=rnd([cpos, cpos_b]), d=round(loop / 2, 3), sev=sev,
                          kind="cap", cap=cref)
                issue("VIA", sev, where, pos, loop,
                      f"{ra}({a_name})->{rb}({b_name}): return via {cref}, +{loop:.2f} mm")

    for s in summary.values():
        s["return_est"] = None if s["broken"] else s["length"] + s["void_extra"] + s["via_extra"]
        s["ratio"] = None if s["return_est"] is None or not s["length"] else s["return_est"] / s["length"]
        s["loop_area"] = s["area_ref"] + s["area_void"] + s["area_via"]
        sev = [i["severity"] for i in issues if i["net"] == s["net"]]
        s["status"] = "FAIL" if "FAIL" in sev else "WARN" if "WARN" in sev else "OK"
        for k in ("sig_by_layer", "under_by_plane"):
            s[k] = {name: round(v, 3) for name, v in s[k].items()}
        for k, v in list(s.items()):
            if isinstance(v, float):
                s[k] = round(v, 3)
    for ov in overlays["vias"]:
        cls = summary[ov["net"]]["cls"]
        if ov.get("sev") and not args.no_class_severity and cls in CLASS_CAP:
            if SEV_RANK[ov["sev"]] < SEV_RANK[CLASS_CAP[cls]]:
                ov["sev"] = CLASS_CAP[cls]

    geo["planes"] = dict(geo["planes"])
    meta = dict(board=os.path.basename(args.board), generated=datetime.datetime.now().isoformat(timespec="minutes"),
                classes=CLASSES, class_file=class_file if class_overrides else None,
                thresholds=dict(max_detour=args.max_detour, max_stitch=args.max_stitch,
                                max_cap=args.max_cap, edge_k=args.edge_k, report_min=args.report_min),
                references={k: dict(layer=v[0], h=round(v[1], 4)) for k, v in refs.items()},
                planes=sorted(planes, key=lambda l: zpos.get(l, 0)), return_band_k=RETURN_BAND_K,
                board_mtime=file_time(args.board), live=False)
    return dict(meta=meta, nets=list(summary.values()), issues=issues, overlays=overlays, geo=geo)


def file_time(path):
    """A file's modification time as an ISO string, to tell when the board was saved after a run."""
    return datetime.datetime.fromtimestamp(os.path.getmtime(path)).isoformat(timespec="seconds")


# --------------------------------------------------------------------------
# Output

def sort_nets(nets):
    order = {"FAIL": 0, "WARN": 1, "OK": 2}
    return sorted(nets, key=lambda s: (order[s["status"]], CLASSES.index(s["cls"]), -(s["ratio"] or 99)))


def report(result, args):
    print("Reference planes:")
    for sig, r in sorted(result["meta"]["references"].items()):
        print(f"  {sig:8s} -> {r['layer']:8s} h={r['h']:.3f} mm")
    print()

    hdr = (f"{'status':6} {'net':30} {'class':12} {'trace':>7} {'return':>7} {'ratio':>6} "
           f"{'loop mm2':>8} {'voids':>5} {'vias':>4} {'stitch':>6}")
    print(hdr)
    print("-" * len(hdr))
    for s in sort_nets(result["nets"]):
        ret = "open" if s["return_est"] is None else f"{s['return_est']:.2f}"
        ratio = "-" if s["ratio"] is None else f"{s['ratio']:.2f}"
        st = "-" if s["worst_stitch"] is None else f"{s['worst_stitch']:.2f}"
        print(f"{s['status']:6} {s['net'][-30:]:30} {s['cls']:12} {s['length']:7.2f} {ret:>7} {ratio:>6} "
              f"{s['loop_area']:8.2f} {s['voids']:5d} {s['transitions']:4d} {st:>6}")
    print("\n(lengths mm; return = trace + void detours + via transfer; loop = estimated current "
          "loop area; stitch = worst distance to return via/cap)\n")

    shown = [i for i in result["issues"] if i["severity"] != "INFO" or args.verbose]
    shown.sort(key=lambda i: (SEV_RANK[i["severity"]], i["net"]))
    if shown:
        print("Issues:")
        for i in shown:
            loc = "" if i["x"] is None else f" @ ({i['x']:.2f}, {i['y']:.2f})"
            capped = f" (capped from {i['raw']})" if i["raw"] != i["severity"] else ""
            print(f"  {i['severity']:4} {i['kind']:8} {i['net']} [{i['layer']}]{loc}: {i['msg']}{capped}")
    nfail, line = summary_line(result)
    print(f"\n{line}")
    return nfail


def summary_line(result):
    nets = result["nets"]
    nf = sum(1 for s in nets if s["status"] == "FAIL")
    nw = sum(1 for s in nets if s["status"] == "WARN")
    return nf, f"{len(nets)} nets checked: {nf} FAIL, {nw} WARN"


def render_html(result):
    template = os.path.join(os.path.dirname(os.path.abspath(__file__)), "return_path_report.html")
    with open(template, encoding="utf-8") as f:
        html = f.read()
    data = json.dumps(result, separators=(",", ":")).replace("</", "<\\/")
    return html.replace("/*__DATA__*/null", data)


def write_html(result, path):
    with open(path, "w", encoding="utf-8") as f:
        f.write(render_html(result))


def write_json(result, path):
    with open(path, "w", encoding="utf-8") as f:
        json.dump({k: result[k] for k in ("meta", "nets", "issues")}, f, indent=2)


# --------------------------------------------------------------------------
# Live report

def serve(args):
    """Serve the HTML report on localhost; its Rerun button re-analyses the saved board in place.

    HTTP runs on worker threads, because a browser can hold idle connections open, while every
    analysis runs here on the main thread, one at a time, since pcbnew isn't thread-safe."""
    import http.server
    import queue
    import threading

    live = {}
    jobs = queue.Queue()

    def run():
        result = analyse(args)
        result["meta"]["live"] = True
        live["html"] = render_html(result).encode("utf-8")
        live["generated"] = result["meta"]["generated"]
        if args.json:
            write_json(result, args.json)
        if args.html:
            write_html(result, args.html)
        line = summary_line(result)[1]
        print(f"[{datetime.datetime.now():%H:%M:%S}] {line}", flush=True)
        return line

    class Handler(http.server.BaseHTTPRequestHandler):
        def send(self, code, body, ctype="application/json"):
            if isinstance(body, dict):
                body = json.dumps(body).encode("utf-8")
            self.send_response(code)
            self.send_header("Content-Type", ctype)
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Cache-Control", "no-store")
            self.end_headers()
            self.wfile.write(body)

        def local(self):
            # Only answer requests addressed to this machine, which guards against DNS rebinding.
            host = (self.headers.get("Host") or "").rsplit(":", 1)[0]
            if host in ("127.0.0.1", "localhost", "[::1]"):
                return True
            self.send(403, b"forbidden", "text/plain")
            return False

        def do_GET(self):
            if not self.local():
                return
            path = self.path.split("?", 1)[0]
            if path in ("/", "/index.html"):
                self.send(200, live["html"], "text/html; charset=utf-8")
            elif path == "/api/status":
                self.send(200, dict(board_mtime=file_time(args.board), generated=live["generated"]))
            else:
                self.send(404, b"not found", "text/plain")

        def do_POST(self):
            if not self.local():
                return
            # The custom header makes a cross-site POST need a CORS preflight, which this server refuses.
            if self.path != "/api/rerun" or self.headers.get("X-Return-Path") != "rerun":
                self.send(404, b"not found", "text/plain")
                return
            job = dict(done=threading.Event())
            jobs.put(job)
            job["done"].wait()
            if "error" in job:
                self.send(200, dict(ok=False, error=job["error"]))
            else:
                self.send(200, dict(ok=True, summary=job["line"]))

        def log_message(self, fmt, *a):  # keep the console to one line per run
            pass

    run()
    server = http.server.ThreadingHTTPServer(("127.0.0.1", args.serve), Handler)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    print(f"Live report on http://127.0.0.1:{args.serve}/  (Ctrl+C to stop)", flush=True)
    try:
        while True:
            try:
                job = jobs.get(timeout=0.5)  # time out now and then so Ctrl+C gets through on Windows
            except queue.Empty:
                continue
            try:
                job["line"] = run()
            except Exception as ex:  # a half-written board, say: keep serving the last good report
                job["error"] = f"{type(ex).__name__}: {ex}"
            finally:
                job["done"].set()
    except KeyboardInterrupt:
        pass


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("board", help=".kicad_pcb file")
    ap.add_argument("--net", action="append", default=[], help="net name glob (repeatable)")
    ap.add_argument("--netclass", action="append", default=[], help="net class (repeatable)")
    ap.add_argument("--include-power", action="store_true", help="also check power-looking nets")
    ap.add_argument("--ref", action="append", default=[], metavar="SIG=REF",
                    help="override reference layer, e.g. In2.Cu=In1.Cu")
    ap.add_argument("--classes", help="JSON of net glob -> class (default: return_path_classes.json by the board)")
    ap.add_argument("--no-class-severity", action="store_true",
                    help="don't cap issue severity on low-speed/static nets")
    ap.add_argument("--step", type=float, default=0.1, help="track sampling step, mm (0.1)")
    ap.add_argument("--res", type=float, default=0.05, help="plane raster resolution, mm (0.05)")
    ap.add_argument("--max-detour", type=float, default=1.0, help="void detour FAIL limit, mm (1.0)")
    ap.add_argument("--report-min", type=float, default=0.2, help="list voids with detour >= this, mm (0.2)")
    ap.add_argument("--max-stitch", type=float, default=2.0, help="signal-via to return-via FAIL limit, mm (2.0)")
    ap.add_argument("--max-cap", type=float, default=5.0, help="signal-via to bridging cap FAIL limit, mm (5.0)")
    ap.add_argument("--terminal-tol", type=float, default=1.0,
                    help="ignore unreferenced trace ends shorter than this (antipads), mm (1.0)")
    ap.add_argument("--edge-k", type=float, default=3.0, help="edge rule multiple of h; 0 disables (3)")
    ap.add_argument("--edge-report", type=float, default=1.0, help="list edge-hugging >= this per net/layer, mm")
    ap.add_argument("--no-fill", action="store_true", help="use zone fills saved in the file")
    ap.add_argument("--json", help="write summary + issues as JSON")
    ap.add_argument("--html", help="write an interactive HTML report (board view + per-net table)")
    ap.add_argument("--serve", nargs="?", type=int, const=8765, metavar="PORT",
                    help="serve a live HTML report on http://127.0.0.1:PORT (8765) with a Rerun button")
    ap.add_argument("--fail", action="store_true", help="exit 1 if any net FAILs")
    ap.add_argument("-v", "--verbose", action="store_true", help="also list INFO items (every via)")
    args = ap.parse_args()

    if args.serve is not None:
        serve(args)
        return
    result = analyse(args)
    nfail = report(result, args)
    if args.json:
        write_json(result, args.json)
    if args.html:
        write_html(result, args.html)
        print(f"HTML report: {args.html}")
    sys.exit(1 if args.fail and nfail else 0)


if __name__ == "__main__":
    main()
