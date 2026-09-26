# KiCad Tools

Scripts for checking and analysing KiCad boards, shared across my projects.

Everything runs under **KiCad's bundled Python**, which already provides `pcbnew` and `numpy`, so
there is nothing to install. Tested with KiCad 10.0.1 on Windows.

## Layout

```
return_path/    return-path checker: plane voids, via stitching, loop area, HTML report
```

## Running

KiCad 10's Python on Windows lives at `%LOCALAPPDATA%\Programs\KiCad\10.0\bin\python.exe`:

```bash
"$LOCALAPPDATA/Programs/KiCad/10.0/bin/python.exe" return_path/return_path_check.py path/to/board.kicad_pcb --html return_path.html
```

In Git Bash, `export MSYS_NO_PATHCONV=1` before passing a net glob such as `--net '/Micro/*'`,
otherwise Bash rewrites it into a Windows path and nothing matches.

Scripts read the board **from disk** — save in pcbnew first. They never write the board back.

## return_path — return-path checker

For every signal net, estimates the path its return current takes through the reference plane and
compares it with the trace. Reference planes come from the board stackup: each copper layer is
referenced to the nearest copper layer (by dielectric thickness) that carries zone fill. Zones are
refilled in memory before analysis, so stale or unfilled zones in the file don't matter.

### What it flags

| Kind | Meaning |
|---|---|
| `VOID` | Trace crosses a hole or gap in its plane. The detour is the shortest path through plane copper between the entry and exit points (A* on a 0.05 mm raster of the fill), minus the trace length over the gap. FAIL above `--max-detour`. |
| `ISLAND` | Copper either side of the gap isn't connected on that layer at all — no local return path. |
| `SPLIT` | The plane under the trace changes net (e.g. GND pour to +3V3 pour); return must go through a capacitor. |
| `TERMINAL` | A trace end (pad/via) sits over a void longer than `--terminal-tol`. |
| `VIA` | Signal changes layers and the reference plane changes. Same-net planes: distance to the nearest stitching via, return ≈ +2 × that. Different-net planes: nearest capacitor bridging them. FAIL above `--max-stitch` / `--max-cap`. |
| `EDGE` | Trace is referenced, but the plane edge is within w/2 + k·h of its centreline (k = `--edge-k`, h = dielectric height), so the return current is crowded. |

### Per-net figures

- **Return length** — trace length + void detours + via transfers.
- **Ratio** — return ÷ trace; 1.00 means the return runs directly under the trace the whole way.
- **Loop area** — trace length × h, plus the plan-view area between each void crossing and its
  detour, plus via-to-return-via distance × plane-to-plane spacing. It's an estimate for comparing
  nets and finding the worst ones, not a field solve.

### Net classes

Nets are classified from their names as `rf`, `clock`, `high-speed`, `analog`, `medium`,
`low-speed` or `static`. Issues on low-speed nets are capped at WARN and on static nets at INFO,
so the FAIL list is the nets where return path actually matters (`--no-class-severity` turns this
off).

Names don't always say what a net does (a `GPIO_0` that is really the codec's MCLK, an unnamed
coupling-cap net), so put a `return_path_classes.json` next to the board to override by net-name
glob. First match wins; see [`return_path/classes.example.json`](return_path/classes.example.json).

### Outputs

- Text summary and issue list on stdout.
- `--html FILE` — self-contained interactive report: each signal layer drawn over its reference
  plane, traces coloured by status or class, void crossings with their detour paths, each signal via
  linked to its nearest return via, and a sortable, filterable per-net table. The page template is
  `return_path/return_path_report.html`.

  Click a trace or a row to select a net and see its **signal path** (solid) and **predicted return
  path**: a translucent, hatched band with a dashed centreline, so it can't be mistaken for copper.
  The band runs in the plane directly under the trace (±3h wide, where about 80 % of the return
  current flows), detours around plane gaps, and at each layer change goes out to the nearest return
  via on one plane and back on the other. Each reference plane has its own colour, and a trace is
  drawn in the colour of the plane carrying its return, so you can see where the return changes
  plane.
- `--json FILE` — summary and issues for scripting.
- `--fail` — exit 1 if any net FAILs, for CI.

### Common options

| Option | Default | |
|---|---|---|
| `--net GLOB` / `--netclass NAME` | all signal nets | Restrict the check (repeatable). |
| `--include-power` | off | Also check nets that look like power rails. |
| `--ref SIG=REF` | from stackup | Override a layer's reference plane, e.g. `In2.Cu=In1.Cu`. |
| `--max-detour` | 1.0 mm | Void detour FAIL limit. |
| `--max-stitch` | 2.0 mm | Signal via to return via FAIL limit. |
| `--max-cap` | 5.0 mm | Signal via to bridging cap FAIL limit. |
| `--edge-k` | 3 | Edge rule multiple of h; 0 disables. |
| `--step` / `--res` | 0.1 / 0.05 mm | Trace sampling step / plane raster resolution. |
| `--no-fill` | off | Use the zone fills saved in the file instead of refilling. |

`--help` lists the rest.

### Limitations

- Models plane copper only. Return current carried by neighbouring traces isn't counted.
- Plan-view geometry: the current spread under a trace and skin/proximity effects aren't modelled.
  Detour lengths are raster shortest paths, within about 3 %.
- Capacitors are recognised by reference designator (`C1`, `C2`, …).
- KiCad's SWIG bindings don't expose the stackup, so it is parsed from the `.kicad_pcb` file.
