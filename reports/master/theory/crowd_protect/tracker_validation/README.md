# Defect-tracker validation (2026-10-22)

**Why.** The channel "topological current" turned out to be a tracking artefact, and cavity escape counts were undercounted, most of all in crowded layers. The cause in both cases was nearest-neighbour linking without motion prediction, which loses fast defects.

## Linkers compared (`trackers.py`)

| Name | Rule | Used in |
|---|---|---|
| **CHAIN** | 10 t.u. frames; each start followed to the nearest same-sign defect within 6 cells; no prediction | `shape_runs.csv` escape, `events.npz`, `kinetics_carriers` |
| **NN** | 2 t.u. frames; mutual nearest neighbour within 3.5 cells of the last position; up to 2 missed frames; no prediction | `exchange_count.py`, `wall_budget.py`, `kink_life.py`, `depth_memory.py` |
| **NEW v1** | Position predicted from the last displacement; global Hungarian assignment; 9-cell gate; up to 2 missed frames | `exchange_count2.py`, `shape_esc2.py`, `escape_events2.py`, `channel/crossing2/tracker2.py` (evening of 22 Oct) |
| **NEW v2** | As v1, but the cost is the smaller of the distances to the predicted and to the last position | Production from the late evening of 22 Oct |

## Tests and results

### 1. Synthetic ground truth (`synth.py` → `synth.txt`)

The simulation is a circle with P = 200:
- a repelling wall layer, sparse (N = 8) or crowded (N = 28);
- bulk +½ with the measured step distribution and persistent directions;
- escapes launched inward, and entries into wall gaps;
- detections rounded to half-cells, with 2% missed.

Results:

| Linker | Recovered / true escape probability |
|---|---|
| NEW v2 | 0.96 (sparse), 0.99 (crowded) |
| NN | 0.61, 0.63 |
| CHAIN | 0.11, 0.04 |

Identity switches per 1000 links:

| Linker | Bulk, sparse | Bulk, crowded | Wall, sparse | Wall, crowded |
|---|---|---|---|---|
| NEW v2 | 10.7 | 2.8 | 2.2 | 0.9 |
| NN | 16.3 | 5.3 | 1.6 | 0.5 |

**NEW v1** recovers escapes equally well (0.96 and 0.98). But in crowded walls it makes 237 wall switches per 1000 links: a single swap gets extrapolated by the prediction and cascades along the ring. v2 removes this.

### 2. Exact closure on real data (`real_tests.py` → `real_tests.txt`)

The charge enclosed by a closed curve changes only by crossings, because pair events are neutral. The curves are:
- d = 4 and d = 8 cells from the wall, in circles and aspect-2 ellipses at ρ = 0.04, 0.09 and 0.14;
- the mid-line of the W24 mixed and symmetric channels.

| Linker | Tracked minus exact (rms) |
|---|---|
| NEW (v1 and v2) | ≤ 1.5 charges in every case |
| NN | up to 8–10 (for example 8.4 at d = 8 in a crowded circle; 9.6 at the mid-line of the mixed channel) |

### 3. Frame decimation (same file)

The escape probability is recomputed using every second frame (4 t.u.).

| Linker | Ratio (half frame rate / full) |
|---|---|
| NEW | 0.94–1.06 |
| NN | 0–0.59 |

### 4. Gate (same file)

NEW gives the same result with 6- and 9-cell gates. A 12-cell gate doubles the escapes in a crowded circle (false links), so 9 is the upper limit.

### 5. Checking movies (`validation_movie.py` → `movies/val_*.mp4`)

OLD (NN) and NEW v2 run side by side on the same frames:
- each track in its own colour, with its trajectory since birth;
- green ring = track born, black cross = track ended;
- red segment = a link longer than 4 cells;
- orange star = an escape detected by that linker.

Each movie is 250 frames (500 t.u.) at 6 fps. The window opens 40 frames before the first escape that NEW finds and OLD misses. Runs:
- sparse and crowded aspect-2 ellipses;
- a crowded circle;
- the W24 mixed channel.

## Status of results built on the older linkers

- **Already corrected:** the escape rates, the per-event outcomes and the exchange hazards.
- **Probably adequate:** analyses that follow slow wall defects (kink statistics, `depth_memory` windows); still to recheck with NEW v2.
- **To recheck:** the wall budget (NN linking; its sign flip), and anything in `kinetics_carriers` / `energetics_general` that uses CHAIN escapes (for example, the escaper energetics sample is biased towards slow escapes in crowded layers).
