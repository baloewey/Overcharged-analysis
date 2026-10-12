# Draft: oval (Cassini) sweep for the server — NOT launched

Prepared 2026-10-03 for the user's approval (decision of 2026-10-02: prepare, don't launch).

## Design (matched variables, user's idea 2)
- **Shape:** b/a = 3.0, 2.5, 2.25, 2.0, 1.8, 1.6, 1.45, 1.3, 1.1 (aspect 1.1–2.5; dense through the crossover seen at
  aspect ≈ 1.15–1.5 in the pilots; b/a ≥ 1.45 convex).
- **Size:** equivalent radius √(A/π) = 1.0, 1.5, 2.0 ℓ_a (inside the region where circles keep the +1 layer, 0.75–2.0 ℓ_a).
- **Linear defect density:** ρ = 2k/P ≈ 0.04, 0.06, 0.09, 0.12 (k rounded; ρ ≤ 0.15 resolution limit).
- 3 seeds; ζ = 0.07; eMode 3; 10k + 100k steps; defects/20, director + flow/1000 (as the circle sweep).
- **324 runs, 2.7 × 10¹² particle-steps ≈ 1.2 days at 96-way**; largest run 124k particles.
- Generator: `runs/cassini/gen_oval_sweep.py OUTROOT` (writes inputs + `jobs_oval_sweep.txt`, largest first).
  Manifest (dry run): `reports/oval_sweep_draft_manifest.txt`.

## Ready on the server
- Branch `cassini` binary built and validated in `~/Overcharged/niamh-overcharged-cassini/` (e63155b); production
  binary untouched. To run: copy `run_queue.sh` with EXE pointing to the new binary (as `runs/cassini/run_queue_cassini.sh`).
- 366 GB free.

## Open choices for the user
1. Add **mirror runs** (−k) on a subset (e.g. 3 shapes × 2 sizes × 2 ρ)? The circle mirror results show the bulk +1 only
   builds up for |k| ≳ 13, so the oval mirror question is "does elongation also block it?".
2. Add a **high-k / high-ρ** column (ρ ≈ 0.14) given that a crowded wall partly restores the layer in pinched ovals?
3. Add **ellipses** as a second family on the server, or keep them on the office machine (weekend pilot running)?
4. Longer runs for the slow cases (vortex reversals at small aspect)?
