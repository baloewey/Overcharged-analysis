# Candidate SI movies (paper 1), 2026-10-15
| File | Content | Source run | Script |
|---|---|---|---|
| S1_circle_overcharged.mp4 | Active circle, k = 6, area 2500: 10 wall +½, bulk +1, single vortex, tilted layer; wall count and bulk charge traces (t 1500–2100) | runs/fineflow/eps1.0_s1 | master/theory/channel/cavity_movie.py |
| S2_passive_swirl_ejection.mp4 | Passive circle under an imposed vortex (A 0.01, λ 2): the ring of 12 ejects two (t 0–600) | runs/swirl_threshold/passive_A0.01_s2 | cavity_movie.py |
| S3_ellipse_leaking.mp4 | Active aspect-2 ellipse: two vortices, layer leaks, bulk +2–3 (t 1500–2100) | runs/fineflow/eps2.0_s1 | cavity_movie.py |
| S4_escape_event.mp4 | One depleting escape, zoomed, with vortex-lobe overlay | runs/fineflow (eps 2.0 s2, event 4) | master/theory/energetics/movie/escape_movie.py |
| S5_channel_flat_limit.mp4 | Winding channel (W 40, +6/+6), flat limit of the wall law | runs/channel/W40_m6_L192_s1 | channel_movie.py |
| X_topological_current_W24.mp4 | Mixed channel (+6/−6, W 24) with accumulated charge transferred vs a symmetric control; for paper 2/3, not paper 1 | runs/channel_mixedW/W24_Qb+6_Qt-6_act_s1 | current_movie.py |
Not tracked in git (mp4); regenerate with the scripts.
