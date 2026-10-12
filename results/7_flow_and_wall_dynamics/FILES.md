# 7_flow_and_wall_dynamics

Mechanism: circulation, wall-defect tilt, wall lattice (pinned / vacant), exchange with the bulk.

Copied by `tools/build_results.py`; regenerate rather than edit. Files marked KEY are the main figures; SUPERSEDED files are kept for the record only.

| File | What it is | Copied from |
|---|---|---|
| `wall_dynamics_ovals_vs_circles.png` | KEY: pinning, upstream drift, drift vs vacancies, k scan | `runs/cassini/fig_wall_dynamics.png` |
| `wall_dynamics_ovals_vs_circles.pdf` | Same, PDF | `runs/cassini/fig_wall_dynamics.pdf` |
| `wall_dynamics_ovals.csv` | Per oval/ellipse run | `runs/cassini/wall_fourier_dyn.csv` |
| `wall_dynamics_circles_k06.csv` | Per circle run, k = 6 (same script) | `runs/cassini/wall_fourier_dyn_circles_k06.csv` |
| `wall_fourier_circles_static.png` | KEY: circle wall lattice: period, occupancy, order, pinning (+k and -k) | `runs/cavity/analysis_fourier/fourier_p1_static.png` |
| `wall_fourier_circles_dynamic.png` | Lattice memory, rearrangement time, sliding speed (none) | `runs/cavity/analysis_fourier/fourier_p2_dynamic.png` |
| `wall_fourier_circles.pdf` | Both pages, PDF | `runs/cavity/analysis_fourier/fourier_maps.pdf` |
| `wall_fourier_circles.csv` | Per circle run | `runs/cavity/analysis_fourier/wall_fourier_all.csv` |
| `flow_mean_circulation.png` | KEY: mean flow and circulation collapse on R/l_a; breakdown above 2.5 l_a | `runs/cavity/analysis_flow/flow_p1_mean_flow.png` |
| `flow_enstrophy.png` | Enstrophy vs R/l_a | `runs/cavity/analysis_flow/flow_p2_enstrophy.png` |
| `wall_defect_tilt.png` | Wall +1/2 lean downstream 60-65 deg (30-40 at high rho) | `runs/cavity/analysis_flow/flow_p3_tilt.png` |
| `flow.pdf` | All flow pages, PDF | `runs/cavity/analysis_flow/flow_maps.pdf` |
| `flow_tilt.csv` | Per run | `runs/cavity/analysis_flow/flow_tilt.csv` |
| `exchange_rates.txt` | Wall -> bulk escape rates: 10-20x faster in elongated ovals | `runs/cavity/analysis_weekend/exchange_rates.txt` |
| `wall_defect_labels_circles_ovals.png` | Where the wall defects go: circles keep them, ovals lose them in 100-200 time units | `runs/cavity/analysis_spatial/bc_labels_k06.png` |
| `wall_defect_labels_circles_ovals.txt` | Same, numbers | `runs/cavity/analysis_spatial/bc_labels_k06.txt` |
| `wall_defect_labels_small_ovals.png` | Same for the small ovals | `runs/night_2026-10-01/bc_labels_P3.png` |
| `wall_defect_labels_lobe_oval.png` | Same for the lobe-shaped oval | `runs/cassini/P4_lobe_ref/bc_labels.png` |
