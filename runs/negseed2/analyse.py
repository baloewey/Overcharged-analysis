#!/usr/bin/env python3
"""negseed2 (pair-creation-cost test): wall -1/2 count (< 4 cells) and interior charge (> 4 cells) vs time, per run; means over 1000-t.u. windows and the last 3000 t.u.
Compare 'base' vs 'full' starts with the sweep value (shape_runs_mirror / mirror_tilt). Writes summary.txt and counts.png."""
import glob, json, os, sys, numpy as np
sys.path.insert(0, "/home/baloewe/Physics/Overcharged/runs/cavity")
from edge_criteria import read_with_angle  # noqa
from wall_fourier import wall_polyline  # noqa
from scipy.spatial import cKDTree
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
SWEEP = {"sparse_circle": 0.27, "sparse_ellipse": -0.41, "small_circle": 1.32, "crowded_circle": 2.03}
out = ["negseed2: extra wall -1/2 (N - 2|k|) and interior charge, windows of 1000 t.u.; last = 3000-6000 t.u.; sweep = mirror sweep value (s01-s03)"]
fig, axs = plt.subplots(1, 4, figsize=(14, 3.4), sharey=False)
for j, case in enumerate(SWEEP):
    ax = axs[j]
    for d in sorted(glob.glob(f"{case}_*_s*")):
        if not os.path.exists(d + "/DONE"): continue
        bc = json.load(open(d + "/input.json"))["BC"][0]; k = bc["customAnch"]; W, s, P = wall_polyline(bc); tree = cKDTree(W)
        t, fr = read_with_angle(d, np.array(bc["Q"][:2], float), 0.0); t = np.array(t)
        n = np.array([np.sum((tree.query(f[:, :2])[0] < 4) & (f[:, 2] < 0)) if len(f) else 0 for f in fr]) - 2 * abs(k)
        qi = np.array([0.5 * np.sum(np.sign(f[tree.query(f[:, :2])[0] >= 4, 2])) if len(f) else 0 for f in fr])
        win = [f"{n[(t >= a) & (t < a + 1000)].mean():5.2f}" for a in range(0, 6000, 1000)]; last = t >= 3000
        out.append(f"{d:28s} k {k}: windows {' '.join(win)} | last extra {n[last].mean():5.2f} interior {qi[last].mean():5.2f} | sweep {SWEEP[case]}")
        sm = np.convolve(n, np.ones(25) / 25, mode="same")
        ax.plot(t, sm, color="C3" if "full" in d else "C0", lw=0.8, label=("Full start" if "full" in d else "Base start") if d.endswith("s1") else None)
    ax.axhline(SWEEP[case], color="k", ls=":", lw=1, label="Sweep value"); ax.set_title(case.replace("_", " ").capitalize()); ax.set_xlabel("Time (t.u.)")
    if j == 0: ax.set_ylabel("Extra wall $-1/2$, $N-2|k|$"); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig("counts.png", dpi=150); fig.savefig("counts.pdf")
open("summary.txt", "w").write("\n".join(out) + "\n"); print("\n".join(out))
