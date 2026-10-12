#!/usr/bin/env python3
"""Is the wall lattice pinned by the anchoring pattern (physical) or by a grid (artefact)?
The anchoring director psi = k phi (+ const) is invariant only under rotations by 2 pi/(2k-2). With K1 != K3 the energy of the wall lattice
depends on its angular position relative to this pattern, with period 2 pi/(2k-2): pinning sites fixed in the lab frame, the same for every seed
and every radius. Grid pinning (period ~1 cell) would put the lattice at a random phase in each run.
For each run: wall +1/2 (within BAND cells of the wall) angular positions phi_j; Z = < (1/N_w) sum_j exp(i (2k-2) phi_j) >_frames.
|Z| ~ 1: lattice of 2k-2 static sites; arg Z / (2k-2) = lattice phase (in rad, defined mod 2 pi/(2k-2)).
Across seeds and densities: mean resultant length Rbar of exp(i arg Z) (1 = identical phase in every run; random phases: ~1/sqrt(n)).
Writes phase_test.txt."""
import glob, json, os, sys
import numpy as np
ROOT = "/home/baloewe/Physics/Overcharged/runs"
sys.path.insert(0, os.path.join(ROOT, "cavity"))
from edge_criteria import read_with_angle  # noqa
BAND = 4.0
H = os.path.dirname(os.path.abspath(__file__)); L = []
def say(s): print(s); L.append(s)
for setname in ("main_act0.07", "main_passive"):
    for k in range(2, 11):
        runs = sorted(glob.glob(f"{ROOT}/cavity/server_main/{setname}/rho*/k{k:02d}/s*/defects.dat.gz"))
        rows = []
        for f in runs:
            run = os.path.dirname(f); j = json.load(open(run + "/input.json")); bc = j["BC"][0]; R = bc["R"]
            if not os.path.exists(run + "/DONE"): continue
            t, fr = read_with_angle(run, np.array(bc["Q"][:2]), 0.2)
            Z = []
            for g in fr[::5]:
                p = g[g[:, 2] > 0]; r = np.hypot(p[:, 0], p[:, 1]); m = R - r < BAND
                if m.sum() == 0: continue
                Z.append(np.exp(1j * (2 * k - 2) * np.arctan2(p[m, 1], p[m, 0])).mean())
            if not Z: continue
            Z = np.mean(Z); rows.append((run.split("/")[-3], run.split("/")[-1], R, abs(Z), np.angle(Z)))
        if len(rows) < 3: continue
        strong = [r for r in rows if r[3] > 0.5]
        def rbar(rs): return abs(np.mean([np.exp(1j * r[4]) for r in rs])) if rs else np.nan
        say(f"{setname:13s} k={k:2d}: {len(rows):2d} runs, |Z| median {np.median([r[3] for r in rows]):.2f}; runs with |Z|>0.5: {len(strong):2d}, "
            f"phase agreement Rbar {rbar(strong):.2f} (random ~{1 / np.sqrt(max(len(strong), 1)):.2f}); "
            f"mean phase (2k-2)*phi0 = {np.angle(np.mean([np.exp(1j * r[4]) for r in strong])) if strong else np.nan:+.2f} rad")
        for r in rows:
            L.append(f"      {r[0]} {r[1]} R {r[2]:.1f}: |Z| {r[3]:.2f}, arg Z {r[4]:+.2f}")
open(os.path.join(H, "phase_test.txt"), "w").write("\n".join(L) + "\n")
