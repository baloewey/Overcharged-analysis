#!/usr/bin/env python3
"""zeta 0.07 partners of act_robust measured with the same function (robust_sets.one = tiltlaw.one + coherent-layer tilt): main sweep circle rho0.06/k06
(R 31.83, R/l_a 1.43; s01-s10, needs flowfield) and shape sweep ellipse eps2.00/Req1.5/k07 (s01-s03). -> zeta_ref.txt"""
import glob, os, numpy as np
from multiprocessing import Pool
import robust_sets as RS
R = "/home/baloewe/Physics/Overcharged/runs"
sets = {"circle_z0.07": sorted(d for d in glob.glob(R + "/cavity/server_main/main_act0.07/rho0.06/k06/s[01][0-9]") if os.path.exists(d + "/flowfield.dat.gz")),
        "ellipse_z0.07": sorted(glob.glob(R + "/shapes_server/ellipse/eps2.00/Req1.5/k07/s0[1-3]"))}
L = ["group | n | wall +1/2 | bulk | tilt mean | coherence C | coherent-layer tilt | frac C>0.9 | |wall v_t| | v_rms | P(sign kept)"]
with Pool(6) as p:
    for nm, runs in sets.items():
        res = [r for r in p.map(RS.one, runs) if r[1] is not None]
        S = np.array([list(r[1][2:11]) + [r[2], r[3]] for r in res], float); m = np.nanmean(S, 0); e = np.nanstd(S, 0) / np.sqrt(max(len(S) - 1, 1))
        L.append(f"{nm:14s} | {len(S)} | {m[0]:5.2f}±{e[0]:.2f} | {m[1]:+.2f}±{e[1]:.2f} | {m[3]:5.1f}±{e[3]:.1f} | {m[4]:.2f} | {m[9]:5.1f}±{e[9]:.1f} | {m[10]:.2f} | {m[5]:.3f} | {m[6]:.3f} | {m[7]:.2f}")
open("zeta_ref.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
