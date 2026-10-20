#!/usr/bin/env python3
"""MPCD density 40 vs 20 (2026-10-20): R 31.83, act 0.07, k 4 and 16; runs/density40 (2 seeds each) vs main sweep rho0.04/k04 and rho0.16/k16 (same R), measured with
robust_sets.one (tiltlaw.one + coherent-layer tilt). -> density40.txt"""
import glob, os, json, numpy as np
from multiprocessing import Pool
import robust_sets as RS
R = "/home/baloewe/Physics/Overcharged/runs"
sets = {"k4 density 40": sorted(glob.glob(R + "/density40/k04_d40_s*")), "k16 density 40": sorted(glob.glob(R + "/density40/k16_d40_s*")),
        "k4 density 20 (main sweep rho0.04)": sorted(d for d in glob.glob(R + "/cavity/server_main/main_act0.07/rho0.04/k04/s[01][0-9]") if os.path.exists(d + "/flowfield.dat.gz")),
        "k16 density 20 (main sweep rho0.16)": sorted(d for d in glob.glob(R + "/cavity/server_main/main_act0.07/rho0.16/k16/s[01][0-9]") if os.path.exists(d + "/flowfield.dat.gz"))}
L = ["set | n | wall +1/2 | bulk charge | coherent-layer tilt (C>0.9) | frac C>0.9 | coherence C | |wall v_t| | v_rms"]
with Pool(8) as p:
    for nm, rr in sets.items():
        res = [x for x in p.map(RS.one, rr) if x[1] is not None]; S = np.array([list(x[1][2:11]) + [x[2], x[3]] for x in res], float); m = np.nanmean(S, 0); e = np.nanstd(S, 0) / np.sqrt(max(len(S) - 1, 1))
        L.append(f"{nm:36s} | {len(S)} | {m[0]:.2f} | {m[1]:+.2f} | {m[9]:.1f} ± {e[9]:.1f} | {m[10]:.2f} | {m[4]:.2f} | {m[5]:.3f} | {m[6]:.3f}")
open("density40.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
