#!/usr/bin/env python3
"""Tilt potential E(tau) beyond first order in K1 - K3 (relax_frank2.run, rigid cores, energy outside the cores) ACROSS CROWDING (user, 2026-10-16:
the earlier analysis was done at one k and R; redo it for more crowded walls). Harmonics E = A1 cos tau + A2 cos 2tau + A3 cos 3tau (kBT, whole ring).
Scans: (i) R = 31.83, k = 3, 4, 5, 6, 8, 10, 12, 15, 20 (spacing 50 -> 5.3 cells); (ii) k = 6 at R = 22.21, 31.83, 44.42 (spacing 14 -> 28: separates size from
spacing). Variants: A (cores fixed, a 1.5, point bulk +1), A-pair (bulk +1 as two +1/2 at +-0.45 R, as split in MPCD), B (cores free, a 1.5).
Prediction with the master's torque balance (drive D sin tau): cos tau* = (|A1| - D)/(4 A2), unstable for D > |A1| - 4A2, slides for D > |A1| + 4A2.
-> k_scan.txt, k_scan.npz"""
import os, sys, time, numpy as np
from multiprocessing import Pool
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
import importlib.util
spec = importlib.util.spec_from_file_location("rf2", os.path.join(H, "relax_frank2.py")); RF = importlib.util.module_from_spec(spec)
src = open(os.path.join(H, "relax_frank2.py")).read().split('if __name__ == "__main__":')[0]; exec(compile(src, "relax_frank2.py", "exec"), RF.__dict__)
cfg0 = RF.config
def cfg_pair(k, R, tau, bulk="point", s=3.0):
    return cfg0(k, R, tau, "pair", s=0.45 * R) if bulk == "pair45" else cfg0(k, R, tau, bulk, s)
RF.config = cfg_pair
def one(job):
    k, R, var = job; kw = dict(k=k, R=R, a=1.5, ntau=16)
    if var == "A_pair": kw["bulk"] = "pair45"
    if var == "B": kw["free_cores"] = True
    t0 = time.time(); r = RF.run(**kw); a1 = RF.harm(r[:, 2]); return dict(k=k, R=R, var=var, A=a1, conv=int(r[:, 5].sum()), sec=time.time() - t0, E=r[:, 2])
if __name__ == "__main__":
    J = [(k, 31.830989, v) for v in ("A", "A_pair", "B") for k in (3, 4, 5, 6, 8, 10, 12, 15, 20)] + [(6, R, v) for v in ("A", "A_pair", "B") for R in (22.21, 44.42)]
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 6) as pool: res = pool.map(one, J, chunksize=1)
    np.savez(os.path.join(H, "k_scan.npz"), res=np.array(res, dtype=object), allow_pickle=True)
    L = ["TILT POTENTIAL vs CROWDING (rigid cores a 1.5, h 0.5, K1 42 / K3 23.1): E = A1 cos tau + A2 cos 2tau + A3 cos 3tau (kBT, whole ring); per defect in brackets",
         "variant | k | R | spacing | A1 | A2 | A3 | |A1| - 4 A2 (onset drive) | |A1| + 4 A2 (sliding) | conv"]
    for r in sorted(res, key=lambda x: (x["var"], x["R"], x["k"])):
        N = 2 * r["k"] - 2; A1, A2, A3 = r["A"]
        L.append(f"{r['var']:6s} | {r['k']:2d} | {r['R']:5.2f} | {2*np.pi*r['R']/N:5.1f} | {A1:+8.2f} ({A1/N:+6.2f}) | {A2:+7.2f} ({A2/N:+6.3f}) | {A3:+6.2f} | {abs(A1)-4*A2:+8.2f} | {abs(A1)+4*A2:+8.2f} | {r['conv']}/16")
    print("\n".join(L)); open(os.path.join(H, "k_scan.txt"), "w").write("\n".join(L) + "\n")
