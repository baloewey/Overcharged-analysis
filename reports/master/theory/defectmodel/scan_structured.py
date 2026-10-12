#!/usr/bin/env python3
"""Defect model, structured (2k-2 held layer + free bulk pair kept out of the wall band), anisotropy pinning B, spiral drive, wall suppression
a_w = 2. Does the effective-boundary loop give a coherent tilt and circulation, and where is the onset in R? k = 6 (B = 90 kBT) and k = 4
(B = 130), R/l_a = 0.3-2.0, spiral drive 0 vs 0.48, v0 = 0.15, f/eta = 1.5, T = 3000 (last half). Writes scan_structured.txt."""
import itertools, os, sys
import numpy as np
from multiprocessing import Pool
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
from defect_langevin import simulate
LA = 22.21
def run(a):
    k, B, Rla, sp = a
    o = simulate(k, Rla * LA, 0.15, 1.5, T=3000.0, seed=1, spiral=sp, structured=True, barrier=B, wall_supp=2.0); h = len(o["t"]) // 2
    c, tl = o["circ"][h:], o["tilt"][h:]
    return (k, Rla, sp, np.mean(np.abs(c)), abs(np.mean(np.sign(c))), np.nanmean(np.abs(tl)), np.nanmean(tl * np.sign(c)))
if __name__ == "__main__":
    grid = [(6, 90.0, r, s) for r in (0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.43, 2.0) for s in (0.0, 0.48)] + [(4, 130.0, r, s) for r in (0.4, 0.6, 0.8, 1.0, 1.43) for s in (0.0, 0.48)]
    with Pool(4) as pool: res = pool.map(run, grid)
    L = [" k  R/l_a spiral | <|circ|> |<sign circ>| <|tilt|> deg  tilt*sign(circ) deg"]
    for r in res: L.append(f"{r[0]:2d} {r[1]:5.2f} {r[2]:5.2f} | {r[3]:8.3f} {r[4]:10.2f} {np.degrees(r[5]):10.1f} {np.degrees(r[6]):+12.1f}")
    print("\n".join(L)); open(os.path.join(H, "scan_structured.txt"), "w").write("\n".join(L) + "\n")
