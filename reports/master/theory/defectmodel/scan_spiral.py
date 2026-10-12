#!/usr/bin/env python3
"""Defect Langevin model + effective-boundary spiral flow (user's request, 2026-10-08). k = 6, R = 31.83 (R/l_a = 1.43), last half of T = 4000.
Grid: self-propulsion v0, Stokeslet f/eta = 1.5, spiral drive (zeta/eta) S in {0, 0.48}, wall suppression of self-propulsion a_w in {0, 2}.
Per run: mean wall count, P(wall = 10), |circulation|, persistence of its sign, tilt x sign(circulation) (tail convention, rad).
Writes scan_spiral.txt."""
import itertools, os, sys
import numpy as np
from multiprocessing import Pool
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
from defect_langevin import simulate
def run(a):
    v0, sp, aw = a
    o = simulate(6, 31.83, v0, 1.5, T=4000.0, seed=1, spiral=sp, wall_supp=aw); h = len(o["t"]) // 2
    nw, c, tl = o["nwall"][h:], o["circ"][h:], o["tilt"][h:]
    return (v0, sp, aw, nw.mean(), np.mean(nw == 10), np.mean(np.abs(c)), abs(np.mean(np.sign(c))), np.nanmean(tl * np.sign(c)))
if __name__ == "__main__":
    grid = list(itertools.product((0.05, 0.1, 0.15, 0.2, 0.3), (0.0, 0.48), (0.0, 2.0)))
    with Pool(4) as pool: res = pool.map(run, grid)
    L = ["v0   spiral a_w | <n_wall> P(n_wall=10) <|circ|> |<sign>| tilt*sign"]
    for r in res: L.append(f"{r[0]:4.2f} {r[1]:5.2f} {r[2]:3.0f} | {r[3]:7.2f} {r[4]:9.2f} {r[5]:8.3f} {r[6]:7.2f} {r[7]:+8.2f}")
    print("\n".join(L)); open(os.path.join(H, "scan_spiral.txt"), "w").write("\n".join(L) + "\n")
