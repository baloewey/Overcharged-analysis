#!/usr/bin/env python3
"""Defect Langevin model, k = 6, R = 31.83: scan self-propulsion v0 and Stokeslet strength f/eta. Per run (T = 4000, last half):
mean wall count (depth < 4), fraction of time with exactly 2k-2 = 10 at the wall, |circulation| (fluid, sample ring), circulation sign
persistence, mean wall tilt (tail vs outward normal) signed with the circulation. Writes scan_v0_f.txt."""
import itertools, os, sys
import numpy as np
from multiprocessing import Pool
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
from defect_langevin import simulate
def run(args):
    v0, fh, seed = args
    o = simulate(6, 31.83, v0, fh, T=4000.0, seed=seed); h = len(o["t"]) // 2
    nw = o["nwall"][h:]; c = o["circ"][h:]; tl = o["tilt"][h:]
    sgn = np.sign(c); tsig = np.nanmean(tl * sgn)
    return (v0, fh, seed, nw.mean(), np.mean(nw == 10), np.mean(np.abs(c)), abs(np.mean(sgn)), tsig, np.mean(np.diff(sgn) != 0))
if __name__ == "__main__":
    grid = list(itertools.product((0.02, 0.05, 0.08, 0.12, 0.16, 0.2, 0.3), (0.0, 0.5, 1.5, 3.0), (1,)))
    with Pool(7) as pool: res = pool.map(run, grid)
    L = ["v0    f/eta  seed | <n_wall>  P(n_wall=10)  <|circ|>  |<sign circ>|  tilt*sign(circ)  sign flips/frame"]
    for r in res: L.append(f"{r[0]:5.2f} {r[1]:5.1f} {r[2]:4d} | {r[3]:7.2f} {r[4]:10.2f} {r[5]:10.4f} {r[6]:10.2f} {r[7]:+12.2f} {r[8]:12.3f}")
    print("\n".join(L)); open(os.path.join(H, "scan_v0_f.txt"), "w").write("\n".join(L) + "\n")
