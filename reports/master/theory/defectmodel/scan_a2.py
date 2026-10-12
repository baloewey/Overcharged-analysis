#!/usr/bin/env python3
"""Structured defect model (k = 6, B = 90 kBT, spiral 0.48, v0 0.15, a_w 2) with and without the second harmonic a2 of the registry energy
(beyond first order in K1 - K3, relax_frank2.py: +4..+7 kBT). Does a locked tilt appear? Per run (T = 3000, last half): |tilt| percentiles,
fraction of time with |tilt| in 45-80 deg, layer slip (deg per 1000 t.u.), circulation sign changes, |circ|. Writes scan_a2.txt."""
import itertools, os, sys
import numpy as np
from multiprocessing import Pool
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
from defect_langevin import simulate
def run(a):
    Rla, a2, sp = a
    o = simulate(6, Rla * 22.21, 0.15, 1.5, T=3000.0, seed=1, spiral=sp, structured=True, barrier=90.0, wall_supp=2.0, a2=a2); h = len(o["t"]) // 2
    tl = np.degrees(o["tilt"][h:]); c = o["circ"][h:]; z = np.array(o["z"][h:])[:, :10]; ph = np.unwrap(np.angle(z), axis=0).mean(1)
    T = o["t"][-1] - o["t"][h]
    return (Rla, a2, sp, *np.percentile(np.abs(tl), [25, 50, 75]), np.mean((np.abs(tl) > 45) & (np.abs(tl) < 80)), np.degrees(abs(ph[-1] - ph[0])) / T * 1000,
            int(np.sum(np.diff(np.sign(c)) != 0)), np.mean(np.abs(c)))
if __name__ == "__main__":
    grid = list(itertools.product((0.6, 0.8, 1.0, 1.43), (0.0, 5.0, 10.0, 20.0), (0.48, 0.24)))
    with Pool(4) as pool: res = pool.map(run, grid)
    L = ["R/l_a  a2  spiral | |tilt| deg p25/p50/p75 | frac 45-80 deg | slip deg/1000tu | sign flips | <|circ|>"]
    for r in res: L.append(f"{r[0]:5.2f} {r[1]:4.0f} {r[2]:5.2f} | {r[3]:5.0f} {r[4]:5.0f} {r[5]:5.0f} | {r[6]:6.2f} | {r[7]:7.0f} | {r[8]:4d} | {r[9]:.3f}")
    print("\n".join(L)); open(os.path.join(H, "scan_a2.txt"), "w").write("\n".join(L) + "\n")
