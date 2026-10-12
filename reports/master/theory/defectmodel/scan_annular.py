#!/usr/bin/env python3
"""Annular-drive particle model (measured effective drive profile, scaled with sin tau / sin 63; no defect Stokeslets), structured
(10 wall +1/2 held, bulk pair kept inside), anisotropy pinning B = 90 kBT and second harmonic a2. Closure test at R/l_a = 1.43 (the size of the
calibration): does the tilt lock, and where? Also drive scale 0.5-1.5 and other sizes (profile scaled with r/R, amplitude with R0/R).
Per run (T = 3000, last half): |tilt| p25/p50/p75 (deg), fraction in 45-80 deg, slip deg/1000 t.u., circulation sign flips, <|circ|>.
Writes scan_annular.txt."""
import itertools, os, sys
import numpy as np
from multiprocessing import Pool
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
from defect_langevin import simulate
def run(a):
    Rla, a2, sc = a
    o = simulate(6, Rla * 22.21, 0.15, 0.0, T=3000.0, seed=1, structured=True, barrier=90.0, wall_supp=2.0, a2=a2, annular=sc); h = len(o["t"]) // 2
    tl = np.degrees(o["tilt"][h:]); c = o["circ"][h:]; z = np.array(o["z"][h:])[:, :10]; ph = np.unwrap(np.angle(z), axis=0).mean(1); T = o["t"][-1] - o["t"][h]
    return (Rla, a2, sc, *np.percentile(np.abs(tl), [25, 50, 75]), np.mean((np.abs(tl) > 45) & (np.abs(tl) < 80)), np.degrees(abs(ph[-1] - ph[0])) / T * 1000,
            int(np.sum(np.diff(np.sign(c)) != 0)), np.mean(np.abs(c)))
if __name__ == "__main__":
    grid = [(1.43, a2, sc) for a2 in (0.0, 5.0, 10.0, 20.0) for sc in (0.5, 1.0, 1.5)] + [(Rla, a2, 1.0) for Rla in (0.6, 0.8, 1.0, 2.0) for a2 in (0.0, 5.0)]
    with Pool(4) as pool: res = pool.map(run, grid)
    L = ["R/l_a  a2  drive | |tilt| deg p25/p50/p75 | frac 45-80 | slip deg/1000tu | flips | <|circ|>"]
    for r in res: L.append(f"{r[0]:5.2f} {r[1]:4.0f} {r[2]:5.2f} | {r[3]:5.0f} {r[4]:5.0f} {r[5]:5.0f} | {r[6]:5.2f} | {r[7]:7.0f} | {r[8]:4d} | {r[9]:.3f}")
    print("\n".join(L)); open(os.path.join(H, "scan_annular.txt"), "w").write("\n".join(L) + "\n")
