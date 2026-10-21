#!/usr/bin/env python3
"""Elastic energy of the wall layer with defect orientation as a slow variable (user, 2026-10-21): the lock valley.
Circle R = 31.83, ring of N = 2k-2 wall +1/2 at the depth-law radius plus a central +1 (as tiltcrowd/crowd_elastic.py). Two coordinates per layer:
  registry shift alpha (rigid rotation of the ring along the wall): the lock tilt is tau_lock = N alpha (the relaxed orientation of a defect at that
  position, ejection/tilted_structures.natural_tilts);
  orientation offset delta: every wall core rotated by delta beyond its relaxed orientation at fixed positions (rigid cores radius 1.5, core director =
  relaxed field + delta/2; field outside relaxed; tilted_structures.energy).
The actual tilt is tau = tau_lock + delta. Elasticity: one constant K = 44 kBT (bare, coupling off) and K1/K3 = 43.1/21.9 (effective, coupling on).
Grid: tau_lock = -180..180 deg in 30-deg steps (alpha = tau_lock/N), delta = -60, -30, 0, 30, 60 deg; k = 4, 6, 10, 16.
Energy per wall defect relative to (alpha = 0, delta = 0). Expect: flat along delta = 0 (zero mode) for one constant; parabolic in delta (twist).
-> lock_valley.txt, lock_valley.npz"""
import os, sys, numpy as np
from multiprocessing import Pool
T_ = "/home/baloewe/Physics/Overcharged/reports/master/theory"
import importlib.util
def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path); m = importlib.util.module_from_spec(spec)
    src = open(path).read().split('\nif __name__ == "__main__":')[0]; exec(compile(src, path, "exec"), m.__dict__); return m
TS = load("ts", T_ + "/ejection/tilted_structures.py")
R = 31.830989; KS = [4, 6, 10, 16]; TL = np.arange(-180, 181, 30); DL = [-60, -30, 0, 30, 60]; ELS = [(44.0, 44.0), (43.1, 21.9)]
def ring(k, al):
    N = 2 * k - 2; rs = R * ((N - 1 + 4) / (3 * N - 1 + 4)) ** (1 / (2 * N)) - 0.61
    return list(rs * np.exp(1j * (2 * np.pi * np.arange(N) / N + al))) + [0j], [0.5] * N + [1.0], N
def one(job):
    k, tl, dl, K1, K3 = job; N = 2 * k - 2; zs, qs, _ = ring(k, np.radians(tl) / N)
    c = [np.radians(dl) / 2] * N + [0.0]; return job, TS.energy(zs, qs, c, R, K1, K3) / N
if __name__ == "__main__":
    jobs = [(k, tl, dl, K1, K3) for K1, K3 in ELS for k in KS for tl in TL for dl in DL]
    with Pool(6) as p: res = dict(p.map(one, jobs, chunksize=1))
    L = ["LOCK VALLEY: energy per wall defect (kBT) relative to tau_lock = 0, delta = 0; rows tau_lock (registry, deg), columns delta (orientation offset, deg)"]
    E = {}
    for K1, K3 in ELS:
        for k in KS:
            M = np.array([[res[(k, tl, dl, K1, K3)] for dl in DL] for tl in TL]); M -= M[list(TL).index(0), DL.index(0)]; E[f"K{K1:g}_{K3:g}_k{k}"] = M
            L.append(f"\nK1 {K1} K3 {K3}, k {k} (N {2*k-2}): columns delta = {DL}")
            for tl, row in zip(TL, M): L.append(f"  tau_lock {tl:+4d}: " + " ".join(f"{v:8.2f}" for v in row))
            tw = (M[:, DL.index(30)] + M[:, DL.index(-30)] - 2 * M[:, DL.index(0)]) / np.radians(30) ** 2
            L.append(f"  valley (delta 0) range {np.ptp(M[:, DL.index(0)]):.2f} kBT per defect; twist stiffness d2E/ddelta2 {tw.min():.1f}-{tw.max():.1f} kBT/rad^2 per defect")
    open("lock_valley.txt", "w").write("\n".join(L) + "\n"); np.savez("lock_valley.npz", TL=TL, DL=DL, **E); print("\n".join(L))
