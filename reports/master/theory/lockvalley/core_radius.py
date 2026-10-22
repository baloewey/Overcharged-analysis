#!/usr/bin/env python3
"""Sensitivity of the lock-valley energetics to the core treatment (user, 2026-10-22): rigid-core radius a = 1.0, 1.5 (default), 2.5 cells,
one constant K = 44, k = 6 and 16 (circle R 31.83, ring 2k-2 + central +1 as lock_valley.py). Twist stiffness (delta = 0, +-30 at tau_lock 0)
and the cost of moving at fixed tilt by a quarter spacing (tau_lock = +-30 with delta = -+30). -> core_radius.txt"""
import numpy as np, importlib.util
from multiprocessing import Pool
T_ = "/home/baloewe/Physics/Overcharged/reports/master/theory"
spec = importlib.util.spec_from_file_location("ts", T_ + "/ejection/tilted_structures.py"); TS = importlib.util.module_from_spec(spec)
exec(compile(open(T_ + "/ejection/tilted_structures.py").read().split('\nif __name__ == "__main__":')[0], "ts", "exec"), TS.__dict__)
R = 31.830989
def one(job):
    k, a, tl, dl = job; N = 2 * k - 2; rs = R * ((N - 1 + 4) / (3 * N - 1 + 4)) ** (1 / (2 * N)) - 0.61
    zs = list(rs * np.exp(1j * (2 * np.pi * np.arange(N) / N + np.radians(tl) / N))) + [0j]; qs = [0.5] * N + [1.0]
    return job, TS.energy(zs, qs, [np.radians(dl) / 2] * N + [0.0], R, 44.0, 44.0, a=a) / N
if __name__ == "__main__":
    jobs = [(k, a, tl, dl) for k in (6, 16) for a in (1.0, 1.5, 2.5) for tl, dl in ((0, 0), (0, 30), (0, -30), (30, -30), (-30, 30))]
    with Pool(6) as p: E = dict(p.map(one, jobs))
    L = ["CORE RADIUS SENSITIVITY (one constant K 44; per wall defect, kBT): twist stiffness d2E/ddelta2 at tau_lock 0 | cost of moving a quarter spacing at fixed tilt (mean of +-)"]
    for k in (6, 16):
        for a in (1.0, 1.5, 2.5):
            e = lambda tl, dl: E[(k, a, tl, dl)]; tw = (e(0, 30) + e(0, -30) - 2 * e(0, 0)) / np.radians(30) ** 2; mv = 0.5 * (e(30, -30) + e(-30, 30)) - e(0, 0)
            L.append(f"k {k:2d} core radius {a:.1f}: twist {tw:6.1f} kBT/rad^2 | fixed-tilt move {mv:5.2f} kBT")
    open("core_radius.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
