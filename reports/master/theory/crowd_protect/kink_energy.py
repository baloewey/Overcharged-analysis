#!/usr/bin/env python3
"""Elastic cost of a wrong-tilt kink in a wall layer vs N (user, 2026-10-22). Ring of N wall +1/2 at the depth-law radius + central +1, circle R 31.83 (P 200),
layer registry alpha = tau/N so the relaxed (lock) tilt of every defect is tau. Baseline: every wall core held at +tau (flow-fixed tilts).
Kink: m adjacent defects held at -tau instead,
  RIGID : at their positions (pure orientation flip, a twist against the lock);
  MOVED : shifted along the wall to their own -tau lock sites (angle -2 tau/N each), then held at -tau.
Rigid cores a = 1.5, field outside relaxed, one constant K 44 (tilted_structures.energy). tau = 30, 55, 65 deg; m = 1, 2; N = 6..28. Cost in kBT.
-> kink_energy.txt, kink_energy.npz"""
import os, numpy as np, importlib.util
from multiprocessing import Pool
T_ = "/home/baloewe/Physics/Overcharged/reports/master/theory"
def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path); m = importlib.util.module_from_spec(spec)
    src = open(path).read().split('\nif __name__ == "__main__":')[0]; exec(compile(src, path, "exec"), m.__dict__); return m
TS = load("ts", T_ + "/ejection/tilted_structures.py")
R = 31.830989; K = 44.0; NS = list(range(6, 29, 2)); TAUS = [30, 55, 65]; MS = [1, 2]
def build(N, tau, m, mode):
    rs = R * ((N - 1 + 4) / (3 * N - 1 + 4)) ** (1 / (2 * N)) - 0.61; t = np.radians(tau); ph = 2 * np.pi * np.arange(N) / N + t / N
    tgt = np.full(N, t)
    if m:
        tgt[:m] = -t
        if mode == "moved": ph[:m] -= 2 * t / N
    zs = list(rs * np.exp(1j * ph)) + [0j]; qs = [0.5] * N + [1.0]; t0 = TS.natural_tilts(zs, qs, R, N)
    return zs, qs, list(np.angle(np.exp(1j * (tgt - t0))) / 2) + [0.0]
def one(job):
    N, tau, m, mode = job; zs, qs, c = build(N, tau, m, mode); return job, TS.energy(zs, qs, c, R, K, K)
if __name__ == "__main__":
    jobs = [(N, tau, 0, "rigid") for N in NS for tau in TAUS] + [(N, tau, m, md) for N in NS for tau in TAUS for m in MS for md in ("rigid", "moved")]
    with Pool(int(os.environ.get("NPROC", 8))) as p: res = dict(p.map(one, jobs, chunksize=1))
    L = ["KINK COST (kBT): m adjacent wall +1/2 held at -tau in a layer held at +tau; one constant K 44, circle P 200, central +1", "N | rho | spacing | " + " | ".join(f"tau {t}: m1 rigid/moved, m2 rigid/moved" for t in TAUS)]
    out = {}
    for N in NS:
        row = []
        for tau in TAUS:
            e0 = res[(N, tau, 0, "rigid")]; v = [res[(N, tau, m, md)] - e0 for m in MS for md in ("rigid", "moved")]; out[f"N{N}_t{tau}"] = v
            row.append(" / ".join(f"{x:6.1f}" for x in v))
        L.append(f"{N:2d} | {N / (2 * np.pi * R):.3f} | {2 * np.pi * R / N:5.1f} | " + " | ".join(row))
    np.savez("kink_energy.npz", NS=NS, TAUS=TAUS, **out); open("kink_energy.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
