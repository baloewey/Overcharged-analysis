#!/usr/bin/env python3
"""Elastic stiffness of the wall layer against spacing (user, 2026-10-15). Circle R = 31.83 (1.43 l_a, as in the MPCD k-scan), wall ring N = 2k-2 +1/2 at
the depth-law radius plus a central +1, k = 3-20 (spacing 50 -> 5.3 cells). Two channels by which a tilt can be resisted elastically:
 (1) REGISTRY (position) channel: a rigid shift of the ring along the wall = tilt through the lock (tau = (2k-2) alpha). Zero mode for K1 = K3; with
     K1 != K3 a pinning potential (first order in Delta = (K1-K3)/2, fluct/anisotropy_pinning.dF, smooth core 1 cell). Stiffness per defect
     K_reg = d^2 E/d tau^2 / N at the minimum (kBT/rad^2).
 (2) ORIENTATION (twist) channel (user: the log repulsion assumes orientations relaxed as fast variables; unrelaxed defects exert torques on each other,
     and likely with the wall): every core rotated by delta beyond its relaxed orientation at FIXED positions (rigid cores radius 1.5 whose director is the
     relaxed field + delta/2; field outside relaxed; ejection/tilted_structures.energy), one-constant K 32.5 and anisotropic 42/23.1.
     Stiffness per defect K_tw = d^2 E/d delta^2 / N (kBT/rad^2), from delta = 0, +-10, +-20 deg.
-> crowd_elastic.txt"""
import os, sys, numpy as np
from multiprocessing import Pool
T_ = "/home/baloewe/Physics/Overcharged/reports/master/theory"; sys.path.insert(0, T_ + "/fluct"); sys.path.insert(0, T_ + "/ejection")
import importlib.util
def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path); m = importlib.util.module_from_spec(spec)
    src = open(path).read(); src = src.split("\nL = []")[0] if "anisotropy" in path else src.split('\nif __name__ == "__main__":')[0]
    exec(compile(src, path, "exec"), m.__dict__); return m
AP = load("ap", T_ + "/fluct/anisotropy_pinning.py"); TS = load("ts", T_ + "/ejection/tilted_structures.py")
R = 31.830989


def ring(k, al=0.0):
    N = 2 * k - 2; rs = R * ((N - 1 + 4) / (3 * N - 1 + 4)) ** (1 / (2 * N)) - 0.61
    return list(rs * np.exp(1j * (2 * np.pi * np.arange(N) / N + al))) + [0j], [0.5] * N + [1.0], N, rs


def reg(k):
    N = 2 * k - 2; per = 2 * np.pi / N; als = np.linspace(-per / 2, per / 2, 25); E = []
    for al in als:
        zs, qs, _, _ = ring(k, al); E.append(AP.dF(zs, qs, R, a=1.0))
    E = np.array(E); i0 = int(np.argmin(E)); tau = N * als
    d2 = (E[min(i0 + 1, 24)] - 2 * E[i0] + E[max(i0 - 1, 0)]) / (tau[1] - tau[0]) ** 2
    return dict(k=k, Kreg=d2 / N, tau_min=np.degrees(tau[i0]), barrier=E.max() - E.min())


def tw(job):
    k, K1, K3 = job; zs, qs, N, rs = ring(k); out = []
    for d in np.radians([-20, -10, 0, 10, 20]):
        c = [d / 2] * N + [0.0]; out.append(TS.energy(zs, qs, c, R, K1, K3))
    E = np.array(out); dd = np.radians(10); d2 = (E[1] - 2 * E[2] + E[3]) / dd ** 2; d2b = (E[0] - 2 * E[2] + E[4]) / (2 * dd) ** 2
    return dict(k=k, K1=K1, K3=K3, Ktw=d2 / N, Ktw20=d2b / N, rs=rs, asym=(E[3] - E[1]) / (2 * dd) / N)


if __name__ == "__main__":
    ks = [3, 4, 5, 6, 8, 10, 12, 15, 20]
    with Pool(6) as pool:
        RG = pool.map(reg, ks); TW = pool.map(tw, [(k, K1, K3) for k in ks for K1, K3 in ((32.5, 32.5), (42.0, 23.1))], chunksize=1)
    L = ["ELASTIC STIFFNESS OF THE WALL LAYER vs SPACING (circle R 31.83; ring 2k-2 + central +1; per wall defect)",
         "k | spacing | registry channel (Delta K, first order): stiffness K_reg (kBT/rad^2 of tilt), minimum at tau, barrier (ring) | twist channel K_tw one-constant / anisotropic (kBT/rad^2; +-10 / +-20 deg) | torque at delta = 0 (kBT/rad, aniso)"]
    for g in RG:
        k = g["k"]; a = [t for t in TW if t["k"] == k and t["K1"] == 32.5][0]; b = [t for t in TW if t["k"] == k and t["K1"] == 42.0][0]
        L.append(f"{k:2d} | {2 * np.pi * a['rs'] / (2 * k - 2):5.1f} | {g['Kreg']:7.1f} at {g['tau_min']:+5.1f} deg, barrier {g['barrier']:6.1f} | {a['Ktw']:7.1f}/{a['Ktw20']:7.1f}  {b['Ktw']:7.1f}/{b['Ktw20']:7.1f} | {b['asym']:+.1f}")
    print("\n".join(L)); open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "crowd_elastic.txt"), "w").write("\n".join(L) + "\n")
