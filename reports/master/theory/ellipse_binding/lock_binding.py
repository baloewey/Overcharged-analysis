#!/usr/bin/env python3
"""Binding of wall +1/2 in ellipses when the tilt is produced by the lock (each defect displaced along the wall, natural core orientation, no clamping).
Sites: the 2k - 2 = 10 homeotropic points of the eMode-2 anchoring, psi_rel(s) = k 2 pi s/P - phi_n(s) = 0 (mod pi) (not evenly spaced in s). Tilted state: each
site displaced along the wall by ds_j such that 2 [psi_rel(s_j + ds_j) - psi_rel(s_j)] = tau (one vortex: same tau everywhere; two vortices: +tau for x < 0, -tau
for x > 0), tau = 63 deg; depth 1.6; bulk pair at +-0.45 A. Binding of site j: energy to move it 1.5 cells inward (natural cores, relaxed; binding.py machinery,
h = 0.35, mean over 4 grid offsets). Reported per site and as tip (|x| > A/2) / flank (|x| < A/4) means, untilted vs tilted. Writes lock_binding.txt."""
import os, sys, time
import numpy as np
from scipy.optimize import brentq
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
os.environ.setdefault("GRID_H", "0.35")
import binding as bd
L = []
def say(s): print(s, flush=True); L.append(s)
TAU = np.radians(63.0); k = 6
def psi_rel_fun(E):
    t, W, s, P, Nn = E.t, E.W, E.s, E.P, E.Nn
    phin = np.unwrap(np.angle(Nn)); psi = k * 2 * np.pi * s / P - phin
    return lambda sq: np.interp(sq % P, s, psi, period=P) + 2 * np.pi * (k - 1) * np.floor(sq / P), P
def sites(E, pattern):
    psi, P = psi_rel_fun(E); out = []
    for n in range(2 * k - 2):
        target = n * np.pi + psi(0.0)            # homeotropic points: psi_rel = psi(0) + n pi (registry phase 0 at s = 0, i.e. at the +x tip)
        sh = brentq(lambda sq: psi(sq) - target, -P / 2, 1.5 * P)
        sg = 0.0 if pattern == "untilted" else (1.0 if pattern == "one vortex" else (1.0 if np.real(E.W[np.searchsorted(E.s, sh % P) % len(E.s)]) < 0 else -1.0))
        sd = brentq(lambda sq: 2 * (psi(sq) - psi(sh)) - sg * TAU, sh - 15, sh + 15) if sg != 0 else sh
        out.append(sd % P)
    return np.array(out)
def config_at(E, svals, inward=None, delta=1.5):
    zs = []
    for j, sq in enumerate(svals):
        i = np.searchsorted(E.s, sq) % len(E.s); z = E.W[i] - bd.DEPTH * E.Nn[i]
        if inward == j: z = z - delta * E.Nn[i]
        zs.append(z)
    zs += [0.45 * E.A + 0j, -0.45 * E.A + 0j]; return zs, [0.5] * (2 * k), np.array(zs[:2 * k - 2])
for eps in (1.0001, 1.5, 2.0):
    offs = [(0.0, 0.0), (0.5 * bd.h, 0.0), (0.0, 0.5 * bd.h), (0.5 * bd.h, 0.5 * bd.h)]; Es = [bd.Ellipse(eps, o) for o in offs]; t0 = time.time()
    for pattern in ("untilted", "one vortex", "two vortices"):
        sv = sites(Es[0], pattern); per = []
        for E in Es:
            zs, qs, zw = config_at(E, sv); E0, _ = E.energy(zs, qs, [None] * 10); dE = []
            for j in range(10):
                zj, qj, _ = config_at(E, sv, inward=j); dE.append(E.energy(zj, qj, [None] * 10)[0] - E0)
            per.append(dE)
        per = np.array(per); dE = per.mean(0); err = per.std(0) / 2; tips = np.abs(zw.real) > Es[0].A / 2; fl = np.abs(zw.real) < Es[0].A / 4
        say(f"eps {eps:.1f} {pattern:12s}: binding (kBT) per site (x/A): " + " ".join(f"{z.real / Es[0].A:+.2f}:{d:.1f}+-{e:.1f}" for z, d, e in zip(zw, dE, err)) +
            f" | tips {dE[tips].mean():.1f}, flanks {dE[fl].mean() if fl.any() else np.nan:.1f}, all-site sd {dE.std():.1f}")
    say(f"   ({time.time() - t0:.0f} s)")
open(os.path.join(H, "lock_binding.txt"), "w").write("\n".join(L) + "\n")
