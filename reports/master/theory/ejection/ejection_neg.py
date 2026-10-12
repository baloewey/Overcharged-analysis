#!/usr/bin/env python3
"""Negative winding (k = -6): does a common tilt make the wall gain two -1/2 (from bulk pair creation)? (user, 2026-10-08)
  A: 2|k| = 12 -1/2 on a ring (passive state). B: 2|k-1| = 14 -1/2 on a ring (the commensurate number = homeotropic points of (k-1) phi) plus two +1/2 at +-s
  (the +1 left in the bulk by the two created pairs). Total charge -6 in both. B has 4 more cores: the elastic energy outside the cores is compared, and the
  core cost of the extra pair, 4 E_core (unknown), must be added to Delta E.
Tilt of a -1/2: arm direction psi = (2/3) theta_rest (mod 2 pi/3); tilt = psi - phi (mod 120 deg). Rotating the core director by c rotates the arms by (2/3) c,
so a common tilt tau needs c_j = (3/2)(tau - tau_j^0) (wrapped, mod pi). B first takes the tilt by a registry shift (lock slope measured numerically), the residual goes on its cores. Bulk +1/2 keep their natural orientation. Rigid cores, relaxed director, one-constant and
anisotropic; s = 0.45 R; core radius 1.5, 2.5. Ring radii: the one-constant equilibrium radius of each ring (minimised numerically in the image energy, the
bulk pair fixed). Writes ejection_neg.txt."""
import os, sys
import numpy as np
from scipy.optimize import minimize_scalar
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
from tilted_structures import energy, fields, analytic_E
L = []
def say(s): print(s, flush=True); L.append(s)
R = 28.2095; k = -6
def ring(N, rr, al=0.0): return list(rr * np.exp(1j * (2 * np.pi * (np.arange(N) + 0.5) / N + al)))
def best_radius(N, extra_z, extra_q):
    f = lambda rr: analytic_E(ring(N, rr) + extra_z, [-0.5] * N + extra_q, R)
    return minimize_scalar(f, bounds=(0.5 * R, R - 0.8), method="bounded").x
def natural_tilts_neg(zs, qs, nw, eps=1e-3):
    out = []
    for j in range(nw):
        zj = zs[j]; p = zj + eps * zj / abs(zj)
        e2, _, _ = fields(np.array([p.real]), np.array([p.imag]), zs, qs, R)
        th_rest2 = e2[0] + np.angle(p - zj)                       # 2 theta_rest: remove own contribution 2q arg = -arg
        arm = (2.0 / 3.0) * (th_rest2 / 2.0)                      # arm direction (one of three)
        out.append(np.angle(np.exp(3j * (arm - np.angle(zj)))) / 3.0)   # tilt mod 120 deg in (-60, 60]
    return np.array(out)
def cores(zs, qs, nw, tau):
    t0 = natural_tilts_neg(zs, qs, nw); d = np.angle(np.exp(3j * (tau - t0))) / 3.0       # tilt change mod 120
    return list(1.5 * d) + [0.0] * (len(zs) - nw), t0
s = 0.45 * R
rA = best_radius(12, [], []) - 0.61; rB = best_radius(14, [s + 0j, -s + 0j], [0.5, 0.5]) - 0.61
say(f"ring radii (one-constant equilibrium, minus the 0.61 offset): A (12 -1/2) depth {R - rA:.2f}, B (14 -1/2 + bulk pair) depth {R - rB:.2f}")
zA = ring(12, rA); qA = [-0.5] * 12; zB = ring(14, rB) + [s + 0j, -s + 0j]; qB = [-0.5] * 14 + [0.5, 0.5]
say(f"natural -1/2 tilts (deg, mod 120): A {np.round(np.degrees(natural_tilts_neg(zA, qA, 12)), 0)}; B {np.round(np.degrees(natural_tilts_neg(zB, qB, 14)), 0)}")
say(f"one-constant image energy (no cores, pi K = {np.pi * 32.5:.0f} kBT): E_B - E_A = {np.pi * 32.5 * (analytic_E(zB, qB, R) - analytic_E(zA, qA, R)):+.0f} kBT")
for lab, K1, K3 in (("one-constant", 32.5, 32.5), ("anisotropic", 42.0, 23.1)):
    for a in (1.5, 2.5):
        EA0 = energy(zA, qA, [0.0] * 12, R, K1, K3, a=a); EB0 = energy(zB, qB, [0.0] * 16, R, K1, K3, a=a); row = []
        for td in (0, 30, 60, 90):
            tau = np.radians(td); cA, _ = cores(zA, qA, 12, tau)
            # B reaches a common tilt by a registry shift (lock), then the residual is put on the cores
            slope = -14.0 / 3.0                                   # lock for the -1/2 ring: d tilt / d alpha = (2/3)(k + 1/2) - 1 = -14/3 (checked numerically)
            al = np.angle(np.exp(3j * (tau - np.radians(60.0)))) / 3 / slope   # natural uniform tilt is 60 deg (mod 120)
            zBt = ring(14, rB, al) + [s + 0j, -s + 0j]; cB, _ = cores(zBt, qB, 14, tau)
            EA = energy(zA, qA, cA, R, K1, K3, a=a); EB = energy(zBt, qB, cB, R, K1, K3, a=a); row.append((td, EA - EA0, EB - EB0, EB - EA))
        say(f"{lab:12s} core {a}: natural states E_B0 - E_A0 = {EB0 - EA0:+.0f} kBT (+ 4 core energies) | common tilt: " +
            " ".join(f"{td}deg: A frustration {fa:+.0f}, B {fb:+.0f}, B-A {d:+.0f};" for td, fa, fb, d in row))
open(os.path.join(H, "ejection_neg.txt"), "w").write("\n".join(L) + "\n")
