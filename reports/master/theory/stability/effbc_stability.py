#!/usr/bin/env python3
"""Linear stability of the untilted wall layer through the effective-boundary loop (2026-10-07 night).
Loop: registry shift alpha of the layer -> tilt tau = (2k-2) alpha (measured lock) -> effective boundary of the bulk theta = phi + c,
c = 90 deg + tau/2 (effective_bc.py) -> extensile active force of the +1 spiral f_phi = (zeta/eta) S sin(tau) / r (spiral_force check)
acting in the smooth region r_c < r < R_e = R - 0.34 l -> axisymmetric Stokes flow v(r) with no slip at R -> the layer at r_s = R - d is
advected: d alpha/dt = v(r_s)/r_s; restoring: anisotropy pinning, d alpha/dt = -(mu E''/(N r_s^2)) alpha (anisotropy_pinning.py:
E'' = registry curvature in kBT/rad^2; N defects of mobility mu, rigid rotation).
Linearised: d alpha/dt = lambda alpha, lambda = (2k-2) (zeta/eta) S h - mu E''/(N r_s^2), h = v(r_s)/(r_s g) for g = 1 (dimensionless
response of the layer to a unit spiral drive). Growth (lambda > 0) = spontaneous tilt + circulation.
Calibration: zeta/eta from the director-field Stokes fit (c = 0.42-0.67), S = 0.876, mu = 0.02, E'' from anisotropy_pinning (k = 4, 6, 8 at
R = 31.83; scaled ~ unchanged with R). Core radius of the bulk +1 (no spiral drive inside): r_c = a l_a with a = 0.25 (assumed; the
measured bulk +1/2 sit at |x| ~ 0.4 A, i.e. the pair is spread over ~ l_a/2). Scan R/l_a and k; also the threshold zeta/eta.
Axisymmetric Stokes: v'' + v'/r - v/r^2 = -g/r on r_c < r < R_e, 0 elsewhere; v(0) = 0, v(R) = 0, v and v' continuous. Solved numerically."""
import os
import numpy as np
from scipy.integrate import solve_bvp
H = os.path.dirname(os.path.abspath(__file__)); LA = 22.21
def response(R, Re, rc, rs):
    # piecewise analytic: homogeneous A r + B/r; particular in the forced band: -(1/2) r ln r
    # unknowns: inner (r<rc): a1 r ; band: a2 r + b2/r - r ln r/2 ; outer: a3 r + b3/r with a3 R + b3/R = 0
    def vp(r): return -0.5 * r * np.log(r)
    def dvp(r): return -0.5 * (np.log(r) + 1)
    M = np.zeros((4, 4)); y = np.zeros(4)
    # unknowns x = [a1, a2, b2, a3]; b3 = -a3 R^2
    # continuity at rc: a1 rc = a2 rc + b2/rc + vp(rc); a1 = a2 - b2/rc^2 + dvp(rc)
    M[0] = [rc, -rc, -1 / rc, 0]; y[0] = vp(rc)
    M[1] = [1, -1, 1 / rc ** 2, 0]; y[1] = dvp(rc)
    # continuity at Re: a2 Re + b2/Re + vp(Re) = a3 (Re - R^2/Re); a2 - b2/Re^2 + dvp(Re) = a3 (1 + R^2/Re^2)
    M[2] = [0, Re, 1 / Re, -(Re - R ** 2 / Re)]; y[2] = -vp(Re)
    M[3] = [0, 1, -1 / Re ** 2, -(1 + R ** 2 / Re ** 2)]; y[3] = -dvp(Re)
    a1, a2, b2, a3 = np.linalg.solve(M, y)
    v = a3 * (rs - R ** 2 / rs) if rs >= Re else a2 * rs + b2 / rs + vp(rs)
    return v / rs
Epp = {4: 2386.0, 6: 4221.0, 8: 4581.0}            # kBT/rad^2, core a = 1 (anisotropy_pinning.txt, active ring, R = 31.83)
mu, S = 0.02, 0.876; L = []
def say(s): print(s); L.append(s)
say("lambda = (2k-2) (zeta/eta) S h - mu E''/(N r_s^2); h = layer advection per unit spiral drive (dimensionless)")
for zeta in (0.42, 0.55, 0.67):
    say(f"\nzeta/eta = {zeta}:")
    for k in (4, 6, 8):
        N = 2 * k - 2; row = []
        for Rla in (0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.5, 2.0, 2.5):
            R = Rla * LA; ell = 2 * np.pi * R / N; d = R * (1 - ((N - 1 + 4) / (3 * N - 1 + 4)) ** (1 / (2 * N))) + 0.61
            rs = R - d; Re = R - 0.34 * ell; rc = 0.25 * LA
            if Re <= rc: row.append(f"R/la {Rla}: no smooth region"); continue
            h = response(R, Re, rc, rs); lam = N * zeta * S * h - mu * Epp[k] / (N * rs ** 2)
            row.append(f"R/la {Rla}: h {h:+.4f}, drive {N * zeta * S * h:+.4f}, pin {mu * Epp[k] / (N * rs ** 2):.4f}, lambda {lam:+.4f}")
        say(f"  k={k}: " + " | ".join(row))
# threshold activity for k = 6 vs R
say("\nthreshold zeta/eta (lambda = 0), k = 6, r_c = 0.25 l_a:")
for Rla in (0.4, 0.5, 0.6, 0.8, 1.0, 1.5, 2.0):
    k = 6; N = 10; R = Rla * LA; ell = 2 * np.pi * R / N; d = R * (1 - ((N - 1 + 4) / (3 * N - 1 + 4)) ** (1 / (2 * N))) + 0.61; rs = R - d; Re = R - 0.34 * ell; rc = 0.25 * LA
    if Re <= rc: say(f"  R/la {Rla}: no smooth region (R_e {Re:.1f} <= r_c {rc:.1f})"); continue
    h = response(R, Re, rc, rs); say(f"  R/la {Rla}: h {h:+.4f}; zeta/eta threshold {mu * Epp[k] / (N * rs ** 2) / (N * S * h) if h > 0 else np.inf:.4f}")
open(os.path.join(H, "effbc_stability.txt"), "w").write("\n".join(L) + "\n")
