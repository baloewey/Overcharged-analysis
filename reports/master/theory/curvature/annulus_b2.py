#!/usr/bin/env python3
"""Link B of curvature_theory.md solved in the annulus (user, 2026-11-23). Axisymmetric texture theta = phi + pi/2 + alpha(r) in the measured azimuthal
flow v(r) (circulation-signed class means from curvature_profiles.npz). Steady director equation WITH advection (D theta/Dt = (v/r) d_phi theta = v/r) and
MPCD's shear susceptibility chi scaling the Jeffery rotation (not the advection):
    D (alpha'' + alpha'/r) = v/r - chi omega/2 + chi lambda e cos 2 alpha,     omega = v' + v/r,  e = e_rtheta = (v' - v/r)/2,  ' = d/dr
Without elasticity: cos 2 alpha_J = [1 + (1/chi - 1)(1 - Lambda)] / lambda, Lambda = |omega|/(2|e|) (chi = 1: alpha_J = theta_L exactly, any curvature).
Parameters (measured, none fitted): D = K/gamma_1 = 0.73 cells^2/t.u. (rot_visc.txt, U 10), lambda 2, chi 0.5.
BCs: alpha'(r_L) = 0 at the layer edge y_L (zero mode: the layer passes no torque); alpha(r_c) = measured alpha at y_c (interior value, data).
Variants: y_L 1.25 / 2.25, y_c 6.25 / 8.25; chi = 1 (advection exactly cancels curvature) as the control.
-> annulus_b2.txt"""
import os, numpy as np
from scipy.integrate import solve_bvp
from scipy.interpolate import UnivariateSpline
H = os.path.dirname(os.path.abspath(__file__)); Z = np.load(os.path.join(H, "curvature_profiles.npz")); YC = Z["YC"]; LA = 22.21; LAM = 2.0; D = 0.73


def classes():
    for lab, cond in (("sparse", Z["rho"] < 0.055), ("crowded", Z["rho"] >= 0.085)):
        for lo, hi in ((0.8, 1.2), (1.2, 1.7), (1.7, 2.2)):
            m = cond & (Z["R"] / LA >= lo) & (Z["R"] / LA < hi)
            if m.sum() >= 3: yield f"{lab} R/l_a {lo}-{hi}", Z["R"][m].mean(), np.nanmean(Z["v"][m], 0), np.nanmean(Z["alpha"][m], 0), int(m.sum())


def solve(R, v, alpha_meas, chi, yL, yc):
    r_grid = R - YC; o = np.argsort(r_grid); sp = UnivariateSpline(r_grid[o], v[o], s=len(v) * (0.01 * np.nanmax(np.abs(v))) ** 2, k=4)
    def rhs(r, Y):
        a, ap = Y; vv, dv = sp(r), sp.derivative()(r); om = dv + vv / r; e = 0.5 * (dv - vv / r)
        return np.vstack([ap, -ap / r + (vv / r - chi * om / 2 + chi * LAM * e * np.cos(2 * a)) / D])
    rc, rL = R - yc, R - yL; ac = np.radians(np.interp(yc, YC, alpha_meas))
    bc = lambda Ya, Yb: np.array([Ya[0] - ac, Yb[1]])
    rr = np.linspace(rc, rL, 200); guess = np.vstack([np.full_like(rr, np.radians(28)), np.zeros_like(rr)])
    s = solve_bvp(rhs, bc, rr, guess, tol=1e-6, max_nodes=20000)
    return s, sp


def alphaJ(R, v, chi, y):
    r_grid = R - YC; o = np.argsort(r_grid); sp = UnivariateSpline(r_grid[o], v[o], s=len(v) * (0.01 * np.nanmax(np.abs(v))) ** 2, k=4)
    r = R - y; vv, dv = sp(r), sp.derivative()(r); om = dv + vv / r; e = 0.5 * (dv - vv / r); c = (vv / r - chi * om / 2) / (-chi * LAM * e)
    return np.degrees(0.5 * np.arccos(c)) if abs(c) <= 1 else np.nan


if __name__ == "__main__":
    YS = [1.25, 2.25, 3.25, 4.25, 5.25, 6.25, 7.25]
    L = ["ANNULUS SOLUTION OF LINK B (D = K/gamma_1 = 0.73, lambda 2; alpha in degrees from the tangent; measured = curvature_profiles class means)",
         "columns per depth y: measured | local alpha_J(chi) without elasticity | annulus solution (alpha'(y_L) = 0, alpha(y_c) = measured)"]
    for name, R, v, am, n in classes():
        L.append(f"\n{name} ({n} runs, mean R {R:.1f})   y: " + "  ".join(f"{y:5.2f}" for y in YS))
        L.append("  measured                 " + "  ".join(f"{np.interp(y, YC, am):+5.1f}" for y in YS))
        for chi in (0.5, 1.0):
            L.append(f"  alpha_J chi {chi}           " + "  ".join(f"{alphaJ(R, v, chi, y):+5.1f}" for y in YS))
            for yL, yc in ((1.25, 8.25), (2.25, 8.25), (1.25, 6.25)):
                s, _ = solve(R, v, am, chi, yL, yc)
                vals = [np.degrees(s.sol(R - y)[0]) if yL <= y <= yc else np.nan for y in YS]
                L.append(f"  annulus chi {chi} yL {yL} yc {yc} {'ok ' if s.success else 'FAIL'} " + "  ".join(f"{x:+5.1f}" if np.isfinite(x) else "   - " for x in vals))
    open(os.path.join(H, "annulus_b2.txt"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
