#!/usr/bin/env python3
"""Radial elastic stiffness of one wall +1/2 vs wall density, against the measured flow response (user, 2026-11-25).
Exact one-constant point-defect energy with wall images (ejection/tilted_structures.analytic_E, x pi K, K 44 kBT): ring of N = 2k-2 wall +1/2 + central +1
in a circle R 31.83; one defect displaced radially inward by delta, the others fixed. kappa_r = d2E/d delta^2 by finite differences at
 (a) the elastic equilibrium depth d* (ring force zero; the depth law) and (b) the MEASURED MPCD depth minus the effective-wall offset 0.61 (the depth law fits MPCD as measured = model + 0.61, 310 passive circles):
     sparse 1.80 -> 1.19, crowded 1.26 -> 0.65; uncorrected values listed for reference.
Densities as the flow-at-core comparison (aspect-2 ellipses k 5 / 14 = rho 0.04 / 0.12): circles k 4 / 12 at R 31.83 (approximation: circle vs ellipse).
Response of an overdamped defect in a harmonic radial well to a flow fluctuation of correlation time tc: velocity/flow ratio s = 1/(1 + kappa mu tc)
(mu = defect mobility, cells/(t.u. kBT/cell); measured free-defect value 0.0054, RESEARCH 2026-10). Measured: s 0.61 (sparse) vs 0.09 (crowded) at d 1-2.
-> radial_stiffness.txt"""
import os, numpy as np, importlib.util
T_ = "/home/baloewe/Physics/Overcharged/reports/master/theory"; MAIN = 'if __name__ == ' + '"__main__":'
spec = importlib.util.spec_from_file_location("ts", T_ + "/ejection/tilted_structures.py"); TS = importlib.util.module_from_spec(spec)
exec(compile(open(T_ + "/ejection/tilted_structures.py").read().split("\n" + MAIN)[0], "ts", "exec"), TS.__dict__)
H = os.path.dirname(os.path.abspath(__file__)); R = 31.830989; K = 44.0


def E(k, depth_all, delta):
    N = 2 * k - 2; rs = R - depth_all; zs = list(rs * np.exp(2j * np.pi * np.arange(N) / N)) + [0j]; qs = [0.5] * N + [1.0]
    zs[0] = (rs - delta) + 0j; return np.pi * K * TS.analytic_E(zs, qs, R)


def dstar(k):
    N = 2 * k - 2; return R - R * ((N - 1 + 4) / (3 * N - 1 + 4)) ** (1 / (2 * N))     # ring force zero (depth law without the 0.61 wall offset)


def kappa(k, d, h=0.05):
    return (E(k, d, h) - 2 * E(k, d, 0) + E(k, d, -h)) / h ** 2, -(E(k, d, h) - E(k, d, -h)) / (2 * h)


if __name__ == "__main__":
    L = ["RADIAL ELASTIC STIFFNESS of one wall +1/2 (one constant K 44, circle R 31.83, ring 2k-2 + central +1; others fixed)",
         "k (rho) | spacing | depth d | kappa_r = d2E/d delta^2 (kBT/cell^2) | force at d (+ = pushes inward) (kBT/cell)"]
    K_ = {}
    for k in (4, 12):
        N = 2 * k - 2; ds = dstar(k)
        for lab, d in (("equilibrium d*", ds), ("measured sparse 1.80 - 0.61", 1.19), ("measured crowded 1.26 - 0.61", 0.65), ("(uncorrected) 1.80", 1.8), ("(uncorrected) 1.26", 1.26)):
            kap, F = kappa(k, d); K_[(k, lab)] = kap
            L.append(f"{k:2d} ({2 * k / (2 * np.pi * R):.3f}) | {2 * np.pi * (R - d) / N:5.1f} | {lab:26s} {d:4.2f} | {kap:8.2f} | {F:+7.2f}")
    ks = K_[(4, "measured sparse 1.80 - 0.61")]; kc = K_[(12, "measured crowded 1.26 - 0.61")]
    L.append(f"\nAt the MEASURED depths (minus the 0.61-cell effective-wall offset of the depth law: measured = model + 0.61): kappa crowded / sparse = {kc:.1f} / {ks:.1f} = {kc / ks:.1f}x. At d*: {K_[(12, 'equilibrium d*')]:.1f} / {K_[(4, 'equilibrium d*')]:.1f} = {K_[(12, 'equilibrium d*')] / K_[(4, 'equilibrium d*')]:.1f}x")
    L.append("Response s = 1/(1 + kappa mu tc) vs measured 0.61 / 0.09 (ratio 6.8):")
    for mu in (0.0054, 0.02):
        for tc in (2.0, 10.0):
            s1, s2 = 1 / (1 + ks * mu * tc), 1 / (1 + kc * mu * tc); L.append(f"  mu {mu}, tc {tc:4.1f}: sparse {s1:.2f}, crowded {s2:.2f}, ratio {s1 / s2:.2f}")
    need = lambda: None
    # kappa mu tc needed for the measured ratio given the stiffness ratio r: (1 + r x)/(1 + x) = 6.8 -> x = 5.8/(r - 6.8)
    r = kc / ks; L.append(f"To get the measured ratio 6.8 from stiffness ratio {r:.1f}, need kappa_sparse mu tc = {5.8 / (r - 6.8):.2f}" if r > 6.8 else
                          f"Stiffness ratio {r:.1f} < 6.8: elasticity alone cannot give the measured response ratio, at any mobility or correlation time (max ratio = stiffness ratio).")
    open(os.path.join(H, "radial_stiffness.txt"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
