#!/usr/bin/env python3
"""Elastic force on one wall +1/2 vs the number of wall defects N (user, 2026-10-22). Analytic one constant K 44 (point charges + wall images; validated against
the relaxed field in drag_energy.py within a few kBT), circle R 31.83 (P 200, rho = N/P approx), N wall +1/2 evenly spaced + central +1 (bulk charge 1).
(1) RING EQUILIBRIUM: radial force on each defect when the whole ring sits at wall distance d (symmetric); its zero = equilibrium depth d*(N); compared with the
    depth-law formula R - R((N+3)/(3N+3))^(1/2N) and with the measured mean resident depth in the shape-sweep circles (events.npz) at the same rho.
(2) SINGLE DEFECT: ring at d*, one defect pulled radially inward by delta, others fixed: restoring force F(delta) = dE/ddelta (> 0 pulls it back to the wall),
    initial stiffness, mean force over 0-6 cells (= E(6)/6), and the line-charge estimate pi^2 K rho / 2.
-> force_vs_N.txt, force_vs_N.npz"""
import csv, numpy as np
R = 31.830989; K = 44.0; P = 2 * np.pi * R; NS = np.arange(4, 31, 2); DEL = np.linspace(0, 8, 161)
def E_probe(zp, others, qo):          # pi K units: probe (1/2) with others + images, self-image
    zo = np.asarray(others); qo = np.asarray(qo)
    return -np.sum(qo * (np.log(np.abs(zp - zo)) + np.log(np.abs(R ** 2 - zp * np.conj(zo))))) - 0.25 * np.log(R ** 2 - abs(zp) ** 2)
def ring_force(N, d, h=1e-4):          # radial force (outward +) on defect 0 when every wall defect is at radius R - d; kBT/cell
    def E(r):                            # total energy of the symmetric ring at radius r, per defect: only probe terms needed with others moved too -> use full sum / N
        z = r * np.exp(2j * np.pi * np.arange(N) / N); q = np.full(N, 0.5); Z = np.r_[z, 0j]; Q = np.r_[q, 1.0]; e = 0.0
        for i in range(N + 1):
            for j in range(N + 1):
                if i != j: e -= Q[i] * Q[j] * np.log(abs(Z[i] - Z[j]))
                e -= Q[i] * Q[j] * np.log(abs(R ** 2 - Z[i] * np.conj(Z[j])))
        return e
    r = R - d; return -np.pi * K * (E(r + h) - E(r - h)) / (2 * h) / N
L = ["ELASTIC FORCE ON ONE WALL +1/2 vs N (analytic one constant K 44; circle R 31.83, P 200; central +1)"]
# measured depths, circles
S = "/home/baloewe/Physics/Overcharged/reports/master/theory/shapes"; Z = np.load(S + "/escape_events/events.npz")["E"]
meta = list(csv.DictReader(open(S + "/escape_events/runs.csv"))); SR = {r["run"]: r for r in csv.DictReader(open(S + "/shape_runs.csv"))}
rho_r = np.array([float(SR[m["run"]]["rho"]) for m in meta]); asp_r = np.array([float(SR[m["run"]]["aspect"]) for m in meta]); esc_r = np.array([float(SR[m["run"]]["esc"]) for m in meta])
lev = np.argmin(np.abs(rho_r[:, None] - np.array([0.04, 0.06, 0.09, 0.12, 0.14])), 1); run = Z[:, 0].astype(int)
meas = [(r0, np.mean(Z[np.isin(run, np.flatnonzero((lev == i) & (asp_r < 1.05))), 5])) for i, r0 in enumerate((0.04, 0.06, 0.09, 0.12, 0.14))]
L.append("measured mean resident depth, sweep circles: " + ", ".join(f"rho {a:.2f}: {b:.2f}" for a, b in meas))
L.append("\nN | rho=N/P | spacing | d* (force zero) | depth-law formula | ring stiffness dF/dd at d* | single defect: F at delta 0.5, 1, 2, 4, 6 | stiffness | E(6) | mean F 0-6 | pi^2 K rho/2")
out = dict(NS=NS, DEL=DEL)
for N in NS:
    ds = np.linspace(0.3, 12, 400); f = np.array([ring_force(N, d) for d in ds]); s = np.flatnonzero(np.diff(np.sign(f)) != 0)
    dstar = ds[s[0]] - f[s[0]] * (ds[s[0] + 1] - ds[s[0]]) / (f[s[0] + 1] - f[s[0]]) if len(s) else np.nan
    kr = (ring_force(N, dstar + 0.01) - ring_force(N, dstar - 0.01)) / 0.02
    form = R - R * ((N + 3) / (3 * N + 3)) ** (1 / (2 * N))
    z = (R - dstar) * np.exp(2j * np.pi * np.arange(N) / N); oth = list(z[1:]) + [0j]; qo = [0.5] * (N - 1) + [1.0]
    E = np.pi * K * np.array([E_probe(z[0] * (R - dstar - dl) / (R - dstar), oth, qo) for dl in DEL]); E -= E[0]; F = np.gradient(E, DEL)
    at = lambda x: np.interp(x, DEL, F)
    out[f"F_N{N}"] = F; out[f"E_N{N}"] = E; out[f"ring_N{N}"] = np.c_[ds, f]
    L.append(f"{N:2d} | {N / P:.3f} | {P / N:5.1f} | {dstar:5.2f} | {form:5.2f} | {kr:+7.1f} | " + " ".join(f"{at(x):5.1f}" for x in (0.5, 1, 2, 4, 6)) +
             f" | {F[1]:.1f}->{(F[3] - F[1]) / (DEL[3] - DEL[1]):+.1f}/cell | {np.interp(6, DEL, E):6.1f} | {np.interp(6, DEL, E) / 6:5.1f} | {np.pi ** 2 * K * N / P / 2:5.1f}")
out["dstar"] = np.array([float(l.split("|")[3]) for l in L[3:]]); out["meas"] = np.array(meas)
np.savez("force_vs_N.npz", **out); open("force_vs_N.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
