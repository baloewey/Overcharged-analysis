#!/usr/bin/env python3
"""Drag energy with relaxing neighbours (2026-10-22; follows the drive x barrier test, whose measured excursion profiles do not collapse on the frozen-layer
drag energy). Analytic one-constant energy (K 44, point charges + wall images; as force_vs_N.py), circle R 31.83 (P 200), ring of N wall +1/2 at the
depth-law radius + central +1. Probe = defect 0 pulled radially inward by delta. FROZEN: other wall defects fixed. RELAXED: at each delta the other N-1
slide along the wall (angles minimised, radius fixed) to minimise the total energy (L-BFGS-B from the previous delta). -> drag_relaxed.txt, drag_relaxed.npz"""
import numpy as np
from scipy.optimize import minimize
R = 31.830989; K = 44.0; DEL = np.linspace(0, 10, 51); NS = [8, 12, 18, 24, 28]
def energy(z, q):
    d = z[:, None] - z[None, :]; np.fill_diagonal(d, 1.0); im = R ** 2 - z[:, None] * np.conj(z[None, :]); qq = q[:, None] * q[None, :]
    e = -np.sum(qq * np.log(np.abs(d))) + np.sum(np.diag(qq) * 0) - np.sum(qq * np.log(np.abs(im)))
    return np.pi * K * e
L = ["DRAG WITH RELAXING NEIGHBOURS (one constant K 44, analytic): E(delta) - E(0), kBT; frozen vs relaxed neighbours; peak restoring force (kBT/cell)"]
out = {"DEL": DEL}
for N in NS:
    rs = R * ((N + 3) / (3 * N + 3)) ** (1 / (2 * N)); ph0 = 2 * np.pi * np.arange(N) / N; q = np.r_[np.full(N, 0.5), 1.0]
    def E(phis, dl):
        z = np.r_[(rs - dl) + 0j, rs * np.exp(1j * phis), 0j]; return energy(z, q)
    fro, rel = [], []; x = ph0[1:].copy()
    for dl in DEL:
        fro.append(E(ph0[1:], dl)); r_ = minimize(lambda p: E(p, dl), x, method="L-BFGS-B"); x = r_.x; rel.append(r_.fun)
    fro = np.array(fro) - fro[0]; rel = np.array(rel) - rel[0]; Ff = np.gradient(fro, DEL); Fr = np.gradient(rel, DEL)
    out[f"frozen_N{N}"] = fro; out[f"relaxed_N{N}"] = rel
    L.append(f"N {N:2d} (rho {N / 200:.2f}): frozen E(2,4,6,10) = " + ", ".join(f"{np.interp(v, DEL, fro):6.1f}" for v in (2, 4, 6, 10)) + f"; peak F {Ff.max():5.1f} | relaxed E(2,4,6,10) = "
             + ", ".join(f"{np.interp(v, DEL, rel):6.1f}" for v in (2, 4, 6, 10)) + f"; peak F {Fr.max():5.1f} at delta {DEL[np.argmax(Fr)]:.1f}; F at 6 cells {np.interp(6, DEL, Fr):5.1f}")
np.savez("drag_relaxed.npz", **out); open("drag_relaxed.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
