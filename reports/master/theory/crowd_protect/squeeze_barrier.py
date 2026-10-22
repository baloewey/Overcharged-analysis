#!/usr/bin/env python3
"""Does squeezing lower the climb? (2026-10-22). Analytic one-constant energy (K 44, point charges + wall images; as force_vs_N.py), circle P 200, ring of
N wall +1/2 at the equilibrium depth + central +1. Probe = defect 0, displaced along the wall towards neighbour 1 by a fraction f of the spacing (others fixed),
then dragged radially inward. Reports E(6) - E(start) and the peak restoring force, and the energy of the squeezed state itself. -> squeeze_barrier.txt"""
import numpy as np
R = 31.830989; K = 44.0; P = 2 * np.pi * R; DEL = np.linspace(0, 8, 161)
def ring_d(N):
    import importlib.util as iu; s = iu.spec_from_file_location("f", "force_vs_N.py"); m = iu.module_from_spec(s)
    src = open("force_vs_N.py").read().split('\nL = ["ELASTIC')[0]; exec(src, m.__dict__); return m
FV = ring_d(0)
L = ["SQUEEZE vs CLIMB (one constant K 44, circle P 200): probe shifted along the wall towards a neighbour by f x spacing, then dragged inward",
     "N | rho | f | squeezed-state energy (kBT, rel. f 0) | E(6) - E(start) | peak restoring force (kBT/cell)"]
for N in (12, 24, 28):
    ds = np.linspace(0.3, 6, 300); fz = np.array([FV.ring_force(N, d) for d in ds]); j = np.flatnonzero(np.diff(np.sign(fz)))[0]; dstar = ds[j]
    r0 = R - dstar; ph = 2 * np.pi * np.arange(N) / N; z = r0 * np.exp(1j * ph); oth = list(z[1:]) + [0j]; qo = [0.5] * (N - 1) + [1.0]; e0 = None
    for f in (0.0, 0.15, 0.3, 0.4):
        zp0 = r0 * np.exp(1j * f * 2 * np.pi / N); E = np.pi * K * np.array([FV.E_probe(zp0 * (r0 - dl) / r0, oth, qo) for dl in DEL]); es = E[0]
        if e0 is None: e0 = es
        E = E - E[0]; F = np.gradient(E, DEL)
        L.append(f"{N:2d} | {N / P:.3f} | {f:.2f} | {es - e0:+6.1f} | {np.interp(6, DEL, E):6.1f} | {F[DEL <= 8].max():5.1f}")
open("squeeze_barrier.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
