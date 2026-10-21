#!/usr/bin/env python3
"""Site selection in the continuum (user, 2026-10-21): does the active one-constant Q-tensor disk put its wall defects at the lock-predicted registry,
as MPCD does? Same estimator as fluct/pinning/ic_test_active.py: wall +1/2 (within 6 cells of the wall; here the defect detector's wall band),
N = 2k-2 of them, Z = mean exp(i N phi), tau = circular mean of (div Q direction - phi), lock-corrected registry D = arg Z - tau (0 = homeotropic
site shifted by tau/N, i.e. the MPCD result). Anchoring theta = k phi in both codes. Final frames of every saved continuum disk field (activity scan
zeta 3-8, imposed-swirl controls), one-constant. -> registry_continuum.txt"""
import glob, numpy as np
import shear_tilt as M
L = ["field | k | wall +1/2 (2k-2) | arg Z (rad) | tau (deg) | D = arg Z - tau (rad) | |Z|"]
for f in sorted(glob.glob("fields_*.npz")):
    F = np.load(f); k = 4 if "k04" in f else 16 if "k16" in f else None
    if k is None: continue
    M.k = k; q1, q2 = F["q1"], F["q2"]; R = float(F["R"]); N = 2 * k - 2
    D = M.defects(q1, q2); w = [(px, py, a) for px, py, q, a in D if q > 0 and np.hypot(px, py) > R - 6]
    if not w: L.append(f"{f} | {k} | 0"); continue
    ph = np.array([np.arctan2(py, px) for px, py, _ in w]); ang = np.array([a for _, _, a in w])
    z = np.mean(np.exp(1j * N * ph)); tau = np.angle(np.mean(np.exp(1j * (ang - ph))))
    L.append(f"{f:34s} | {k:2d} | {len(w):2d} ({N}) | {np.angle(z):+.2f} | {np.degrees(tau):+6.1f} | {np.angle(np.exp(1j * (np.angle(z) - tau))):+.2f} | {abs(z):.2f}")
open("registry_continuum.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
