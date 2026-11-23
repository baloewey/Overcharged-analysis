#!/usr/bin/env python3
"""Elastic cost of the collective (registry) tilt WITH an inner flow-aligned ring as a boundary condition (user, 2026-11-23: in the runs the director just
inside the layer sits at ~ the Leslie angle to the flow; tilt_energy_landscape.py left the interior free, which made registry a zero mode).
Annulus r_i <= r < R (R 31.83): wall anchoring theta = k phi held on a ghost band outside R; inner circle (r < r_i = R - D) held at the spiral
theta = phi + pi/2 + alpha_i (winding +1). Wall ring of N = 2k-2 +1/2 at the depth-law radius rs, rotated by al (registry), positions FIXED, orientations FREE.
Field theta = Theta_s + psi, Theta_s = sum_j (1/2) arg(z - z_j) (ring defects only), psi smooth (initial psi = phi: correct winding at both boundaries).
One-constant nematic energy on a cut-free lattice form E = (K/2) sum_bonds sin^2(theta_i - theta_j) (-> (K/2) int |grad theta|^2; theta mod pi), bonds
touching a core disk (r < 1 of a defect) excluded; h = 0.5; L-BFGS. Tilt of each defect after relaxation: tail direction 2 c_j from the local core phase (theta - arg(z - z_j)/2 on a
ring of radius 1.25, free cells only) minus the outward normal.
Scan: k 4/6/8/10/12/16, D 4/6/10 cells (inner ring at r_i = R - D), alpha_i 30 deg (theta_L) and 25 deg (measured MPCD band), registry offsets spanning
one period (tau_target -150..180 step 15). -> per config: E(tau) (relative), tau* (minimum), stiffness; inner_bc.txt"""
import os, numpy as np, scipy.sparse as sp
from scipy.optimize import minimize
from multiprocessing import Pool
H = os.path.dirname(os.path.abspath(__file__)); R = 31.830989; K = 44.0; h = 0.5; RC = 1.0
x = np.arange(-R - 1.5, R + 1.5 + h / 2, h); X, Y = np.meshgrid(x, x, indexing="ij"); Z = X + 1j * Y; RR = np.abs(Z); PH = np.angle(Z); n = len(x)
idx = np.arange(X.size).reshape(X.shape)
bonds = np.r_[np.c_[idx[:-1, :].ravel(), idx[1:, :].ravel()], np.c_[idx[:, :-1].ravel(), idx[:, 1:].ravel()]]


def solve(k, D, alpha_i, al):
    N = 2 * k - 2; rs = R * ((N - 1 + 4) / (3 * N - 1 + 4)) ** (1 / (2 * N)) - 0.61; zs = rs * np.exp(1j * (2 * np.pi * np.arange(N) / N + al))
    Ths = sum(0.5 * np.angle(Z - zj) for zj in zs).ravel()
    wallc = (RR >= R) & (RR < R + 1.5); inner = RR < R - D; free = (RR < R) & ~inner
    fixed_theta = np.where(wallc, k * PH, np.where(inner, PH + np.pi / 2 + alpha_i, 0)).ravel()
    known = (wallc | inner).ravel(); use = (wallc | inner | free).ravel()
    core = np.zeros(X.size, bool)
    for zj in zs: core |= (np.abs(Z - zj) < RC).ravel()
    b = bonds[use[bonds[:, 0]] & use[bonds[:, 1]] & ~core[bonds[:, 0]] & ~core[bonds[:, 1]]]
    fi = np.flatnonzero(free.ravel()); pos = -np.ones(X.size, int); pos[fi] = np.arange(len(fi))
    Pm = sp.csr_matrix((np.ones(len(fi)), (fi, np.arange(len(fi)))), shape=(X.size, len(fi)))
    base = np.where(known, fixed_theta, Ths)                                   # fixed cells carry the total theta; free cells: Theta_s + psi
    def f(v):
        th = base + Pm @ v; d = th[b[:, 0]] - th[b[:, 1]]; E = 0.5 * K * np.sum(np.sin(d) ** 2)
        g = 0.5 * K * np.sin(2 * d); gr = np.bincount(b[:, 0], g, X.size) - np.bincount(b[:, 1], g, X.size); return E, Pm.T @ gr
    v0 = PH.ravel()[fi].copy()
    res = minimize(f, v0, jac=True, method="L-BFGS-B", options=dict(maxiter=40000, gtol=1e-9, ftol=1e-15, maxcor=30))
    th = (base + Pm @ res.x).reshape(X.shape); tilts = []
    for zj in zs:                                                              # local core phase: theta ~ (1/2) arg(z - z_j) + c_j near the core -> tail at 2 c_j
        ang = np.linspace(-np.pi, np.pi, 48, endpoint=False); px = zj + 1.25 * np.exp(1j * ang)
        ii = np.clip(np.round((px.real - x[0]) / h).astype(int), 0, n - 1); jj = np.clip(np.round((px.imag - x[0]) / h).astype(int), 0, n - 1)
        okc = free[ii, jj]                                                     # only free (relaxed) cells, never the fixed wall band
        c2 = np.angle(np.mean(np.exp(2j * (th[ii, jj][okc] - 0.5 * ang[okc])))) if okc.any() else np.nan
        tilts.append(np.angle(np.exp(1j * (c2 - np.angle(zj)))))
    tl = np.array(tilts); tau = np.angle(np.mean(np.exp(1j * tl)))
    return res.fun, float(np.degrees(tau)), float(np.abs(np.mean(np.exp(1j * tl)))), res.success


def job(j):
    k, D, ai, tt = j; N = 2 * k - 2; e, tau, coh, ok = solve(k, D, np.radians(ai), np.radians(tt) / N); return (k, D, ai, tt, e, tau, coh, ok)


if __name__ == "__main__":
    KS = (4, 6, 8, 10, 12, 16); DS = (4, 6, 10, 25); AIS = (30, 25); TTS = range(-150, 181, 15)
    jobs = [(k, D, ai, tt) for k in KS for D in DS for ai in AIS for tt in TTS]
    with Pool(int(os.environ.get("NPROC", 10))) as p: res = p.map(job, jobs, chunksize=1)
    L = ["COLLECTIVE REGISTRY TILT WITH A FLOW-ALIGNED INNER RING (annulus, one constant K 44, wall theta = k phi, inner spiral at r_i = R - D, angle alpha_i)",
         "per (k, D, alpha_i): tau_target -> (measured layer tilt, E - E_min per defect); tau* = measured tilt at the minimum; 2 alpha_i = flat-wall prediction"]
    for k in KS:
        N = 2 * k - 2; a = 2 * np.pi * R / N
        for ai in AIS:
            for D in DS:
                g = sorted([r for r in res if r[0] == k and r[1] == D and r[2] == ai], key=lambda r: r[3]); E = np.array([r[4] for r in g]); i = int(np.argmin(E))
                # stiffness: per-defect energy rise 20 deg either side of the minimum (quadratic fit over measured tilts within +-35 deg of tau*)
                tt_ = np.array([r[5] for r in g]); ee = (E - E.min()) / N; m = np.abs(tt_ - g[i][5]) < 35
                cf = np.polyfit(np.radians(tt_[m] - g[i][5]), ee[m], 2) if m.sum() >= 3 else [np.nan] * 3
                L.append(f"k {k:2d} (spacing ~{a:4.1f}) alpha_i {ai} D {D:2d}: tau* = {g[i][5]:+6.1f} (2 alpha_i = {2 * ai}; coh {g[i][6]:.2f}); stiffness/defect {2 * cf[0]:5.1f} kBT/rad^2 | " +
                         "  ".join(f"{r[5]:+.0f}:{(r[4] - E.min()) / N:.2f}" for r in g) + ("" if all(r[7] for r in g) else "  [some not converged]"))
    open(os.path.join(H, "inner_bc.txt"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
