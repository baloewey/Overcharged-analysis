#!/usr/bin/env python3
"""Registry/tilt potential E(tau) of the active wall layer beyond first order in K1 - K3 (user's question, 2026-10-08).
Director theta = Theta + psi: Theta = exact one-constant field of the defects (same-sign images; Dirichlet anchoring theta = k phi on |z| = R,
lc.c), psi single-valued, psi = 0 on the wall, relaxed by minimising the anisotropic Frank energy
  F = 1/2 int w [ Kbar |grad theta|^2 + Delta (splay^2 - bend^2) ] dA,  splay^2 - bend^2 = -cos2theta (tx^2 - ty^2) - 2 sin2theta tx ty,
Kbar = (K1 + K3)/2, Delta = (K1 - K3)/2, w = smooth core cutoff prod_i (1 - exp(-|z - z_i|^2 / a^2)) (as anisotropy_pinning.py).
Grid: nodes with spacing h inside the disk; central differences (sparse); L-BFGS with the analytic gradient.
Configuration: k = 6, R = 31.83, 2k - 2 = 10 wall +1/2 at the active layer radius, rigidly shifted by alpha = tau/(2k - 2) (registry; tau = tilt of
the unrelaxed field), bulk +1 at the centre (one q = 1 core, or two +1/2 at +-s).
Output per tau: F0 (unrelaxed, = first order), F (relaxed), the relaxed tilt of the wall defects tau + 2 <psi(z_j)>; harmonics of E(tau).
Checks: K1 = K3 -> relaxed E(tau) flat (zero mode); F0 harmonics = first-order result (-44 cos tau at k = 6, a = 1)."""
import os, sys, time
import numpy as np
import scipy.sparse as sp
from scipy.optimize import minimize
H = os.path.dirname(os.path.abspath(__file__))
def setup(R, h):
    x = np.arange(-R, R + h / 2, h); X, Y = np.meshgrid(x, x, indexing="ij"); n = len(x)
    inside = np.hypot(X, Y) < R - 0.5 * h
    idx = -np.ones(X.shape, int); idx[inside] = np.arange(inside.sum())
    nI = inside.sum(); full = X.size
    P = sp.csr_matrix((np.ones(nI), (np.flatnonzero(inside.ravel()), np.arange(nI))), shape=(full, nI))
    I = sp.identity(n, format="csr"); D1 = sp.diags([-1, 1], [-1, 1], shape=(n, n)) / (2 * h)
    Dx = sp.kron(D1, I, format="csr"); Dy = sp.kron(I, D1, format="csr")
    return X, Y, inside, P, Dx, Dy
def fields(X, Y, zs, qs, R):
    z = X + 1j * Y; e2 = np.ones(z.shape, complex); gx = np.zeros(z.shape); gy = np.zeros(z.shape)
    for zi, q in zip(zs, qs):
        a = z - zi; a = np.where(np.abs(a) < 1e-9, 1e-9, a); b = 1 - z * np.conj(zi) / R ** 2
        e2 *= ((a / np.abs(a)) * (b / np.abs(b))) ** (2 * q)
        wv = 1 / a - (np.conj(zi) / R ** 2) / b
        gx += q * wv.imag; gy += q * wv.real
    return np.angle(e2), gx, gy                                       # 2 Theta (mod 2 pi), grad Theta
def energy_fn(Th2, Gx, Gy, w, P, Dx, Dy, h, Kb, De):
    Th2, Gx, Gy, w = Th2.ravel(), Gx.ravel(), Gy.ravel(), w.ravel()
    def f(psi):
        pf = P @ psi; gx = Gx + Dx @ pf; gy = Gy + Dy @ pf; C = np.cos(Th2 + 2 * pf); S = np.sin(Th2 + 2 * pf)
        dens = Kb * (gx * gx + gy * gy) + De * (-C * (gx * gx - gy * gy) - 2 * S * gx * gy)
        F = 0.5 * h * h * np.sum(w * dens)
        dgx = h * h * w * (Kb * gx + De * (-C * gx - S * gy)); dgy = h * h * w * (Kb * gy + De * (C * gy - S * gx))
        dpsi = h * h * w * De * (S * (gx * gx - gy * gy) - 2 * C * gx * gy)
        g = P.T @ (Dx.T @ dgx + Dy.T @ dgy + dpsi)
        return F, g
    return f
def config(k, R, tau, bulk="point", s=3.0):
    N = 2 * k - 2; rl = R * ((N - 1 + 4) / (3 * N - 1 + 4)) ** (1 / (2 * N)) - 0.61; al = tau / (2 * k - 2)
    zs = list(rl * np.exp(1j * (2 * np.pi * np.arange(N) / N + al))); qs = [0.5] * N
    if bulk == "point": zs.append(0j); qs.append(1.0)
    else: zs += [s + 0j, -s + 0j]; qs += [0.5, 0.5]
    return zs, qs
def run(k=6, R=31.830989, K1=42.0, K3=23.1, h=0.5, a=1.0, ntau=16, bulk="point"):
    X, Y, inside, P, Dx, Dy = setup(R, h); out = []
    for tau in np.linspace(0, 2 * np.pi, ntau, endpoint=False):
        zs, qs = config(k, R, tau, bulk)
        Th2, Gx, Gy = fields(X, Y, zs, qs, R)
        w = inside.astype(float)
        for zi in zs: w *= 1 - np.exp(-np.abs(X + 1j * Y - zi) ** 2 / a ** 2)
        f = energy_fn(Th2, Gx, Gy, w, P, Dx, Dy, h, 0.5 * (K1 + K3), 0.5 * (K1 - K3))
        psi0 = np.zeros(P.shape[1]); F0, _ = f(psi0)
        res = minimize(f, psi0, jac=True, method="L-BFGS-B", options=dict(maxiter=5000, gtol=1e-7))
        pf = (P @ res.x).reshape(X.shape)
        # relaxed tilt of the wall defects: tail direction shifts by 2 psi at the core (bilinear interpolation)
        from scipy.ndimage import map_coordinates
        x0 = X[0, 0]; zw = np.array(zs[:2 * k - 2])
        psi_at = map_coordinates(pf, [(zw.real - x0) / h, (zw.imag - x0) / h], order=1)
        out.append((tau, F0, res.fun, tau + 2 * np.mean(psi_at), res.nit, res.success))
    return np.array(out, dtype=float)
def harm(t, E):
    c = np.fft.rfft(E) / len(E) * 2; return c[1].real, c[2].real, c[3].real
if __name__ == "__main__":
    L = []
    def say(s): print(s, flush=True); L.append(s)
    for lab, kw in (("check K1 = K3 (zero mode), h 0.5, a 1", dict(K1=32.5, K3=32.5)),
                    ("K1 42, K3 23.1, h 0.5, a 1", dict()),
                    ("K1 42, K3 23.1, h 0.35, a 1", dict(h=0.35)),
                    ("K1 42, K3 23.1, h 0.5, a 2", dict(a=2.0)),
                    ("K1 42, K3 23.1, h 0.5, a 1, bulk = two +1/2 at +-3", dict(bulk="pair"))):
        t0 = time.time(); r = run(**kw); t, F0, F, trel = r[:, 0], r[:, 1], r[:, 2], r[:, 3]
        a0 = harm(t, F0); a1 = harm(t, F)
        say(f"{lab}: unrelaxed harmonics cos tau {a0[0]:+.2f}, cos 2tau {a0[1]:+.2f}, cos 3tau {a0[2]:+.2f} | relaxed: cos tau {a1[0]:+.2f}, "
            f"cos 2tau {a1[1]:+.2f}, cos 3tau {a1[2]:+.2f} kBT | relaxation lowers F by {np.mean(F0 - F):.1f} kBT; relaxed tilt - registry tilt "
            f"(rms) {np.degrees(np.sqrt(np.mean(np.angle(np.exp(1j * (trel - t))) ** 2))):.1f} deg; converged {int(r[:, 5].sum())}/{len(r)}; {time.time() - t0:.0f} s")
        np.savetxt(os.path.join(H, "E_tau_" + lab.split(":")[0].replace(" ", "_").replace(",", "").replace("=", "").replace("+-", "pm").replace("/", "") + ".txt"), r,
                   header="tau F_unrelaxed F_relaxed tilt_relaxed nit success")
    open(os.path.join(H, "relax_frank.txt"), "w").write("\n".join(L) + "\n")
