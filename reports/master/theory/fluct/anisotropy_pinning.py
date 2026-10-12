#!/usr/bin/env python3
"""Registry pinning from splay-bend anisotropy (first order in Delta = (K1 - K3)/2).
One-constant solution in the disk (Dirichlet anchoring theta = k phi + theta0 on |z| = R), exact for any defect positions:
  e^{2 i theta} = e^{2 i theta0} prod_i [ (z - z_i)/|z - z_i| * (1 - z conj(z_i)/R^2)/|1 - z conj(z_i)/R^2| ]^{2 q_i}
  grad theta = sum_i q_i [ (Im, Re)(1/(z - z_i)) + (Im, Re)(-(conj(z_i)/R^2)/(1 - z conj(z_i)/R^2)) ]
Its one-constant energy is independent of a rigid rotation alpha of the lattice (zero mode). With K1 = K + Delta, K3 = K - Delta the
first-order correction is dF = (Delta/2) int (splay^2 - bend^2) dA, splay^2 - bend^2 = -cos(2 theta)(theta_x^2 - theta_y^2) - 2 sin(2 theta) theta_x theta_y,
integrated over the disk with a smooth core cutoff 1 - exp(-|z - z_i|^2/a^2) per defect. dF(alpha) is the pinning potential of the lattice registry.
Configurations: passive ring N = 2k (q = 1/2) and active ring N = 2k - 2 plus a central +1, at the depth-law radius (+0.61 offset).
Also: the head tilt of each wall +1/2 (from the regular part of theta at its position) to locate the minimum relative to 'heads at the wall'.
K1 = 42.0, K3 = 23.1 kBT (Sec. frank) -> Delta = 9.45 kBT. Predicted collective rms: R / sqrt(dF''(alpha_min)) (kBT = 1).
Writes anisotropy_pinning.txt and anisotropy_pinning.png."""
import json, os
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
H = os.path.dirname(os.path.abspath(__file__))
DELTA = (42.0 - 23.1) / 2
def field(zs, qs, R, X, Y, th0=0.0):
    z = X + 1j * Y; e2 = np.exp(2j * th0) * np.ones_like(z); gx = np.zeros(z.shape); gy = np.zeros(z.shape)
    for zi, q in zip(zs, qs):
        a = z - zi; b = 1 - z * np.conj(zi) / R ** 2
        e2 *= ((a / np.abs(a)) * (b / np.abs(b))) ** (2 * q)
        w = 1 / a + (-(np.conj(zi) / R ** 2) / b)
        gx += q * w.imag * (-1) * (-1); gy += q * w.real
        # d/dx arg f = Im(f'/f), d/dy arg f = Re(f'/f)
    return e2, gx, gy
def dF(zs, qs, R, h=0.2, a=1.0):
    x = np.arange(-R + h / 2, R, h); X, Y = np.meshgrid(x, x); r = np.hypot(X, Y)
    m = r < R
    for zi in zs: m &= np.abs(X + 1j * Y - zi) > 0.05                    # avoid the singular point itself
    w = np.ones(m.sum())
    for zi in zs: w *= 1 - np.exp(-np.abs(X[m] + 1j * Y[m] - zi) ** 2 / a ** 2)   # smooth core cutoff (radius a)
    e2, gx, gy = field(zs, qs, R, X[m], Y[m])
    c2, s2 = e2.real, e2.imag
    dens = -c2 * (gx ** 2 - gy ** 2) - 2 * s2 * gx * gy
    return 0.5 * DELTA * (w * dens).sum() * h * h
def tilts(zs, qs, R, eps=1e-3):
    """head direction of each +1/2: theta_rest = theta minus its own (1/2) arg(z - z_j); head = 2 theta_rest; tilt = head - outward normal."""
    out = []
    for j, (zj, qj) in enumerate(zip(zs, qs)):
        if qj != 0.5 or abs(zj) < 1: continue
        p = zj + eps * zj / abs(zj)
        e2, _, _ = field(zs, qs, R, np.array([p.real]), np.array([p.imag]))
        own = ((p - zj) / abs(p - zj))                        # e^{i arg(p - z_j)} for q = 1/2 -> exponent 2q = 1
        e2r = e2[0] / own
        head = np.angle(e2r)                                  # 2 theta_rest
        out.append(np.angle(np.exp(1j * (head - np.angle(zj)))))
    return np.array(out)
L = []; fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
R = 31.830989
for kind, ci in (("passive", 0), ("active", 1)):
    for k in (4, 6, 8):
        N = 2 * k if kind == "passive" else 2 * k - 2
        Qb = k - N / 2                                           # bulk charge at the centre (0 or 1)
        rs = R * ((N - 1 + 4 * Qb) / (3 * N - 1 + 4 * Qb)) ** (1 / (2 * N)) - 0.61
        per = 2 * np.pi / (2 * k - 2)
        als = np.linspace(0, per, 25)
        for a in (0.5, 1.0, 2.0):
            E = []
            for al in als:
                zs = list(rs * np.exp(1j * (2 * np.pi * np.arange(N) / N + al))); qs = [0.5] * N
                if Qb: zs.append(0j); qs.append(1.0)
                E.append(dF(zs, qs, R, a=a))
            E = np.array(E); E -= E.min(); i0 = E.argmin()
            # curvature at the minimum (periodic finite difference)
            d2 = (E[(i0 + 1) % 24] - 2 * E[i0] + E[(i0 - 1) % 24]) / (als[1] - als[0]) ** 2 if 0 < i0 < 24 else (E[1] - 2 * E[0] + E[23]) / (als[1] - als[0]) ** 2
            zs = list(rs * np.exp(1j * (2 * np.pi * np.arange(N) / N + als[i0]))); qs = [0.5] * N
            if Qb: zs.append(0j); qs.append(1.0)
            tl = tilts(zs, qs, R)
            sig = R / np.sqrt(d2) if d2 > 0 else np.inf
            L.append(f"{kind:7s} k={k} N={N:2d} core a={a}: pinning barrier {E.max():.2f} kBT, minimum at (2k-2)*alpha = {(2 * k - 2) * als[i0]:+.2f} rad, "
                     f"curvature {d2:.0f} kBT/rad^2 -> collective rms {sig:.2f} cells; head tilts at the minimum (deg): "
                     + " ".join(f"{np.degrees(t):+.0f}" for t in tl))
            print(L[-1])
            if a == 1.0: ax[ci].plot((2 * k - 2) * als, E, "o-", ms=3, label=f"k = {k}, N = {N}")
    ax[ci].set_xlabel(r"registry $(2k-2)\,\alpha$ (rad)"); ax[ci].set_ylabel(r"$\Delta F$ ($k_BT$)"); ax[ci].legend(fontsize=8)
    ax[ci].set_title(f"({'ab'[ci]}) {kind} ring, first order in $K_1-K_3$, core 1 cell (smooth)", fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(H, "anisotropy_pinning.png"), dpi=150)
open(os.path.join(H, "anisotropy_pinning.txt"), "w").write("\n".join(L) + "\n")
