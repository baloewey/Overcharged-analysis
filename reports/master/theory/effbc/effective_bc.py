#!/usr/bin/env python3
"""Effective boundary condition seen by the bulk through the wall layer (user's question, 2026-10-07).
Exact one-constant director (disk, Dirichlet anchoring theta = k phi as in lc.c with O = 0):
   theta(z) = sum_i q_i [arg(z - z_i) + arg(1 - z conj(z_i)/R^2)] (+ const fixed by the anchoring; exact for any positions).
Layer: N = 2k - 2 wall +1/2 at the depth-law radius, at registry alpha (sites phi_j = 2 pi j/N + alpha; tilt tau = (2k-2) alpha;
alpha = 0 is the measured registry, tails at the wall). Bulk: +1 at the centre (or two +1/2 at +-a).
On circles of radius r inside the layer: theta - phi = c + modulation. Reported: the mean offset c (pi/2 = planar/tangential,
0 = homeotropic/radial, else spiral) and the rms modulation vs depth delta = r_ring - r in units of the spacing l = 2 pi r_ring / N,
compared with exp(-2 pi delta / l). Writes effective_bc.txt and effective_bc.png."""
import os
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
H = os.path.dirname(os.path.abspath(__file__)); L = []
def say(s): print(s); L.append(s)
def theta(z, zs, qs, R):
    th = np.zeros(z.shape)
    for zi, q in zip(zs, qs):
        th += q * (np.angle(z - zi) + np.angle(1 - z * np.conj(zi) / R ** 2))
    return th
R = 31.830989
fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
for k in (3, 6, 10):
    N = 2 * k - 2; rs = R * ((N - 1 + 4) / (3 * N - 1 + 4)) ** (1 / (2 * N)) - 0.61; ell = 2 * np.pi * rs / N
    for tau_deg in (0, 30, 63, 90):
        al = np.radians(tau_deg) / (2 * k - 2)
        zs = list(rs * np.exp(1j * (2 * np.pi * np.arange(N) / N + al))) + [0j]; qs = [0.5] * N + [1.0]
        rows = []
        for dl in np.linspace(0.05, 2.0, 40) * ell:
            r = rs - dl
            if r < 2: break
            ph = np.linspace(0, 2 * np.pi, 720, endpoint=False); z = r * np.exp(1j * ph)
            d = np.angle(np.exp(2j * (theta(z, zs, qs, R) - ph))) / 2          # theta - phi mod pi, in (-pi/2, pi/2]
            c = np.angle(np.mean(np.exp(2j * d))) / 2; mod = np.sqrt(np.mean(np.angle(np.exp(2j * (d - c))) ** 2)) / 2
            rows.append((dl / ell, c, mod))
        rows = np.array(rows)
        cmid = rows[len(rows) // 2, 1]
        # depth where the modulation falls below 0.1 rad (about 6 deg)
        dsm = rows[np.argmax(rows[:, 2] < 0.1), 0] if (rows[:, 2] < 0.1).any() else np.nan
        say(f"k={k:2d} N={N:2d} spacing l={ell:5.1f}  tilt {tau_deg:3d} deg: effective offset c = {np.degrees(cmid) % 180:6.1f} deg "
            f"(90 = planar), modulation < 0.1 rad beyond depth {dsm:.2f} l ({dsm * ell:.1f} cells)")
        if tau_deg == 0:
            ax[0].semilogy(rows[:, 0], rows[:, 2], "o-", ms=3, label=f"k = {k} (N = {N})")
        if k == 6: ax[1].plot(rows[:, 0], np.degrees(rows[:, 1]) % 180, "-", label=f"tilt {tau_deg} deg")
xx = np.linspace(0.05, 2, 50); ax[0].semilogy(xx, 0.5 * np.exp(-2 * np.pi * xx), "k--", label=r"$\propto e^{-2\pi\delta/\ell}$")
ax[0].set_xlabel(r"depth inside the ring, $\delta/\ell$"); ax[0].set_ylabel("rms of $\\theta-\\phi-c$ (rad)"); ax[0].legend(fontsize=7)
ax[0].set_title("(a) Site pattern fades with depth", fontsize=10)
ax[1].set_xlabel(r"depth inside the ring, $\delta/\ell$"); ax[1].set_ylabel(r"offset $c=\langle\theta-\phi\rangle$ (deg)"); ax[1].legend(fontsize=7)
ax[1].set_title("(b) Effective anchoring offset, k = 6", fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(H, "effective_bc.png"), dpi=150)
open(os.path.join(H, "effective_bc.txt"), "w").write("\n".join(L) + "\n")
