#!/usr/bin/env python3
"""Flat-wall version of the analytic registry pinning (user, 2026-10-16): which untilted state does splay != bend prefer at a flat wall, at first order
in Delta = (K1 - K3)/2? (The circle calculation fluct/anisotropy_pinning gives tau = 0 = radial arm at the wall; the Landau-de Gennes flat wall relaxes
to tau = 180 deg.) One +1/2 per anchoring period a at depth d above the wall y = 0; wall anchoring theta_w(x) = pi/2 + pi x / a (strong); far field free.
Exact one-constant field: theta = theta0 + Im F(z) + h(x, y),  F = (1/2)[log sin(pi(z - z0)/a) - log sin(pi(z - conj z0)/a)] (row of +1/2 at z0 + j a with
opposite-sign images: uniform far field), h harmonic, periodic, decaying (e^{-|q|y}), cancelling the ripple of Im F on the wall so theta(x, 0) = theta_w(x).
Registry shift x0 -> uniform rotation by pi x0/a (zero mode for K1 = K3). First-order energy dF(x0) = (Delta/2) int (splay^2 - bend^2) w dA over one period,
smooth core cutoff w = 1 - exp(-r^2/c^2). Tilt = angle of the radial arm (2 theta_rest) from the outward normal (-y). K1 42, K3 23.1 -> Delta 9.45 kBT.
-> flat_pinning.txt"""
import os, numpy as np
H = os.path.dirname(os.path.abspath(__file__)); DELTA = (42.0 - 23.1) / 2


def setup(a, d, x0, hh=0.1, ymax=None):
    ymax = ymax or max(4 * a, 25.0); x = np.arange(0, a, hh) + hh / 2; y = np.arange(hh / 2, ymax, hh); X, Y = np.meshgrid(x, y, indexing="ij"); Z = X + 1j * Y
    z0 = x0 + 1j * d; k = np.pi / a
    F = 0.5 * (np.log(np.sin(k * (Z - z0))) - np.log(np.sin(k * (Z - np.conj(z0)))))
    dF = 0.5 * k * (1 / np.tan(k * (Z - z0)) - 1 / np.tan(k * (Z - np.conj(z0))))           # F'(z); theta_x = Im F', theta_y = Re F'
    # wall values of Im F (continuous branch), ripple relative to the linear winding
    xw = np.arange(0, a, hh / 4); zw = xw + 0j
    fw = np.unwrap(2 * np.imag(0.5 * (np.log(np.sin(k * (zw - z0))) - np.log(np.sin(k * (zw - np.conj(z0))))))) / 2
    slope = (fw[-1] - fw[0]) / (xw[-1] - xw[0]); sgn = np.sign(slope) if slope != 0 else 1
    rip = fw - sgn * k * xw; rip -= rip.mean()
    c = np.fft.rfft(rip) / len(rip); q = 2 * np.pi * np.arange(len(c)) / a
    h = np.zeros(X.shape); hx = np.zeros(X.shape); hy = np.zeros(X.shape)
    for m in range(1, min(len(c), 400)):
        e = np.exp(1j * q[m] * X - q[m] * Y); amp = 2 * c[m]
        h -= np.real(amp * e); hx -= np.real(1j * q[m] * amp * e); hy -= np.real(-q[m] * amp * e)
    base = np.imag(F) + h; thw = np.pi / 2 + sgn * k * X[:, 0]                                 # anchoring along the bottom row (y = hh/2 ~ wall)
    th0 = 0.5 * np.angle(np.mean(np.exp(2j * (thw - base[:, 0]))))                             # constant fixing theta = theta_w (mod pi) on the wall
    theta = th0 + base
    # make the branch of Im F consistent with fw on the wall: fix by matching theta at the bottom row to the anchoring modulo pi
    tx = np.imag(dF) + hx; ty = np.real(dF) + hy
    return X, Y, theta, tx, ty, z0, sgn, slope


def energy(a, d, x0, core=1.0):
    X, Y, th, tx, ty, z0, sgn, slope = setup(a, d, x0)
    dx = (X - z0.real + a / 2) % a - a / 2; r2 = dx ** 2 + (Y - z0.imag) ** 2; w = 1 - np.exp(-r2 / core ** 2)
    dens = -np.cos(2 * th) * (tx ** 2 - ty ** 2) - 2 * np.sin(2 * th) * tx * ty
    hh = X[1, 0] - X[0, 0]; return 0.5 * DELTA * np.sum(w * dens) * hh * hh, sgn, slope


def tilt(a, d, x0, eps=1e-3):
    X, Y, th, tx, ty, z0, sgn, slope = setup(a, d, x0)
    # regular part at the core: theta - (1/2) arg(z - z0), evaluated a small distance below the core along -y
    i = np.argmin(np.abs(X[:, 0] - x0)); j = np.argmin(np.abs(Y[0] - (d - 0.5)))
    zz = X[i, j] + 1j * Y[i, j]; rest = th[i, j] - 0.5 * np.angle(zz - z0)
    return np.degrees(np.angle(np.exp(1j * (2 * rest + np.pi / 2))))                      # radial arm 2 theta_rest vs outward normal (-pi/2)


if __name__ == "__main__":
    L = ["FLAT-WALL REGISTRY PINNING, first order in Delta = 9.45 kBT (K1 > K3). For each (a, d): dF(x0) over one period; minimum -> preferred tilt.",
         "a | d | core | winding check | dF range (kBT per defect) | tilt at the minimum (deg) | tilt at the maximum | stiffness d2F/dtau2 (kBT/rad^2)"]
    for a in (8.0, 12.0, 20.0, 30.0, 40.0):
        for d in (1.0, 1.5, 2.0):
            x0s = np.linspace(0, a, 25, endpoint=False); E = []; T = []
            for x0 in x0s:
                e, sgn, slope = energy(a, d, x0); E.append(e); T.append(tilt(a, d, x0))
            E = np.array(E); T = np.array(T); i0 = int(np.argmin(E)); i1 = int(np.argmax(E)); dtau = 2 * np.pi / 25
            d2 = (E[(i0 + 1) % 25] - 2 * E[i0] + E[(i0 - 1) % 25]) / dtau ** 2
            L.append(f"{a:4.0f} | {d:.1f} | 1.0 | slope*a/pi {slope * a / np.pi:+.2f} | {E.max() - E.min():7.2f} | {T[i0]:+7.1f} | {T[i1]:+7.1f} | {d2:8.2f}")
            print(L[-1], flush=True)
    open(os.path.join(H, "flat_pinning.txt"), "w").write("\n".join(L) + "\n")
