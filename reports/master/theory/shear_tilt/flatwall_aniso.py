#!/usr/bin/env python3
"""Flat no-slip wall in simple shear with splay != bend. Same as flatwall.py plus the cubic Landau-de Gennes term
  F3 = (L3/2) int Q_kl d_k Q_ab d_l Q_ab   ->  for Q = S(nn - I/2):  F3 = (L3 S^3 / 2) int [ (n.grad theta)^2 - (n_perp.grad theta)^2 ]
so K3 (bend) = 2 L1 S^2 + L3 S^3, K1 (splay) = 2 L1 S^2 - L3 S^3. With Q = [[q1, q2], [q2, -q1]], G_ij = d_i q1 d_j q1 + d_i q2 d_j q2:
  H_q1 = L1 lap q1 + bulk - (L3/2)(G_xx - G_yy) + L3 div(Q grad q1),   H_q2 = L1 lap q2 + bulk - L3 G_xy + L3 div(Q grad q2)   (H = -1/2 dF/dq).
Two ghost rows at the wall (anchoring) and two at the top (zero gradient). L1 is rescaled so that (K1 + K3)/2 = 2 S^2 (same mean K as the one-constant runs).
--check: finite-difference test that H is the functional derivative; measured K1, K3 from small splay/bend modulations."""
import argparse, json
import numpy as np
from scipy import ndimage
S0, h = 0.876, 0.5; c = 2 / S0 ** 2
def ops(nx, ny):
    def dx(f): return (np.roll(f, -1, 0) - np.roll(f, 1, 0)) / (2 * h)
    def dy(f):
        g = np.zeros_like(f); g[:, 1:-1] = (f[:, 2:] - f[:, :-2]) / (2 * h); return g
    def lap(f):
        g = np.zeros_like(f); g[:, 1:-1] = (np.roll(f, 1, 0) + np.roll(f, -1, 0) - 2 * f)[:, 1:-1] / h ** 2 + (f[:, 2:] + f[:, :-2] - 2 * f[:, 1:-1]) / h ** 2; return g
    return dx, dy, lap
def Hfield(q1, q2, L1, L3, ops_):
    dx, dy, lap = ops_; tr = 2 * (q1 ** 2 + q2 ** 2)
    g1x, g1y, g2x, g2y = dx(q1), dy(q1), dx(q2), dy(q2); Gxx = g1x ** 2 + g2x ** 2; Gyy = g1y ** 2 + g2y ** 2; Gxy = g1x * g1y + g2x * g2y
    # (Q grad q)_l = sum_k Q_kl d_k q :  x-comp = q1 d_x q + q2 d_y q ;  y-comp = q2 d_x q - q1 d_y q
    d1 = dx(q1 * g1x + q2 * g1y) + dy(q2 * g1x - q1 * g1y); d2 = dx(q1 * g2x + q2 * g2y) + dy(q2 * g2x - q1 * g2y)
    H1 = L1 * lap(q1) + q1 - c * tr * q1 - 0.5 * L3 * (Gxx - Gyy) + L3 * d1
    H2 = L1 * lap(q2) + q2 - c * tr * q2 - L3 * Gxy + L3 * d2
    return H1, H2
def energy(q1, q2, L1, L3, ops_):
    dx, dy, lap = ops_; g1x, g1y, g2x, g2y = dx(q1), dy(q1), dx(q2), dy(q2); m = slice(2, -2)
    f = L1 * (g1x ** 2 + g1y ** 2 + g2x ** 2 + g2y ** 2) + L3 * (q1 * (g1x ** 2 + g2x ** 2 - g1y ** 2 - g2y ** 2) + 2 * q2 * (g1x * g1y + g2x * g2y)) \
        - (q1 ** 2 + q2 ** 2) + c * (q1 ** 2 + q2 ** 2) ** 2
    return np.sum(f[:, m]) * h * h
def params(ratio):          # K1/K3 = ratio at mean K = 2 S^2 (L1 = 1 for ratio 1)
    Km = 2 * S0 ** 2; K1 = 2 * Km * ratio / (1 + ratio); K3 = 2 * Km / (1 + ratio); L3 = (K3 - K1) / (2 * S0 ** 3); L1 = (K1 + K3) / (4 * S0 ** 2); return L1, L3, K1, K3
def run(lam, s, ratio=1.0, a_s=20.0, H=30.0, T0=1500.0, T1=3000.0, dt=0.02):
    L1, L3, K1, K3 = params(ratio); nx, ny = int(round(a_s / h)), int(round(H / h)) + 4
    x = (np.arange(nx) + 0.5) * h; y = (np.arange(ny) - 2) * h; X, Y = np.meshgrid(x, y, indexing="ij"); O = ops(nx, ny)
    thw = np.pi / 2 + np.pi * x / a_s; qw1, qw2 = 0.5 * S0 * np.cos(2 * thw), 0.5 * S0 * np.sin(2 * thw)
    w = np.exp(-np.maximum(Y, 0) / 3.0); th = w * (np.pi / 2 + np.pi * X / a_s); q1, q2 = 0.5 * S0 * np.cos(2 * th), 0.5 * S0 * np.sin(2 * th)
    def bc(q1, q2):
        q1[:, :2] = qw1[:, None]; q2[:, :2] = qw2[:, None]; q1[:, -2:] = q1[:, -3:-2]; q2[:, -2:] = q2[:, -3:-2]; return q1, q2
    vx = s * np.maximum(Y, 0); Exy = Oxy = 0.5 * s; dx, dy, lap = O
    def step(q1, q2, lam_c, flow):
        H1, H2 = Hfield(q1, q2, L1, L3, O)
        if flow: d1 = -vx * dx(q1) + 2 * Oxy * q2 + H1; d2 = -vx * dx(q2) + lam_c * Exy - 2 * Oxy * q1 + H2
        else: d1, d2 = H1, H2
        return bc(q1 + dt * d1, q2 + dt * d2)
    q1, q2 = bc(q1, q2)
    for _ in range(int(T0 / dt)): q1, q2 = step(q1, q2, 0, False)
    import flatwall as F
    o0 = F.measure(q1[:, 1:-1], q2[:, 1:-1], X[:, 1:-1] , Y[:, 1:-1], a_s)
    for _ in range(int(T1 / dt)): q1, q2 = step(q1, q2, lam * S0, True)
    return o0, F.measure(q1[:, 1:-1], q2[:, 1:-1], X[:, 1:-1], Y[:, 1:-1], a_s)
def check():
    rng = np.random.default_rng(0); nx, ny = 24, 28; O = ops(nx, ny); L1, L3 = 0.8, -0.5
    q1 = 0.3 + 0.05 * rng.standard_normal((nx, ny)); q2 = 0.1 + 0.05 * rng.standard_normal((nx, ny))
    H1, H2 = Hfield(q1, q2, L1, L3, O); i, j = 7, 13; e = 1e-6
    for comp, Hc in ((0, H1), (1, H2)):
        a, b = q1.copy(), q2.copy(); (a if comp == 0 else b)[i, j] += e; Ep = energy(a, b, L1, L3, O)
        a, b = q1.copy(), q2.copy(); (a if comp == 0 else b)[i, j] -= e; Em = energy(a, b, L1, L3, O)
        print(f"component q{comp + 1}: -1/2 dF/dq (finite diff, per cell area) {-(Ep - Em) / (2 * e) / (2 * h * h):+.6f}   H {Hc[i, j]:+.6f}")
    for ratio in (1.0, 1.82):          # measured K from small modulations theta = eps sin(k x) / sin(k y) around a director along x
        L1, L3, K1t, K3t = params(ratio); n = 64; xx = (np.arange(n) + 0.5) * h; X, Y = np.meshgrid(xx, xx, indexing="ij"); O2 = ops(n, n); kk = 2 * np.pi / (n * h); eps = 0.02
        out = []
        for th in (eps * np.sin(kk * X), eps * np.sin(kk * Y)):
            q1 = 0.5 * S0 * np.cos(2 * th); q2 = 0.5 * S0 * np.sin(2 * th); E = energy(q1, q2, L1, L3, O2) - energy(0.5 * S0 + 0 * q1, 0 * q2, L1, L3, O2)
            area = (n * h) * ((n - 4) * h); out.append(E / (0.5 * eps ** 2 * kk ** 2 * 0.5 * area))
        print(f"ratio {ratio}: bend (theta(x), director || x) K = {out[0]:.3f} (target K3 {K3t:.3f});  splay (theta(y)) K = {out[1]:.3f} (target K1 {K1t:.3f})")
if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--check", action="store_true"); a = ap.parse_args()
    if a.check: check()
