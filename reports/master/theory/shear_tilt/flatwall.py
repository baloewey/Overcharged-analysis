#!/usr/bin/env python3
"""Local limit: flat no-slip wall in simple shear v_x = s*y (y = distance from the wall), periodic in x with one site spacing a_s (one +1/2 per period).
Anchoring at the wall theta_w(x) = pi/2 + pi x / a_s (homeotropic at x = 0, rotates by pi per period -> one +1/2 screens it); far boundary y = H: zero
normal gradient (free, the zero mode can rotate). Same Q dynamics as shear_tilt.py (L = a = Gamma = 1, S0 0.876, lam_c = lam S0, one-way flow).
Output: tilt of the +1/2 (tail angle - outward normal (-y), positive downstream), its depth, its registry x/a_s, the far-field director angle alpha
(from the wall direction x), and the director angle at the defect depth away from the core."""
import argparse, json
import numpy as np
from scipy import ndimage
S0, h = 0.876, 0.5; c = 2 / S0 ** 2
def run(lam, s, a_s=20.0, H=40.0, T0=300.0, T1=1500.0, dt=0.04):
    nx, ny = int(round(a_s / h)), int(round(H / h)) + 2
    x = (np.arange(nx) + 0.5) * h; y = (np.arange(ny) - 1) * h          # row 0 = ghost wall row (y = -h/2 .. ), rows 1.. interior
    X, Y = np.meshgrid(x, y, indexing="ij")
    thw = np.pi / 2 + np.pi * x / a_s; qw1, qw2 = 0.5 * S0 * np.cos(2 * thw), 0.5 * S0 * np.sin(2 * thw)
    # initial: one +1/2 at (a_s/2... placed at x = 0 (homeotropic point), depth 2.5) with its image -> field theta = th_w-like; simple start: anchoring
    # director extended into the bulk decays to planar (theta = 0) with depth; defect nucleates by relaxation
    w = np.exp(-Y / 3.0); th0 = np.angle(np.exp(1j * 0)); th = np.where(True, w * (np.pi / 2 + np.pi * X / a_s) + (1 - w) * 0.0, 0)
    q1, q2 = 0.5 * S0 * np.cos(2 * th), 0.5 * S0 * np.sin(2 * th)
    def bc(q1, q2):
        q1[:, 0] = qw1; q2[:, 0] = qw2; q1[:, -1] = q1[:, -2]; q2[:, -1] = q2[:, -2]; return q1, q2
    def lap(f): return (np.roll(f, 1, 0) + np.roll(f, -1, 0) - 2 * f) / h ** 2 + (np.roll(f, 1, 1) + np.roll(f, -1, 1) - 2 * f) / h ** 2
    def gx(f): return (np.roll(f, -1, 0) - np.roll(f, 1, 0)) / (2 * h)
    def gy(f): return (np.roll(f, -1, 1) - np.roll(f, 1, 1)) / (2 * h)
    vx = s * np.maximum(Y, 0); Exy = 0.5 * s; Oxy = 0.5 * s                      # dvx/dy = s: E_xy = s/2, Omega_xy = (d_y v_x - d_x v_y)/2 = s/2
    def step(q1, q2, lam_c, flow):
        tr = 2 * (q1 ** 2 + q2 ** 2); H1 = lap(q1) + q1 - c * tr * q1; H2 = lap(q2) + q2 - c * tr * q2
        if flow: d1 = -vx * gx(q1) + 2 * Oxy * q2 + H1; d2 = -vx * gx(q2) + lam_c * Exy - 2 * Oxy * q1 + H2
        else: d1, d2 = H1, H2
        return bc(q1 + dt * d1, q2 + dt * d2)
    q1, q2 = bc(q1, q2)
    for _ in range(int(T0 / dt)): q1, q2 = step(q1, q2, 0, False)
    out0 = measure(q1, q2, X, Y, a_s)
    for _ in range(int(T1 / dt)): q1, q2 = step(q1, q2, lam * S0, True)
    return out0, measure(q1, q2, X, Y, a_s)
def measure(q1, q2, X, Y, a_s):
    S = 2 * np.hypot(q1, q2); inn = Y > 0.3
    lab, n = ndimage.label((S < 0.5 * S0) & inn, structure=np.ones((3, 3)))
    dvx = (np.roll(q1, -1, 0) - np.roll(q1, 1, 0)) / (2 * h) + (np.roll(q2, -1, 1) - np.roll(q2, 1, 1)) / (2 * h)
    dvy = (np.roll(q2, -1, 0) - np.roll(q2, 1, 0)) / (2 * h) - (np.roll(q1, -1, 1) - np.roll(q1, 1, 1)) / (2 * h)
    far = np.angle(np.mean((q1 + 1j * q2)[:, -12:-2])) / 2
    res = dict(ndef=int(n), far_alpha=float(np.degrees(far)))
    if n:
        i = np.argmax(ndimage.sum(S0 - S, lab, range(1, n + 1))) + 1; m = lab == i
        # periodic-safe centroid in x
        ang = np.angle(np.sum(np.exp(2j * np.pi * X[m] / a_s) * (S0 - S[m]))); cx = (ang % (2 * np.pi)) * a_s / (2 * np.pi); cy = float(np.sum(Y[m] * (S0 - S[m])) / np.sum(S0 - S[m]))
        dxp = (X - cx + a_s / 2) % a_s - a_s / 2; near = np.hypot(dxp, Y - cy) < 1.6
        tail = np.arctan2(dvy[near].mean(), dvx[near].mean()); tilt = np.angle(np.exp(1j * (tail + np.pi / 2)))
        row = np.argmin(abs(Y[0] - cy)); away = np.abs(dxp[:, row]) > a_s / 4
        loc = np.angle(np.mean((q1 + 1j * q2)[away, row])) / 2
        res.update(tilt=float(np.degrees(tilt)), depth=cy, registry=float(dxp[0, 0] * 0 + ((cx + a_s / 2) % a_s - a_s / 2) / a_s), alpha_at_depth=float(np.degrees(loc)))
    return res
if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--lam", type=float, default=2); ap.add_argument("--s", type=float, default=0.05); a = ap.parse_args()
    print(run(a.lam, a.s))
