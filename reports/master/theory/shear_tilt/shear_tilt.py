#!/usr/bin/env python3
"""Does the no-slip shear set the tilt of the wall +1/2? Passive 2D Beris-Edwards Q-tensor dynamics in the anchored disk with an IMPOSED flow (one-way, no backflow).
  dQ/dt + v.grad Q = lam_c E + (Omega Q - Q Omega) + Gamma H,   H = L lap Q + a Q - c tr(Q^2) Q,   Q = [[q1, q2], [q2, -q1]],  Q = S (nn - I/2)
  Omega_ij = (d_j v_i - d_i v_j)/2 (director co-rotates with the fluid), E = symmetric strain rate. In simple shear the director obeys
  dtheta/dt = -gdot/2 + (lam_c / (2 S)) gdot cos 2theta, so lam_c = lam * S_eq gives cos 2theta_L = 1/lam as in MPCD (tumble parameter lam).
  L = a = Gamma = 1 (core ~ 1 cell), S_eq = 0.876, h = 0.5, R = 31.83, k = 6; anchoring theta = k phi held on the ghost ring r >= R.
Flow v_theta(r) = U * shape(r): shape from the measured passive swirl run (profile_swirl_threshold_passive_A0.01_s2.txt, normalised to max 1), variants:
  noslip  : measured shape (shear everywhere, no-slip at the wall);
  slip    : Omega(r) = v/r frozen at its value at r_cut = R - 6 for r > r_cut  -> rigid rotation of the wall band (no shear there, wall slips);
  wallonly: rigid rotation for r < r_cut, measured shape beyond                -> shear only in the wall band (no-slip kept).
Start: analytic ring of N wall +1/2 (12 = passive, or 10 + bulk pair = ejected), relaxed T0 without flow, then flow for T1.
Output (stdout + JSON): per sampled time, wall +1/2 count, tilt coherence, mean tilt signed with the flow (tail angle - outward normal, like defects.dat), ring rotation."""
import json, os, sys, argparse
import numpy as np
from scipy import ndimage
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(os.path.dirname(H), "ejection"))
from tilted_structures import fields  # noqa
R, k, h, S0 = 31.830989, 6, 0.5, 0.876; Lc, a, G = 1.0, 1.0, 1.0; c = 2 * a / S0 ** 2
x = np.arange(-R - 2, R + 2 + h / 2, h); X, Y = np.meshgrid(x, x, indexing="ij"); Rr = np.hypot(X, Y); PH = np.arctan2(Y, X)
inside = Rr < R; ghost = (~inside) & (Rr < R + 2)
qb1, qb2 = 0.5 * S0 * np.cos(2 * k * PH), 0.5 * S0 * np.sin(2 * k * PH)
def shape(variant):
    d = np.loadtxt(os.path.join(H, "profile_swirl_threshold_passive_A0.01_s2.txt")); r, v = np.r_[0, d[:, 0], R], np.r_[0, d[:, 1], 0]; v = v / v.max()
    rr = np.linspace(0, R, 2000); vv = np.interp(rr, r, v); rc = R - 6; Om = vv / np.maximum(rr, 1e-9); Oc = np.interp(rc, rr, Om)
    if variant == "slip": Om = np.where(rr > rc, Oc, Om)
    elif variant == "wallonly": Om = np.where(rr < rc, Oc, Om)
    return rr, Om
def setup_flow(U, variant):
    rr, Om = shape(variant); W = U * np.interp(Rr, rr, Om)          # angular velocity field Omega(r)
    W = np.where(inside, W, 0.0); vx, vy = -W * Y, W * X
    d = lambda f, ax: np.gradient(f, h, axis=ax)
    dvx_dx, dvx_dy, dvy_dx, dvy_dy = d(vx, 0), d(vx, 1), d(vy, 0), d(vy, 1)
    Exx, Exy = 0.5 * (dvx_dx - dvy_dy), 0.5 * (dvx_dy + dvy_dx); Oxy = 0.5 * (dvx_dy - dvy_dx)   # Omega_xy = (d_y v_x - d_x v_y)/2
    return vx, vy, Exx, Exy, Oxy
def lap(f): return (np.roll(f, 1, 0) + np.roll(f, -1, 0) + np.roll(f, 1, 1) + np.roll(f, -1, 1) - 4 * f) / h ** 2
def grad(f): return (np.roll(f, -1, 0) - np.roll(f, 1, 0)) / (2 * h), (np.roll(f, -1, 1) - np.roll(f, 1, 1)) / (2 * h)
def step(q1, q2, dt, lam_c, fl):
    trQ2 = 2 * (q1 ** 2 + q2 ** 2); H1 = Lc * lap(q1) + a * q1 - c * trQ2 * q1; H2 = Lc * lap(q2) + a * q2 - c * trQ2 * q2
    if fl is None: d1, d2 = G * H1, G * H2
    else:
        vx, vy, Exx, Exy, Oxy = fl; g1x, g1y = grad(q1); g2x, g2y = grad(q2)
        # (Omega Q - Q Omega) with Omega = [[0, Oxy], [-Oxy, 0]]: xx-component = 2 Oxy q2, xy-component = -2 Oxy q1
        d1 = -(vx * g1x + vy * g1y) + lam_c * Exx + 2 * Oxy * q2 + G * H1
        d2 = -(vx * g2x + vy * g2y) + lam_c * Exy - 2 * Oxy * q1 + G * H2
    q1 = q1 + dt * d1; q2 = q2 + dt * d2
    q1[ghost] = qb1[ghost]; q2[ghost] = qb2[ghost]; q1[~(inside | ghost)] = 0; q2[~(inside | ghost)] = 0
    return q1, q2
def defects(q1, q2):
    S = 2 * np.sqrt(q1 ** 2 + q2 ** 2); lab, n = ndimage.label((S < 0.5 * S0) & inside)
    out = []
    if n == 0: return out
    cms = ndimage.center_of_mass(S0 - S, lab, range(1, n + 1)); g1x, g1y = grad(q1); g2x, g2y = grad(q2); dvx, dvy = g1x + g2y, g2x - g1y   # div Q
    for cx, cy in cms:
        px, py = x[0] + cx * h, x[0] + cy * h; th = np.linspace(0, 2 * np.pi, 48, endpoint=False); rl = 2.5
        ix, iy = (px + rl * np.cos(th) - x[0]) / h, (py + rl * np.sin(th) - x[0]) / h
        a1 = ndimage.map_coordinates(q1, [ix, iy], order=1); a2 = ndimage.map_coordinates(q2, [ix, iy], order=1)
        w = np.angle(np.exp(1j * np.diff(np.r_[np.arctan2(a2, a1), np.arctan2(a2[0], a1[0])]))).sum() / (4 * np.pi)
        m = np.hypot(X - px, Y - py) < 1.6; ang = np.arctan2(dvy[m].mean(), dvx[m].mean())
        out.append((px, py, round(w * 2) / 2, ang))
    return out
def observe(q1, q2, sgn):
    D = defects(q1, q2); wall = [(px, py, ang) for px, py, q, ang in D if q > 0 and np.hypot(px, py) > R - 6]
    if not wall: return dict(nwall=0, nbulk_p=sum(1 for d in D if d[2] > 0), nbulk_m=sum(1 for d in D if d[2] < 0))
    phi = np.array([np.arctan2(py, px) for px, py, _ in wall]); tau = np.angle(np.exp(1j * (np.array([w[2] for w in wall]) - phi))) * sgn
    z = np.mean(np.exp(1j * tau))
    return dict(nwall=len(wall), nbulk_p=sum(1 for d in D if d[2] > 0) - len(wall), nbulk_m=sum(1 for d in D if d[2] < 0),
                coh=float(abs(z)), tilt=float(np.degrees(np.angle(z))), tilts=[float(t) for t in np.degrees(tau)], phis=[float(p) for p in phi])
def run(lam, U, variant, ring, T0=300.0, T1=3000.0, dt=0.04, every=100.0):
    zs = list((R - 2.6) * np.exp(2j * np.pi * (np.arange(ring) + 0.5) / ring)); qs = [0.5] * ring
    if ring == 2 * k - 2: zs = list((R - 2.6) * np.exp(2j * np.pi * np.arange(ring) / ring)) + [0.45j * R, -0.45j * R]; qs = [0.5] * (ring + 2)
    th2, _, _ = fields(X, Y, zs, qs, R); q1 = np.where(inside, 0.5 * S0 * np.cos(th2), qb1); q2 = np.where(inside, 0.5 * S0 * np.sin(th2), qb2)
    for _ in range(int(T0 / dt)): q1, q2 = step(q1, q2, dt, 0, None)
    obs = [dict(t=0.0, **observe(q1, q2, np.sign(U) or 1))]; fl = setup_flow(U, variant); lam_c = lam * S0; t = 0.0
    while t < T1 - 1e-9:
        for _ in range(int(every / dt)): q1, q2 = step(q1, q2, dt, lam_c, fl)
        t += every; obs.append(dict(t=t, **observe(q1, q2, np.sign(U) or 1)))
    return obs
if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--lam", type=float, default=2); ap.add_argument("--U", type=float, default=0.03)
    ap.add_argument("--variant", default="noslip"); ap.add_argument("--ring", type=int, default=12); ap.add_argument("--T1", type=float, default=3000)
    ap.add_argument("--T0", type=float, default=300); ap.add_argument("--out", default=None)
    A = ap.parse_args(); obs = run(A.lam, A.U, A.variant, A.ring, T0=A.T0, T1=A.T1)
    for o in obs: print(o["t"], o["nwall"], o.get("nbulk_p"), o.get("nbulk_m"), round(o.get("coh", np.nan), 2), round(o.get("tilt", np.nan), 1), flush=True)
    if A.out: json.dump(dict(args=vars(A), obs=obs), open(A.out, "w"))
