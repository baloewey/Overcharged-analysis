#!/usr/bin/env python3
"""Active continuum model in ellipses: Q-tensor (one-constant LdG, S0 0.876, L = a = Gamma = 1) + Stokes with backflow and extensile activity (as active_bf.py),
in an ellipse of semi-axes A = R_eq sqrt(eps), B = R_eq / sqrt(eps) (same area as the circle R_eq = 31.83), arclength-uniform anchoring (MPCD eMode 2):
theta_w(s) = 2 pi k s / P from the tip (A, 0) (homeotropic there), held on a ghost band of 2 cells outside the wall; no-slip Stokes (clamped streamfunction).
Optional additive noise on Q (amplitude noise * sqrt(dt) per step) so that the dynamics can switch between states.
Observables every 100 t.u.: wall +1/2 (distance to the wall < 4), bulk +1/2 / -1/2, tilt coherence and mean tilt (tail - local outward normal), circulation
<v_theta>/v_rms, two-vortex flag (the halves x < 0 and x > 0 rotate oppositely about their centroids).
Usage: ellipse_active.py --eps 2 --zeta 0.3 [--lam 2 --k 6 --T1 8000 --noise 0.005 --seed 1 --out file.json]"""
import argparse, json, time
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spla
from scipy import ndimage
from scipy.spatial import cKDTree
S0, h = 0.876, 0.5; Lc = a_ = G = 1.0; c_ = 2 / S0 ** 2; REQ = 31.830989
class Geo:
    def __init__(self, eps, k):
        self.A, self.B = REQ * np.sqrt(eps), REQ / np.sqrt(eps)
        x = np.arange(-self.A - 4, self.A + 4 + h / 2, h); y = np.arange(-self.B - 4, self.B + 4 + h / 2, h); self.x, self.y = x, y
        self.X, self.Y = np.meshgrid(x, y, indexing="ij")
        t = np.linspace(0, 2 * np.pi, 4001)[:-1]; W = np.column_stack([self.A * np.cos(t), self.B * np.sin(t)]); seg = np.hypot(*np.diff(np.vstack([W, W[:1]]), axis=0).T)
        s = np.concatenate([[0], np.cumsum(seg)[:-1]]); self.P = seg.sum(); self.W = W; self.tree = cKDTree(W)
        nrm = np.column_stack([np.cos(t) / self.A, np.sin(t) / self.B]); nrm /= np.hypot(*nrm.T)[:, None]; self.nang = np.arctan2(nrm[:, 1], nrm[:, 0])
        pts = np.column_stack([self.X.ravel(), self.Y.ravel()]); d, j = self.tree.query(pts); d = d.reshape(self.X.shape); j = j.reshape(self.X.shape)
        self.inside = (self.X / self.A) ** 2 + (self.Y / self.B) ** 2 < 1; self.ghost = (~self.inside) & (d < 2.0); self.depth = np.where(self.inside, d, -d)
        thw = 2 * np.pi * k * s[j] / self.P; self.qb1 = np.where(self.ghost, 0.5 * S0 * np.cos(2 * thw), 0); self.qb2 = np.where(self.ghost, 0.5 * S0 * np.sin(2 * thw), 0)
        self.thw = thw
        n0, n1 = self.X.shape; I0, I1 = sp.identity(n0), sp.identity(n1); D0 = sp.diags([1, -2, 1], [-1, 0, 1], shape=(n0, n0)) / h ** 2; D1 = sp.diags([1, -2, 1], [-1, 0, 1], shape=(n1, n1)) / h ** 2
        LAP = (sp.kron(D0, I1) + sp.kron(I0, D1)).tocsr(); self.idx = np.flatnonzero(self.inside.ravel())
        self.solve = spla.factorized((LAP @ LAP).tocsr()[self.idx][:, self.idx].tocsc())
def lap(f): return (np.roll(f, 1, 0) + np.roll(f, -1, 0) + np.roll(f, 1, 1) + np.roll(f, -1, 1) - 4 * f) / h ** 2
def dx(f): return (np.roll(f, -1, 0) - np.roll(f, 1, 0)) / (2 * h)
def dy(f): return (np.roll(f, -1, 1) - np.roll(f, 1, 1)) / (2 * h)
def molecular(q1, q2):
    tr = 2 * (q1 ** 2 + q2 ** 2); return Lc * lap(q1) + a_ * q1 - c_ * tr * q1, Lc * lap(q2) + a_ * q2 - c_ * tr * q2
def flow(g, q1, q2, zeta, lam_c):
    zeta = zeta * g.zmask if hasattr(g, 'zmask') else zeta
    if getattr(g, 'actmode', 'all') == 'core':                     # activity only at wall-defect cores: S < 0.6 within 4 cells of the wall, dilated by 1 cell
        core = ndimage.binary_dilation((2 * np.hypot(q1, q2) < 0.6) & g.inside & (g.depth < 4), iterations=2); zeta = zeta * core
    h1, h2 = molecular(q1, q2); h1 = np.where(g.inside, h1, 0); h2 = np.where(g.inside, h2, 0); asym = 2 * (q1 * h2 - q2 * h1)
    sxx = -zeta * q1 - lam_c * h1; sxy = -zeta * q2 - lam_c * h2 + asym; syx = -zeta * q2 - lam_c * h2 - asym; syy = zeta * q1 + lam_c * h1
    fx = dx(sxx) + dy(sxy); fy = dx(syx) + dy(syy)
    if getattr(g, 'swirl', 0.0): fx = fx + g.swirl_fx; fy = fy + g.swirl_fy
    psi = np.zeros(g.X.size); psi[g.idx] = g.solve((dx(fy) - dy(fx)).ravel()[g.idx]); psi = psi.reshape(g.X.shape)
    vx, vy = dy(psi), -dx(psi); return vx, vy, 0.5 * (dx(vx) - dy(vy)), 0.5 * (dy(vx) + dx(vy)), 0.5 * (dy(vx) - dx(vy))
def step(g, q1, q2, dt, lam_c, fl, noise, rng):
    H1, H2 = molecular(q1, q2); vx, vy, Exx, Exy, Oxy = fl
    d1 = -(vx * dx(q1) + vy * dy(q1)) + lam_c * Exx + 2 * Oxy * q2 + G * H1; d2 = -(vx * dx(q2) + vy * dy(q2)) + lam_c * Exy - 2 * Oxy * q1 + G * H2
    q1 = q1 + dt * d1; q2 = q2 + dt * d2
    if noise: q1 = q1 + noise * np.sqrt(dt) * rng.standard_normal(q1.shape) * g.inside; q2 = q2 + noise * np.sqrt(dt) * rng.standard_normal(q2.shape) * g.inside
    q1 = np.where(g.ghost, g.qb1, np.where(g.inside, q1, 0)); q2 = np.where(g.ghost, g.qb2, np.where(g.inside, q2, 0)); return q1, q2
def observe(g, q1, q2, vx, vy):
    S = 2 * np.hypot(q1, q2); lab, n = ndimage.label((S < 0.5 * S0) & g.inside); o = dict(nwall=0, nbp=0, nbm=0)
    if n:
        cms = ndimage.center_of_mass(S0 - S, lab, range(1, n + 1)); dvx, dvy = dx(q1) + dy(q2), dx(q2) - dy(q1); taus = []
        for cx, cy in cms:
            px, py = g.x[0] + cx * h, g.y[0] + cy * h; th = np.linspace(0, 2 * np.pi, 48, endpoint=False)
            a1 = ndimage.map_coordinates(q1, [(px + 2.5 * np.cos(th) - g.x[0]) / h, (py + 2.5 * np.sin(th) - g.y[0]) / h], order=1)
            a2 = ndimage.map_coordinates(q2, [(px + 2.5 * np.cos(th) - g.x[0]) / h, (py + 2.5 * np.sin(th) - g.y[0]) / h], order=1)
            q = round(np.angle(np.exp(1j * np.diff(np.r_[np.arctan2(a2, a1), np.arctan2(a2[0], a1[0])]))).sum() / (4 * np.pi) * 2) / 2
            d, j = g.tree.query([px, py]); inside = (px / g.A) ** 2 + (py / g.B) ** 2 < 1
            if q * g.ksign > 0 and d < 4 and inside:
                m = np.hypot(g.X - px, g.Y - py) < 1.6; tail = np.arctan2(dvy[m].mean(), dvx[m].mean()); taus.append(np.angle(np.exp(1j * (tail - g.nang[j])))); o["nwall"] += 1
            elif q > 0: o["nbp"] += 1
            elif q < 0: o["nbm"] += 1
            if q * g.ksign > 0 and d < 4 and inside and q < 0: taus.pop()   # no tilt for -1/2 walls
        if taus: z = np.mean(np.exp(1j * np.array(taus))); o["coh"] = float(abs(z)); o["tilt"] = float(np.degrees(np.angle(z)))
    m = g.inside & (g.depth > 3); r = np.maximum(np.hypot(g.X, g.Y), 0.5); vt = (g.X * vy - g.Y * vx) / r; vr = np.sqrt(np.mean((vx ** 2 + vy ** 2)[m]))
    o["circ"] = float(vt[m].mean() / max(vr, 1e-12)); o["vrms"] = float(vr)
    def rot(sel):
        cx, cy = g.X[sel].mean(), g.Y[sel].mean(); ddx, ddy = g.X[sel] - cx, g.Y[sel] - cy; return np.mean((ddx * vy[sel] - ddy * vx[sel]) / np.maximum(np.hypot(ddx, ddy), 0.5))
    o["two"] = int(rot(m & (g.X < 0)) * rot(m & (g.X > 0)) < 0); return o
def run(eps, zeta, lam=2.0, k=6, T0=600.0, T1=8000.0, dt=0.02, every=100.0, noise=0.005, seed=1, log=None, nowall=0.0, bandonly=0.0, actmode='all', swirl=0.0):
    g = Geo(eps, k); g.ksign = 1 if k > 0 else -1; rng = np.random.default_rng(seed)
    if nowall > 0: g.zmask = np.clip((g.depth - nowall) / 2.0, 0, 1)   # zeta = 0 within nowall cells of the wall, full beyond nowall + 2
    if bandonly > 0: g.zmask = 1 - np.clip((g.depth - bandonly) / 2.0, 0, 1)   # activity only within bandonly cells of the wall (ramp over 2)
    g.actmode = actmode; g.swirl = swirl
    if swirl:
        r = np.maximum(np.hypot(g.X, g.Y), 1e-9); m = (r < 0.6 * g.B) & g.inside; g.swirl_fx = -swirl * m * g.Y / r; g.swirl_fy = swirl * m * g.X / r
    q1 = np.where(g.inside | g.ghost, 0.5 * S0 * np.cos(2 * g.thw), 0); q2 = np.where(g.inside | g.ghost, 0.5 * S0 * np.sin(2 * g.thw), 0)   # nearest-wall director everywhere
    zero = (0, 0, 0, 0, 0)
    for _ in range(int(T0 / dt)): q1, q2 = step(g, q1, q2, dt, 0, zero, 0, rng)
    q1 += 0.01 * rng.standard_normal(q1.shape) * g.inside; q2 += 0.01 * rng.standard_normal(q2.shape) * g.inside
    lam_c = lam * S0; obs = [dict(t=0.0, **observe(g, q1, q2, np.zeros_like(q1), np.zeros_like(q1)))]; t = 0.0; t0 = time.time()
    while t < T1 - 1e-9:
        for _ in range(int(every / dt)): fl = flow(g, q1, q2, zeta, lam_c); q1, q2 = step(g, q1, q2, dt, lam_c, fl, noise, rng)
        t += every; o = observe(g, q1, q2, fl[0], fl[1]); o["t"] = t; obs.append(o)
        if log: print(f"eps {eps} zeta {zeta} t {t:.0f} ({time.time() - t0:.0f} s): wall {o['nwall']} bulk +{o['nbp']}/-{o['nbm']} coh {o.get('coh', np.nan):.2f} tilt {o.get('tilt', np.nan):+.0f} circ {o['circ']:+.2f} two {o['two']}", file=log, flush=True)
    return obs, (q1, q2)
if __name__ == "__main__":
    import sys
    ap = argparse.ArgumentParser(); ap.add_argument("--eps", type=float, default=1.0); ap.add_argument("--zeta", type=float, default=0.3); ap.add_argument("--lam", type=float, default=2.0)
    ap.add_argument("--k", type=int, default=6); ap.add_argument("--T0", type=float, default=600); ap.add_argument("--T1", type=float, default=8000); ap.add_argument("--noise", type=float, default=0.005)
    ap.add_argument("--seed", type=int, default=1); ap.add_argument("--nowall", type=float, default=0.0); ap.add_argument("--req", type=float, default=31.830989); ap.add_argument("--bandonly", type=float, default=0.0); ap.add_argument("--actmode", default="all"); ap.add_argument("--swirl", type=float, default=0.0); ap.add_argument("--out", default=None); a = ap.parse_args(); REQ = a.req; globals()['REQ'] = a.req
    obs, (q1, q2) = run(a.eps, a.zeta, a.lam, a.k, a.T0, a.T1, noise=a.noise, seed=a.seed, log=sys.stdout, nowall=a.nowall, bandonly=a.bandonly, actmode=a.actmode, swirl=a.swirl)
    if a.out: json.dump(dict(args=vars(a), obs=obs), open(a.out, "w")); np.savez_compressed(a.out.replace(".json", "_state.npz"), q1=q1, q2=q2)
