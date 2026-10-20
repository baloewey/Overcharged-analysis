#!/usr/bin/env python3
"""Q-tensor bench with backflow and activity (extends shear_tilt.py). Same Q dynamics; the flow is now SOLVED each few steps from Stokes with no-slip:
  eta lap^2 psi = d_x f_y - d_y f_x,  v = (d_y psi, -d_x psi),  psi = 0 on the ghost ring (clamped: no-slip, staircase wall),
  f = div(sigma) + f_body,  sigma = -zeta Q (active, extensile for zeta > 0) - lam_c H + (Q H - H Q) (backflow),  f_body = optional azimuthal A within 0.6 R.
eta = 1 (rotational viscosity of the Q model gamma_1 = 2 S^2 / Gamma = 1.53). K = 2 L S^2 = 1.535 -> active length sqrt(K/zeta).
Validation (--validate): uniform azimuthal force f in r < R gives v_theta = f r (R - r) / (3 eta)."""
import argparse, json, os, sys
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spla
import shear_tilt as M
X, Y, h, R, inside, ghost = M.X, M.Y, M.h, M.R, M.inside, M.ghost
n = X.shape[0]; I1 = sp.identity(n, format="csr"); D2 = sp.diags([1, -2, 1], [-1, 0, 1], shape=(n, n)) / h ** 2
LAP = (sp.kron(D2, I1) + sp.kron(I1, D2)).tocsr(); idx = np.flatnonzero(inside.ravel())
BI = (LAP @ LAP).tocsr()[idx][:, idx].tocsc(); SOLVE = spla.factorized(BI)
def dx(f): return (np.roll(f, -1, 0) - np.roll(f, 1, 0)) / (2 * h)
def dy(f): return (np.roll(f, -1, 1) - np.roll(f, 1, 1)) / (2 * h)
def stokes(fx, fy, eta=1.0):
    rhs = (dx(fy) - dy(fx)).ravel()[idx] / eta; psi = np.zeros(X.size); psi[idx] = SOLVE(rhs); psi = psi.reshape(X.shape)
    return dy(psi), -dx(psi)
def forces(q1, q2, zeta, lam_c, A_body):
    trQ2 = 2 * (q1 ** 2 + q2 ** 2); h1 = M.Lc * M.lap(q1) + M.a * q1 - M.c * trQ2 * q1; h2 = M.Lc * M.lap(q2) + M.a * q2 - M.c * trQ2 * q2
    h1 = np.where(inside, h1, 0); h2 = np.where(inside, h2, 0); asym = 2 * (q1 * h2 - q2 * h1)            # (QH - HQ)_xy
    sxx = -zeta * q1 - lam_c * h1; sxy = -zeta * q2 - lam_c * h2 + asym; syx = -zeta * q2 - lam_c * h2 - asym; syy = zeta * q1 + lam_c * h1
    fx = dx(sxx) + dy(sxy); fy = dx(syx) + dy(syy)
    if A_body:
        m = np.hypot(X, Y) < SWIRL_FRAC * R; r = np.maximum(np.hypot(X, Y), 1e-9); fx = fx - A_body * m * Y / r; fy = fy + A_body * m * X / r
    return fx, fy
def flowfields(vx, vy):
    Exx, Exy, Oxy = 0.5 * (dx(vx) - dy(vy)), 0.5 * (dy(vx) + dx(vy)), 0.5 * (dy(vx) - dx(vy)); return vx, vy, Exx, Exy, Oxy
def circulation(vx, vy):
    r = np.maximum(np.hypot(X, Y), 1e-9); vt = (X * vy - Y * vx) / r; m = inside & (r < R - 3)
    return float(vt[m].mean() / max(np.sqrt(np.mean((vx ** 2 + vy ** 2)[m])), 1e-12)), float(np.sqrt(np.mean((vx ** 2 + vy ** 2)[m])))
SWIRL_FRAC = 0.6   # body force acts inside SWIRL_FRAC * R (default as before; > 1 = everywhere in the disk), set from active_variants via env SWIRL_FRAC
def run(lam, zeta, A_body=0.0, ring=12, T0=300.0, T1=6000.0, dt=0.04, every=250.0, nflow=5, eta=1.0, seed=1):
    rng = np.random.default_rng(seed)
    zs = list((R - 2.6) * np.exp(2j * np.pi * (np.arange(ring) + 0.5) / ring)); qs = [0.5] * ring
    if ring == 2 * M.k - 2: zs = list((R - 2.6) * np.exp(2j * np.pi * np.arange(ring) / ring)) + [0.45j * R, -0.45j * R]; qs = [0.5] * (ring + 2)
    th2, _, _ = M.fields(X, Y, zs, qs, R); q1 = np.where(inside, 0.5 * M.S0 * np.cos(th2), M.qb1); q2 = np.where(inside, 0.5 * M.S0 * np.sin(th2), M.qb2)
    for _ in range(int(T0 / dt)): q1, q2 = M.step(q1, q2, dt, 0, None)
    q1 = q1 + 0.01 * rng.standard_normal(q1.shape) * inside; q2 = q2 + 0.01 * rng.standard_normal(q2.shape) * inside   # break symmetry
    lam_c = lam * M.S0; obs = []; t = 0.0; fl = None; i = 0
    while t < T1 - 1e-9:
        for _ in range(int(every / dt)):
            if i % nflow == 0: vx, vy = stokes(*forces(q1, q2, zeta, lam_c, A_body), eta=eta); fl = flowfields(vx, vy)
            q1, q2 = M.step(q1, q2, dt, lam_c, fl); i += 1
        t += every; circ, vrms = circulation(vx, vy); o = M.observe(q1, q2, np.sign(circ) or 1); o.update(t=t, circ=circ, vrms=vrms); obs.append(o)
    return obs, (q1, q2, vx, vy)
if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--validate", action="store_true"); ap.add_argument("--lam", type=float, default=2.0)
    ap.add_argument("--zeta", type=float, default=0.01); ap.add_argument("--A", type=float, default=0.0); ap.add_argument("--ring", type=int, default=12)
    ap.add_argument("--T1", type=float, default=6000); ap.add_argument("--out", default=None)
    a = ap.parse_args()
    if a.validate:
        r = np.maximum(np.hypot(X, Y), 1e-9); vx, vy = stokes(-Y / r * inside, X / r * inside); vt = (X * vy - Y * vx) / r
        for rr in (5, 10, 15, 20, 25, 28, 30):
            m = inside & (abs(r - rr) < 0.5); print(rr, round(vt[m].mean(), 3), "analytic", round(rr * (R - rr) / 3, 3))
        sys.exit()
    obs, _ = run(a.lam, a.zeta, a.A, a.ring, T1=a.T1)
    for o in obs: print(o["t"], o["nwall"], o.get("nbulk_p"), o.get("nbulk_m"), round(o.get("coh", np.nan), 2), round(o.get("tilt", np.nan), 1), round(o["circ"], 2), f"{o['vrms']:.2e}", flush=True)
    if a.out: json.dump(dict(args=vars(a), obs=obs), open(a.out, "w"))
