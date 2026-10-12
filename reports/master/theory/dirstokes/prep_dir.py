#!/usr/bin/env python3
"""Full active force from the measured director field (no defect idealisation), per flow interval.
Q = S (n n - I/2) on the cell grid (directorfield.dat: t, i, j, k, nx, ny, nz, S; snapshots), active force density f = -div Q
(extensile; the overall factor zeta_eff and the sign are left to the fit, the physical sign is +). The flow output is averaged over
(t - dT, t], so f is averaged over the two snapshots at t - dT and t. Derivatives: centred differences between cells inside the cavity,
one-sided where a neighbour is outside, zero if both are. Split: wall band (cell centre within BAND = 4 cells of the wall) and bulk.
Writes <tag>.npz: x, y (cell centres inside, relative to the centre), fxw, fyw, fxb, fyb (intervals x cells), vx, vy (measured
interval-averaged flow at the same cells), t, A, B; and <tag>_grid.npz with the forces on the full grid for the Stokes solver."""
import gzip, json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); RUNS = "/home/baloewe/Physics/Overcharged/runs"
sys.path.insert(0, os.path.join(RUNS, "cavity")); from wall_fourier import wall_polyline, project  # noqa
BAND, NF = 4.0, 40
RUNLIST = [("k02", "cavity/Rscan_R31.83_seeds/s02/k02"), ("k06", "cavity/Rscan_R31.83_seeds/s02/k06"), ("k10", "cavity/Rscan_R31.83_seeds/s02/k10")] + \
          [(f"e{e.replace('.', '')}_{s}", f"cassini/W3_ellipse/eps{e}_k06/{s}") for e in ("1.2", "1.5", "2.0") for s in ("s01", "s02")]
def read_frames(fn, ncol):
    rows = np.loadtxt((l for l in gzip.open(fn, "rt") if len(l.split()) == ncol and l.split()[0][0].isdigit()))
    ts = np.unique(rows[:, 0]); return ts, rows
for tag, rel in RUNLIST:
    run = os.path.join(RUNS, rel); bc = json.load(open(run + "/input.json"))["BC"][0]; Q = bc["Q"]
    W, s, P = wall_polyline(bc); A, B = np.abs(W[:, 0]).max(), np.abs(W[:, 1]).max()
    td, Rd = read_frames(run + "/directorfield.dat.gz", 8); tf, Rf = read_frames(run + "/flowfield.dat.gz", 7)
    nx_, ny_ = int(Rd[:, 1].max()) + 1, int(Rd[:, 2].max()) + 1
    X = np.arange(nx_) + 0.5 - Q[0]; Y = np.arange(ny_) + 0.5 - Q[1]; XX, YY = np.meshgrid(X, Y, indexing="ij")
    inside = (XX / A) ** 2 + (YY / B) ** 2 < 1
    dw, _ = project(W, s, np.column_stack([XX.ravel(), YY.ravel()])); dw = dw.reshape(XX.shape)
    band = inside & (dw < BAND)
    def Qgrid(t):
        r = Rd[Rd[:, 0] == t]; G = np.zeros((nx_, ny_, 3))
        i, j = r[:, 1].astype(int), r[:, 2].astype(int); nxv, nyv, S = r[:, 4], r[:, 5], r[:, 7]
        G[i, j, 0] = S * (nxv * nxv - 0.5); G[i, j, 1] = S * nxv * nyv; G[i, j, 2] = -G[i, j, 0]
        return G
    def deriv(F, axis):
        out = np.zeros_like(F); m = inside
        fwd = np.roll(m, -1, axis); bwd = np.roll(m, 1, axis); Ff = np.roll(F, -1, axis); Fb = np.roll(F, 1, axis)
        both = m & fwd & bwd; out[both] = (Ff[both] - Fb[both]) / 2
        of = m & fwd & ~bwd; out[of] = Ff[of] - F[of]
        ob = m & ~fwd & bwd; out[ob] = F[ob] - Fb[ob]
        return out
    def force(G):
        fx = -(deriv(G[..., 0], 0) + deriv(G[..., 1], 1)); fy = -(deriv(G[..., 1], 0) + deriv(G[..., 2], 1))
        return fx * inside, fy * inside
    dT = tf[1] - tf[0]
    pick = [t for t in np.linspace(tf[len(tf) // 5], tf[-1], NF) ]
    pick = [tf[np.argmin(np.abs(tf - p))] for p in pick]
    vcell = inside & (dw > 1.0)
    out = dict(fxw=[], fyw=[], fxb=[], fyb=[], vx=[], vy=[], t=[])
    gw, gb = [], []
    for t in pick:
        t0 = t - dT
        if t0 not in set(td) or t not in set(td): t0 = td[np.argmin(np.abs(td - t0))]; t1 = td[np.argmin(np.abs(td - t))]
        else: t1 = t
        G = 0.5 * (Qgrid(t0) + Qgrid(t1)); fx, fy = force(G)
        r = Rf[Rf[:, 0] == t]; V = np.zeros((nx_, ny_, 2)); V[r[:, 1].astype(int), r[:, 2].astype(int)] = r[:, 4:6]
        out["fxw"].append((fx * band)[vcell]); out["fyw"].append((fy * band)[vcell])
        out["fxb"].append((fx * (inside & ~band))[vcell]); out["fyb"].append((fy * (inside & ~band))[vcell])
        out["vx"].append(V[..., 0][vcell]); out["vy"].append(V[..., 1][vcell]); out["t"].append(t)
        gw.append(np.stack([fx * band, fy * band], -1)); gb.append(np.stack([fx * (inside & ~band), fy * (inside & ~band)], -1))
    np.savez(os.path.join(HERE, f"{tag}.npz"), x=XX[vcell], y=YY[vcell], A=A, B=B, **{k: np.array(v) for k, v in out.items()})
    np.savez(os.path.join(HERE, f"{tag}_grid.npz"), X=X, Y=Y, wall=np.array(gw), bulk=np.array(gb), A=A, B=B)
    fw = np.hypot(out["fxw"], out["fyw"]).mean(); fb = np.hypot(out["fxb"], out["fyb"]).mean()
    print(f"{tag}: grid {nx_}x{ny_}, A {A:.2f} B {B:.2f}, dT {dT:g}, {len(pick)} intervals, |f| mean wall band {fw:.4f} bulk {fb:.4f}, |v| mean {np.hypot(out['vx'], out['vy']).mean():.3f}")
