#!/usr/bin/env python3
"""Initial director with seeded defects for the MPCD solver option oDist = 10 (DIRFILE, branch seeddir) (user, 2026-10-12).
Field (as binding.py): theta = sum_i q_i arg(z - z_i) + psi; psi relaxed (anisotropic Frank, K1 42 / K3 23.1) with the eMode-2 anchoring k 2 pi s/P on the wall and constant
core values (cores relaxed). Then sampled at the MPCD cell centres and written as "ix iy theta" for dirInitFile.
Defects: N_wall +1/2 evenly in arclength at depth 1.6 (rotated by an offset in units of the spacing), bulk defects either given explicitly (x,y,q in cavity-centred coordinates)
or n_plus bulk +1/2 / n_minus bulk -1/2 spread evenly on the inner ellipse rho = 0.5 / 0.25. The total charge must equal k.
Usage: seed_director.py RUNDIR --wall 10 [--nplus 2 --nminus 0 | --bulk "x,y,q;x,y,q"] [--offset 0.5]
  reads RUNDIR/input.json (geometry, k), writes RUNDIR/dirinit.dat and sets oDist 10 + dirInitFile in RUNDIR/input.json; also RUNDIR/dirinit_seed.txt (defect list) and dirinit.png."""
import os, sys, json, argparse, numpy as np
TH = "/home/baloewe/Physics/Overcharged/reports/master/theory"; sys.path.insert(0, TH + "/ellipse_binding"); sys.path.insert(0, TH + "/energetics"); os.environ.setdefault("GRID_H", "0.5")
import binding as bd  # noqa
from energetics_general import EllipseAB  # noqa
import scipy.sparse as sp
from scipy.optimize import minimize
from scipy.interpolate import RegularGridInterpolator
def relaxed_theta(E, zs, qs, K1=42.0, K3=23.1):
    """binding.Ellipse.energy with all cores free; returns theta on E's grid (nan outside the cavity) and the energy."""
    Z = E.Z; e2 = np.ones(Z.shape, complex); gx = np.zeros(Z.shape); gy = np.zeros(Z.shape)
    for zi, q in zip(zs, qs):
        dz = Z - zi; dz = np.where(np.abs(dz) < 1e-9, 1e-9, dz); e2 *= (dz / np.abs(dz)) ** (2 * q); w = 1 / dz; gx += q * w.imag; gy += q * w.real
    Th2 = np.angle(e2).ravel(); Gx, Gy = gx.ravel(), gy.ravel()
    th_s_wall = np.zeros(len(E.W))
    for zi, q in zip(zs, qs): th_s_wall += q * np.unwrap(np.angle(E.W - zi))
    g = bd.k * 2 * np.pi * E.s / E.P - th_s_wall; g = np.unwrap(2 * g) / 2
    psi_b = np.zeros(Z.size); psi_b[E.band_idx] = g[E.band_near]
    core = -np.ones(Z.shape, int)
    for j, zi in enumerate(zs): core[(np.abs(Z - zi) < bd.a_core) & E.inside] = j
    free = E.inside & (core < 0); nF = free.sum()
    Pm = sp.csr_matrix((np.ones(nF), (np.flatnonzero(free.ravel()), np.arange(nF))), shape=(Z.size, nF))
    rows = np.flatnonzero((core >= 0).ravel()); Qm = sp.csr_matrix((np.ones(len(rows)), (rows, core.ravel()[rows])), shape=(Z.size, len(zs)))
    w = free.ravel().astype(float); Kb, De = 0.5 * (K1 + K3), 0.5 * (K1 - K3); h = bd.h
    def f(v):
        pf = psi_b + Pm @ v[:nF] + Qm @ v[nF:]
        gxx = Gx + E.Dx @ pf; gyy = Gy + E.Dy @ pf; C = np.cos(Th2 + 2 * pf); S = np.sin(Th2 + 2 * pf)
        dens = Kb * (gxx ** 2 + gyy ** 2) + De * (-C * (gxx ** 2 - gyy ** 2) - 2 * S * gxx * gyy); F = 0.5 * h * h * np.sum(w * dens)
        dgx = h * h * w * (Kb * gxx + De * (-C * gxx - S * gyy)); dgy = h * h * w * (Kb * gyy + De * (C * gyy - S * gxx)); dps = h * h * w * De * (S * (gxx ** 2 - gyy ** 2) - 2 * C * gxx * gyy)
        gf = E.Dx.T @ dgx + E.Dy.T @ dgy + dps; return F, np.concatenate([Pm.T @ gf, Qm.T @ gf])
    res = minimize(f, np.zeros(nF + len(zs)), jac=True, method="L-BFGS-B", options=dict(maxiter=30000, gtol=1e-8))
    pf = psi_b + Pm @ res.x[:nF] + Qm @ res.x[nF:]; theta = (Th2 / 2 + pf).reshape(Z.shape)
    return theta, res.fun
def main():
    ap = argparse.ArgumentParser(); ap.add_argument("rundir"); ap.add_argument("--wall", type=int, required=True); ap.add_argument("--nplus", type=int, default=0)
    ap.add_argument("--nminus", type=int, default=0); ap.add_argument("--bulk", default=""); ap.add_argument("--offset", type=float, default=0.5); ap.add_argument("--depth", type=float, default=1.6)
    a = ap.parse_args(); fn = os.path.join(a.rundir, "input.json"); J = json.load(open(fn)); bc = J["BC"][0]; k = bc["customAnch"]; bd.k = k
    ainv = bc.get("aInv", [1, 1, 1]); B = float(bc["R"]); A = B * float(ainv[0] or 1); E = EllipseAB(A, B); Q = bc["Q"]; nx, ny = J["domain"][0], J["domain"][1]
    zs, qs = [], []
    for j in range(a.wall):
        z, _ = E.site((j + a.offset) / a.wall, a.depth); zs.append(z); qs.append(0.5 * np.sign(k))
    if a.bulk:
        for item in a.bulk.split(";"):
            x, y, q = map(float, item.split(",")); zs.append(x + 1j * y); qs.append(q)
    for n_, q_, rho in ((a.nplus, 0.5, 0.5), (a.nminus, -0.5, 0.25)):
        for j in range(n_): phi = 2 * np.pi * (j + 0.25) / n_; zs.append(rho * A * np.cos(phi) + 1j * rho * B * np.sin(phi)); qs.append(q_)
    if abs(sum(qs) - k) > 1e-9: sys.exit(f"total charge {sum(qs)} != k = {k}: adjust --wall/--nplus/--nminus/--bulk")
    theta, Eel = relaxed_theta(E, np.array(zs), np.array(qs))
    # MPCD cell (ix, iy) centre in cavity coordinates: (ix + 0.5 - Q0, iy + 0.5 - Q1); interpolate cos/sin of 2 theta (nematic) on E's grid
    xg, yg = E.X[:, 0], E.Y[0, :]; c2 = np.where(np.isfinite(theta) & E.inside, np.cos(2 * theta), 0.0); s2 = np.where(np.isfinite(theta) & E.inside, np.sin(2 * theta), 0.0)
    ic = RegularGridInterpolator((xg, yg), c2, bounds_error=False, fill_value=0.0); is_ = RegularGridInterpolator((xg, yg), s2, bounds_error=False, fill_value=0.0)
    IX, IY = np.meshgrid(np.arange(nx), np.arange(ny), indexing="ij"); px = IX + 0.5 - Q[0]; py = IY + 0.5 - Q[1]
    inside = (px / A) ** 2 + (py / B) ** 2 < 1.0
    th_cell = 0.5 * np.arctan2(is_((px, py)), ic((px, py))); th_cell[~inside] = 0.5 * np.pi
    with open(os.path.join(a.rundir, "dirinit.dat"), "w") as fo:
        for ix_, iy_, t_ in zip(IX[inside], IY[inside], th_cell[inside]): fo.write(f"{ix_} {iy_} {t_:.6f}\n")
    with open(os.path.join(a.rundir, "dirinit_seed.txt"), "w") as fo:
        fo.write(f"# seeded defects (cavity-centred x y q), k {k}, wall {a.wall}, bulk +1/2 {sum(1 for z, q in zip(zs, qs) if q > 0) - a.wall}, -1/2 {sum(1 for q in qs if q < 0)}; Frank energy {Eel:.1f} kBT\n")
        for z, q in zip(zs, qs): fo.write(f"{z.real:.3f} {z.imag:.3f} {q:+.1f}\n")
    for s in J.get("species", []): s["oDist"] = 10
    J["dirInitFile"] = os.path.abspath(os.path.join(a.rundir, "dirinit.dat")); json.dump(J, open(fn, "w"), indent=1)
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(6, 6 * B / A + 0.5)); m = inside
        ax.quiver(px[m], py[m], np.cos(th_cell[m]), np.sin(th_cell[m]), headwidth=0, headlength=0, headaxislength=0, pivot="mid", scale=60, width=0.002)
        for z, q in zip(zs, qs): ax.plot(z.real, z.imag, "o" if q > 0 else "^", color="#c2452d" if q > 0 else "#1f5fa8", ms=7)
        tt = np.linspace(0, 2 * np.pi, 400); ax.plot(A * np.cos(tt), B * np.sin(tt), "k-"); ax.set_aspect("equal"); ax.set_title(f"seeded: wall {a.wall}, bulk charge {sum(qs) - 0.5 * a.wall * np.sign(k):+.1f} (k {k})")
        fig.savefig(os.path.join(a.rundir, "dirinit.png"), dpi=120, bbox_inches="tight")
    except Exception as ex: print("plot skipped:", ex)
    print(f"wrote {a.rundir}/dirinit.dat ({inside.sum()} cells); {len(zs)} defects, Frank energy {Eel:.1f} kBT")
if __name__ == "__main__": main()
