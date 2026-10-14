#!/usr/bin/env python3
"""Elastic shielding by the wall layer and how it degrades in elongated cavities (user, 2026-10-14; cf. the 'bulk defects' potential' of the Japan poster).
VERSION 1, elastic potential felt by one bulk +1/2 with everything else fixed (one-constant, K = 32.5 kBT):
  measured frames (fine-output runs, area 2500, k 6): wall layer and all other defects fixed at their measured positions; one bulk +1/2 (the probe) is moved over
  a 2-cell grid of the cavity interior (> 2 cells from the wall), U(x) = E(probe at x) - min.
  Circle: analytic one-constant energy with wall images (strong anchoring theta = k phi; image of charge q at R^2/conj(z) with the SAME sign; validated against the
  grid solver), E = -pi K [ sum_{i != j} q_i q_j ln|z_i - z_j| + sum_{i,j} q_i q_j ln|R - z_i conj(z_j)/R| ].
  Ellipse: no closed form; numerical one-constant energy (ellipse_binding.Ellipse.energy, K1 = K3, cores relaxed, natural tilts).
  Cases: circle (full layer of 10); ellipse aspect 2 with an idealised full layer (10 evenly spaced in arclength, the other bulk +1/2 at -0.45 A); ellipse aspect 2,
  measured (depleted) layer. Maps averaged over several frames (circle, measured ellipse). Shielding measure: spread (s.d.) of U over the interior more than 6 cells
  from the wall and more than 4 cells from any fixed bulk defect.
VERSION 2, potential of mean force from the simulations: U_pmf(x) = -ln(rho_+(x)/<rho_+>) from the time-averaged +1/2 density (1-cell bins, Gaussian smoothing 1 cell),
  in units of an effective temperature; circle vs aspect-2 ellipse (same runs).
Usage: shielding.py -> shielding.txt, fig_shielding.png (+pdf)"""
import os, sys, glob, json, numpy as np
from multiprocessing import Pool
from scipy.spatial import cKDTree
from scipy.ndimage import gaussian_filter
R_ = "/home/baloewe/Physics/Overcharged/runs"; H = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault("GRID_H", "0.5")
sys.path.insert(0, R_ + "/cavity"); sys.path.insert(0, os.path.join(H, "..", "ellipse_binding")); sys.path.insert(0, os.path.join(H, "..", "energetics"))
from edge_criteria import read_with_angle  # noqa
import binding as B  # noqa
import energetics as En  # noqa
K = 32.5; GRID = 2.0; RUNS = {"circle": [f"{R_}/fineflow/eps1.0_s1", f"{R_}/fineflow2/eps1.0_s3"], "ellipse": [f"{R_}/fineflow/eps2.0_s1", f"{R_}/fineflow/eps2.0_s2"]}


def e_disk(zs, qs, R):
    zs = np.asarray(zs); qs = np.asarray(qs); d = np.abs(zs[:, None] - zs[None, :]); np.fill_diagonal(d, 1.0)
    img = np.abs(R - zs[:, None] * np.conj(zs[None, :]) / R)
    return -np.pi * K * (np.sum(np.outer(qs, qs) * np.log(d)) + np.sum(np.outer(qs, qs) * np.log(img)))


def frames(run, n=8):
    bc = json.load(open(run + "/input.json"))["BC"][0]; W, s, P = En.wall_polyline(bc); tree = cKDTree(W)
    t, fr = read_with_angle(run, np.array(bc["Q"][:2], float), 0.5); out = []
    for f in fr[::max(len(fr) // (4 * n), 1)]:
        if not len(f) or abs(0.5 * np.sign(f[:, 2]).sum() - 6) > 1e-6: continue
        d = tree.query(f[:, :2])[0]; bulkp = np.flatnonzero((f[:, 2] > 0) & (d > 6))
        if len(bulkp): out.append((f, d, bulkp[0]))
        if len(out) >= n: break
    return out


def grid_points(inside_fn, xr, yr):
    X, Y = np.meshgrid(np.arange(-xr, xr + 1e-9, GRID), np.arange(-yr, yr + 1e-9, GRID), indexing="ij"); Z = (X + 1j * Y).ravel()
    return X, Y, Z, inside_fn(Z)


def map_circle(args):
    f, d, ip = args; R = 2500 ** 0.5 / np.pi ** 0.5; zs = f[:, 0] + 1j * f[:, 1]; qs = 0.5 * np.sign(f[:, 2])
    X, Y, Z, ok = grid_points(lambda z: np.abs(z) < R - 2, R, R); U = np.full(Z.shape, np.nan)
    for n in np.flatnonzero(ok):
        z = zs.copy(); z[ip] = Z[n]; U[n] = e_disk(z, qs, R)
    others = np.delete(np.flatnonzero((f[:, 2] != 0) & (d > 4)), np.flatnonzero(np.flatnonzero((f[:, 2] != 0) & (d > 4)) == ip))
    return X, Y, U.reshape(X.shape), zs[others], zs[d <= 4]


def map_ellipse(args):
    zs, qs, ip = args; E = B.Ellipse(2.0); B.k = 6; A, Bb = E.A, E.B
    X, Y, Z, ok = grid_points(lambda z: (z.real / (A - 2)) ** 2 + (z.imag / (Bb - 2)) ** 2 < 1, A, Bb); U = np.full(Z.shape, np.nan)
    nw = int(np.sum(np.array([np.min(np.abs(E.W - zz)) for zz in zs]) <= 4))
    order = sorted(range(len(zs)), key=lambda j: np.min(np.abs(E.W - zs[j])))   # wall defects first (energy() treats the first entries as wall cores)
    zz = [zs[j] for j in order]; qq = [qs[j] for j in order]; ipp = order.index(ip)
    for n in np.flatnonzero(ok):
        z = list(zz); z[ipp] = Z[n]
        if min(abs(Z[n] - z[j]) for j in range(len(z)) if j != ipp) < 2.5: continue
        U[n] = E.energy(z, qq, [None] * nw, K1=K, K3=K)[0]
    return X, Y, U.reshape(X.shape)


def pmf(runs, shape):
    H2 = None
    for run in runs:
        bc = json.load(open(run + "/input.json"))["BC"][0]; t, fr = read_with_angle(run, np.array(bc["Q"][:2], float), 0.5)
        A = bc["R"] * (bc.get("aInv", [1])[0] or 1); Bb = bc["R"]; ex = np.arange(-A - 1, A + 1.01, 1.0); ey = np.arange(-Bb - 1, Bb + 1.01, 1.0)
        h = sum(np.histogram2d(f[f[:, 2] > 0, 0], f[f[:, 2] > 0, 1], [ex, ey])[0] for f in fr if len(f))
        H2 = h if H2 is None else H2 + h
    Xc, Yc = np.meshgrid(0.5 * (ex[1:] + ex[:-1]), 0.5 * (ey[1:] + ey[:-1]), indexing="ij"); ins = (Xc / A) ** 2 + (Yc / Bb) ** 2 < 1
    rho = gaussian_filter(np.where(ins, H2, 0.0), 1.0) / np.maximum(gaussian_filter(ins.astype(float), 1.0), 1e-9); rho[~ins] = np.nan
    return Xc, Yc, -np.log(rho / np.nanmean(rho)), A, Bb


if __name__ == "__main__":
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    L = ["ELASTIC SHIELDING: potential felt by one bulk +1/2 with everything else fixed (one-constant K = 32.5 kBT); PMF from the measured +1/2 density."]
    cf = [x for r in RUNS["circle"] for x in frames(r)]
    ef = [x for r in RUNS["ellipse"] for x in frames(r)]
    E0 = B.Ellipse(2.0); A, Bb = E0.A, E0.B
    ideal = ([E0.site(i / 10, 1.6)[0] for i in range(10)] + [0.45 * A + 0j, -0.45 * A + 0j], [0.5] * 12, 10)
    meas = [(list(f[:, 0] + 1j * f[:, 1]), list(0.5 * np.sign(f[:, 2])), ip) for f, d, ip in ef]
    with Pool(int(os.environ.get("NPROC", 6))) as pool:
        C = pool.map(map_circle, cf); EI = pool.map(map_ellipse, [ideal]); EM = pool.map(map_ellipse, meas)
    def spread(X, Y, U, inner, avoid):
        m = inner(X + 1j * Y) & np.isfinite(U)
        for z in avoid: m &= np.abs(X + 1j * Y - z) > 4
        return np.nanstd(U[m])
    Rc = (2500 / np.pi) ** 0.5
    sc = [spread(X, Y, U, lambda z: np.abs(z) < Rc - 6, oth) for X, Y, U, oth, _ in C]
    si = spread(*EI[0], lambda z: (z.real / (A - 6)) ** 2 + (z.imag / (Bb - 6)) ** 2 < 1, [-0.45 * A])
    sm = [spread(*m_, lambda z: (z.real / (A - 6)) ** 2 + (z.imag / (Bb - 6)) ** 2 < 1, [zz for zz, q in zip(mm[0], mm[1]) if abs(zz) < 0.9 * A and zz != mm[0][mm[2]]])
          for m_, mm in zip(EM, meas)]
    nwall = [int(np.sum(d <= 4)) for f, d, ip in ef]
    L.append(f"interior spread of U (s.d., kBT; > 6 cells from the wall, > 4 from fixed bulk defects): circle (measured, {len(C)} frames) median {np.median(sc):.1f}; "
             f"ellipse aspect 2, idealised full layer {si:.1f}; ellipse aspect 2, measured layer ({len(EM)} frames, wall defects {min(nwall)}-{max(nwall)}) median {np.median(sm):.1f}")
    # figure: (a) circle analytic, frame-averaged (relative), (b) ellipse ideal, (c) ellipse measured (one representative frame), (d) PMF circle, (e) PMF ellipse
    fig = plt.figure(figsize=(15, 8.4)); gs = fig.add_gridspec(2, 3)
    def show(ax, X, Y, U, title, pts=(), wall=None, vmax=60):
        U = U - np.nanmin(U); im = ax.pcolormesh(X, Y, np.minimum(U, vmax), cmap="magma_r", shading="nearest", vmin=0, vmax=vmax); ax.set_aspect("equal"); ax.axis("off")
        if wall is not None: ax.plot(wall.real, wall.imag, "k-", lw=1)
        for z, c in pts: ax.plot(z.real, z.imag, "o", ms=4, color=c, mec="w", mew=0.5)
        ax.set_title(title, fontsize=10); return im
    th = np.linspace(0, 2 * np.pi, 400); Wc = Rc * np.exp(1j * th); We = A * np.cos(th) + 1j * Bb * np.sin(th)
    Uc = np.nanmean([u - np.nanmin(u) for _, _, u, _, _ in C], 0)
    im = show(fig.add_subplot(gs[0, 0]), C[0][0], C[0][1], Uc, "(a) Circle, full layer (analytic, frame average)", [(z, "#1f5fa8") for z in C[0][4]], Wc)
    show(fig.add_subplot(gs[0, 1]), *EI[0], "(b) Ellipse 2, idealised full layer", [(z, "#1f5fa8") for z in ideal[0][:10]] + [(-0.45 * A, "#c2452d")], We)
    j = int(np.argsort(sm)[len(sm) // 2]); f, d, ip = ef[j]; zs = f[:, 0] + 1j * f[:, 1]
    show(fig.add_subplot(gs[0, 2]), *EM[j], f"(c) Ellipse 2, measured layer ({int(np.sum(d <= 4))} wall defects)",
         [(z, "#1f5fa8") for z in zs[(d <= 4) & (f[:, 2] > 0)]] + [(z, "#c2452d") for z in np.delete(zs, ip)[(np.delete(d, ip) > 4) & (np.delete(f[:, 2], ip) > 0)]]
         + [(z, "#e08a3c") for z in zs[f[:, 2] < 0]], We)
    cb = fig.colorbar(im, ax=fig.axes[:3], fraction=0.02, pad=0.01); cb.set_label("elastic potential of one bulk $+1/2$ ($k_BT$)")
    for i, (key, ttl) in enumerate((("circle", "(d) Circle: potential of mean force of $+1/2$"), ("ellipse", "(e) Ellipse 2: potential of mean force of $+1/2$"))):
        Xp, Yp, Up, A_, B_ = pmf(RUNS[key], key); ax = fig.add_subplot(gs[1, i]); im2 = ax.pcolormesh(Xp, Yp, np.clip(Up, -3, 3), cmap="RdBu_r", shading="nearest", vmin=-3, vmax=3)
        ax.plot(A_ * np.cos(th), B_ * np.sin(th), "k-", lw=1); ax.set_aspect("equal"); ax.axis("off"); ax.set_title(ttl, fontsize=10)
        inner = (Xp / (A_ - 6)) ** 2 + (Yp / (B_ - 6)) ** 2 < 1; L.append(f"PMF {key}: interior s.d. {np.nanstd(Up[inner]):.2f} kBT_eff; wall band minimum {np.nanmin(Up):.2f}")
    cb2 = fig.colorbar(im2, ax=fig.axes[-2:], fraction=0.03, pad=0.01); cb2.set_label(r"$-\ln(\rho_+/\langle\rho_+\rangle)$")
    ax = fig.add_subplot(gs[1, 2]); ax.axis("off")
    ax.text(0, 0.95, "Blue dots: wall $+1/2$ (fixed)\nRed dots: other bulk $+1/2$ (fixed)\nOrange: bulk $-1/2$ (fixed)\n\nTop row: one-constant elastic energy of\none bulk $+1/2$ moved over the cavity,\nall other defects fixed.\nBottom row: $-\\ln$ of the time-averaged\n$+1/2$ density from the simulations.", va="top", fontsize=9)
    fig.savefig(os.path.join(H, "fig_shielding.png"), dpi=200, bbox_inches="tight")
    print("\n".join(L)); open(os.path.join(H, "shielding.txt"), "w").write("\n".join(L) + "\n")
