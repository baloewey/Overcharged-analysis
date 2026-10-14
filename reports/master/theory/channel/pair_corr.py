#!/usr/bin/env python3
"""Pair correlations of bulk defects: does the overcharge make the excess +1/2 a correlated, spaced-out liquid? (user, 2026-10-14)
Geometries: periodic box (neutral active turbulence, the control), winding channels (periodic in x), cavities (circles and ellipses; bulk = more than 4 cells
from the wall). For each group: g_ab(r) for ++, --, +- bulk pairs (minimum image where periodic), normalised by pairs between frames 200 t.u. apart
(same one-body density profile, no correlation; removes both the confinement shape and inhomogeneous densities); mean nearest-like-neighbour distance of +1/2 vs a Poisson gas at the same density.
Frames: second half of each run, every 5th frame. Usage: pair_corr.py -> pair_corr.txt, fig_pair_corr.png (+pdf)"""
import os, sys, glob, json, numpy as np
from multiprocessing import Pool
from matplotlib.path import Path
from scipy.spatial import cKDTree
sys.path.insert(0, "/home/baloewe/Physics/Overcharged/runs/cavity"); sys.path.insert(0, "/home/baloewe/Physics/Overcharged/reports/master/theory/energetics")
from edge_criteria import read_with_angle  # noqa
import energetics as En  # noqa
H = os.path.dirname(os.path.abspath(__file__)); R_ = "/home/baloewe/Physics/Overcharged/runs"
RB = np.arange(0, 40.5, 1.0); BAND = 4.0; rng = np.random.default_rng(3)


def geometry(run):
    J = json.load(open(run + "/input.json")); BC = J.get("BC", []); Lx, Ly = J["domain"][:2]
    if not BC or J.get("domainWalls", 0) == 1 and len(BC) == 0:
        return dict(kind="box", Lx=Lx, Ly=Ly, Q=np.zeros(2), inside=lambda p: np.ones(len(p), bool), px=Lx, py=Ly, area=Lx * Ly)
    if len(BC) >= 4 and "customAnch" in BC[2] and BC[2].get("eMode", 0) == 4 or (len(BC) >= 4 and abs(BC[2]["Q"][1] - BC[3]["Q"][1]) > 0 and BC[2].get("aInv", [0, 1, 0])[1] == 1):
        y0, y1 = BC[2]["Q"][1], BC[3]["Q"][1]
        return dict(kind="channel", Q=np.zeros(2), inside=lambda p: (p[:, 1] - y0 >= BAND) & (y1 - p[:, 1] >= BAND), px=Lx, py=None,
                    box=(0, Lx, y0 + BAND, y1 - BAND), area=Lx * (y1 - y0 - 2 * BAND))
    bc = BC[0]; W, s, P = En.wall_polyline(bc); tree = cKDTree(W); path = Path(W)
    ins = lambda p: path.contains_points(p) & (tree.query(p)[0] >= BAND)
    return dict(kind="cavity", Q=np.array(bc["Q"][:2], float), inside=ins, px=None, py=None, box=(W[:, 0].min(), W[:, 0].max(), W[:, 1].min(), W[:, 1].max()))


def dists(A, B, g, same):
    if len(A) == 0 or len(B) == 0: return np.zeros(0)
    dx = A[:, None, 0] - B[None, :, 0]; dy = A[:, None, 1] - B[None, :, 1]
    if g["px"]: dx = (dx + g["px"] / 2) % g["px"] - g["px"] / 2
    if g["py"]: dy = (dy + g["py"] / 2) % g["py"] - g["py"] / 2
    d = np.hypot(dx, dy)
    return d[np.triu_indices(len(A), 1)] if same else d.ravel()


def uniform_pdf(g, n=4000):
    if g["kind"] == "box": pts = rng.uniform([0, 0], [g["Lx"], g["Ly"]], (n, 2))
    else:
        x0, x1, y0, y1 = g["box"]; pts = np.zeros((0, 2))
        while len(pts) < n:
            c = rng.uniform([x0, y0], [x1, y1], (4 * n, 2)); pts = np.vstack([pts, c[g["inside"](c)]])
        pts = pts[:n]
    dd = dists(pts, pts, g, True); h = np.histogram(dd, RB)[0].astype(float); return h / len(dd), pts   # fraction of ALL pairs per bin


def one(run):
    g = geometry(run); pu, upts = uniform_pdf(g)
    t, fr = read_with_angle(run, g["Q"], 0.5); H3 = np.zeros((3, len(RB) - 1)); NP = np.zeros(3); nn, nnp, dens = [], [], []
    HX = np.zeros((3, len(RB) - 1)); NX = np.zeros(3); LAG = 100        # cross-frame reference: same one-body densities, no correlation (frames 200 t.u. apart)
    sel = lambda f: (lambda b: (f[b & (f[:, 2] > 0), :2], f[b & (f[:, 2] < 0), :2]))(g["inside"](f[:, :2]))
    for n_, f in enumerate(fr[::5]):
        if not len(f): continue
        P_, M_ = sel(f)
        for i, (A, B, same) in enumerate(((P_, P_, True), (M_, M_, True), (P_, M_, False))):
            d = dists(A, B, g, same); H3[i] += np.histogram(d, RB)[0]; NP[i] += len(d)
        j = 5 * n_ + LAG
        if j < len(fr) and len(fr[j]):
            P2, M2 = sel(fr[j])
            for i, (A, B) in enumerate(((P_, P2), (M_, M2), (P_, M2))):
                d = dists(A, B, g, False); HX[i] += np.histogram(d, RB)[0]; NX[i] += len(d)
        if len(P_) >= 2:
            d = dists(P_, P_, g, False).reshape(len(P_), len(P_)); np.fill_diagonal(d, np.inf); nn.extend(d.min(1))
            q = upts[rng.choice(len(upts), len(P_), replace=False)]; dq = dists(q, q, g, False).reshape(len(q), len(q)); np.fill_diagonal(dq, np.inf); nnp.extend(dq.min(1))
        dens.append((len(P_), len(M_)))
    gr = np.array([H3[i] / np.maximum(NP[i] * pu, 1e-12) for i in range(3)])
    area = len(upts) and (g.get("area") or None)
    return dict(run=os.path.relpath(run, R_), kind=g["kind"], H=H3, NP=NP, pu=pu, gr=gr, HX=HX, NX=NX, nn=np.mean(nn) if nn else np.nan, nnp=np.mean(nnp) if nnp else np.nan,
                npm=np.mean(dens, 0) if dens else np.array([np.nan, np.nan]))


GROUPS = {
    "periodic box, neutral (zeta 0.07)": sorted(glob.glob(f"{R_}/activeLength/L200/act_0.070*")),
    "circles R_eq 2 l_a (bulk +1)": sorted(glob.glob(f"{R_}/shapes_server/circle/circle/Req2.0/k*/s*")),
    "ellipses 2.0-2.5, R_eq 2 l_a": sorted(glob.glob(f"{R_}/shapes_server/ellipse/eps2.00/Req2.0/k*/s*") + glob.glob(f"{R_}/shapes_server/ellipse/eps2.50/Req2.0/k*/s*")),
    "channels W 40-60, q 0.031": sorted(glob.glob(f"{R_}/channel/W40_m6_L192_s*") + glob.glob(f"{R_}/channel/W60_m6_L192_s*")),
    "narrow channels W 16-24": sorted(glob.glob(f"{R_}/channel_narrow/W16_m6_s*") + glob.glob(f"{R_}/channel_narrow/W24_m6_s*")),
    "homeotropic channel W 40 (no winding)": sorted(glob.glob(f"{R_}/channel/W40_m0_L192_s*")),
}

if __name__ == "__main__":
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    runs = [r for v in GROUPS.values() for r in v]
    with Pool(int(os.environ.get("NPROC", 8))) as pool: res = dict(zip(runs, pool.map(one, runs, chunksize=1)))
    rc = 0.5 * (RB[1:] + RB[:-1]); L = ["PAIR CORRELATIONS of bulk defects (more than 4 cells from walls; second half, every 5th frame). g(r) = same-frame pairs / pairs between frames 200 t.u. apart (removes the inhomogeneous one-body density).",
                                        "group | runs | mean bulk +1/2, -1/2 per frame | g++ g-- g+- for r < 4 | for 4-8 | for 8-16 | position of the g++ maximum (cells) | <nn distance +1/2> / Poisson"]
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.4)); cols = ["0.4", "#1f5fa8", "#c2452d", "#7aa36b", "#e08a3c", "#7a3b8f"]
    for gi, (name, rr) in enumerate(GROUPS.items()):
        rr = [res[r] for r in rr if r in res]
        if not rr: continue
        H3 = sum(x["H"] for x in rr)
        NPu = sum(x["NP"][:, None] * x["HX"] / np.maximum(x["NX"][:, None], 1) for x in rr)   # expected same-frame counts from the cross-frame pair distribution
        G = H3 / np.maximum(NPu, 1e-12)
        def band(a, b): m = (rc >= a) & (rc < b); return H3[:, m].sum(1) / np.maximum(NPu[:, m].sum(1), 1e-12)
        b1, b2, b3 = band(0, 4), band(4, 8), band(8, 16); npm = np.nanmean([x["npm"] for x in rr], 0)
        pk = rc[np.argmax(np.where(rc < 30, G[0], 0))]; nnr = np.nanmean([x["nn"] / x["nnp"] for x in rr])
        L.append(f"{name:40s} | {len(rr):2d} | {npm[0]:5.1f} {npm[1]:5.1f} | {b1[0]:.2f} {b1[1]:.2f} {b1[2]:.2f} | {b2[0]:.2f} {b2[1]:.2f} {b2[2]:.2f} | {b3[0]:.2f} {b3[1]:.2f} {b3[2]:.2f} | {pk:4.1f} | {nnr:.2f}")
        for i in range(3):
            ok = NPu[i] > 20 * (NPu[i].sum() / len(NPu[i])) * 0 + 1
            ax[i].plot(rc, np.where(NPu[i] > 1, G[i], np.nan), "-", color=cols[gi % len(cols)], lw=1.8, label=name)
    for i, t in enumerate(("(a) $+1/2$ with $+1/2$", "(b) $-1/2$ with $-1/2$", "(c) $+1/2$ with $-1/2$")):
        ax[i].axhline(1, color="0.7", lw=0.8, ls=":"); ax[i].set_xlabel("separation r (cells)"); ax[i].set_ylabel("pair correlation g(r)"); ax[i].set_title(t, fontsize=10)
        ax[i].set_xlim(0, 30); [ax[i].spines[s_].set_visible(False) for s_ in ("top", "right")]
    ax[2].set_yscale("log")
    fig.legend(*ax[0].get_legend_handles_labels(), fontsize=8, frameon=False, ncol=3, loc="lower center"); fig.tight_layout(rect=(0, 0.13, 1, 1))
    fig.savefig(os.path.join(H, "fig_pair_corr.png"), dpi=200)
    print("\n".join(L)); open(os.path.join(H, "pair_corr.txt"), "w").write("\n".join(L) + "\n")
