#!/usr/bin/env python3
"""Vorticity around +1/2 defects (user, 2026-10-11, from the movies: escapers sit near zero vorticity once detached; before, in signed vorticity, as if the wall suppressed the
lobe on its side). Fine-flow runs. Vorticity from the flow frame on the lattice, smoothed (Gaussian sigma 1 cell, fluid cells inside the cavity only), signed by the instantaneous
cavity circulation (sign of the mean wall-band tangential flow; > 0 = co-rotating with the cavity). Sampled (discs of radius 1.5): core; WALL side and BULK side (centre +- 2.5
cells along the inward normal); LEFT and RIGHT lobes (centre +- 2.5 cells perpendicular to the head; left = head rotated +90 deg). Lobe asymmetry: right - left (sign-free:
for a free +1/2 the two lobes are equal and opposite). Rows: depleting escapers (lags -30..+6 frames), resident controls (same frames), bulk +1/2 (depth > 8), and all +1/2 by depth.
Usage: vortdef.py TAG RUNDIR [...] -> vortdef_<TAG>.txt, vortdef_<TAG>.npy, fig_vortdef_<TAG>.png"""
import os, sys, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from multiprocessing import Pool
from scipy.spatial import cKDTree
from scipy.ndimage import gaussian_filter
from matplotlib.path import Path
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
import energetics as En  # noqa
from selfflow import load_flow  # noqa
LAGS = list(range(-30, 4))
def one(run):
    J, bc, k, t, fr, ev = En.events(run, max_events=50); Q = bc["Q"]; sg = 1 if k > 0 else -1
    W, s, P = En.wall_polyline(bc); tree = cKDTree(W); Wc = W[:, 0] + 1j * W[:, 1]; Tn = np.roll(Wc, -1) - np.roll(Wc, 1); Tn /= np.abs(Tn); Nn = -1j * Tn
    if np.mean((Nn * np.conj(Wc)).real) < 0: Nn = -Nn; Tn = -Tn
    rows, ft, idx = load_flow(run); ix = rows[:, 1].astype(int); iy = rows[:, 2].astype(int); nx, ny = ix.max() + 1, iy.max() + 1
    XX, YY = np.meshgrid(np.arange(nx) + 0.5 - Q[0], np.arange(ny) + 0.5 - Q[1], indexing="ij"); inside = Path(W).contains_points(np.c_[XX.ravel(), YY.ravel()]).reshape(XX.shape)
    dW, jW = tree.query(np.c_[XX.ravel(), YY.ravel()]); dW = dW.reshape(XX.shape); jW = jW.reshape(XX.shape); band = inside & (dW > 1) & (dW < 5)
    cache = {}
    def field(tt):
        t0 = ft[np.argmin(np.abs(ft - tt))]
        if t0 in cache: return cache[t0]
        ii = idx[t0]; U = np.zeros((nx, ny)); V = np.zeros((nx, ny)); U[ix[ii], iy[ii]] = rows[ii, 4]; V[ix[ii], iy[ii]] = rows[ii, 5]
        w8 = gaussian_filter(inside.astype(float), 1.0) + 1e-9; Us = gaussian_filter(np.where(inside, U, 0), 1.0) / w8; Vs = gaussian_filter(np.where(inside, V, 0), 1.0) / w8
        om = np.gradient(Vs, axis=0) - np.gradient(Us, axis=1); om[~inside] = np.nan
        gam = np.mean(U[band] * Tn[jW[band]].real + V[band] * Tn[jW[band]].imag); cache.clear(); cache[t0] = (om * np.sign(gam), gam); return cache[t0]
    def disc(om, c, r=1.5):
        m = (np.hypot(XX - c[0], YY - c[1]) < r) & np.isfinite(om); return om[m].mean() if m.any() else np.nan
    def sample(p, ang, m):
        om, gam = field(t[m]); d, j = tree.query(p); nin = -Nn[j]; nv = np.array([nin.real, nin.imag]); hd = -np.exp(1j * ang); lf = 1j * hd; lv = np.array([lf.real, lf.imag])
        sg_ = np.sign(gam) or 1.0; lw = float(np.dot(lv, nv) < 0)     # 1 if the LEFT lobe is on the wall side
        return (d, disc(om, p), disc(om, p - 2.5 * nv), disc(om, p + 2.5 * nv), disc(om, p + 2.5 * lv), disc(om, p - 2.5 * lv), abs(gam), sg_, lw, abs(np.dot(lv, nv)))
    def angle_at(p, m):
        f = fr[m]; return f[int(np.argmin(np.hypot(*(f[:, :2] - p).T))), 3]
    def ctv(p, ang):
        _, j = tree.query(p); return np.cos(np.angle(np.exp(1j * (ang - np.angle(Nn[j])))))
    out = []; rng = np.random.default_rng(9); A = np.abs(W[:, 0]).max()
    for q, e in enumerate(ev):
        n = e["n"]; ch = dict(e["chain"]); c = ch[n]
        for m in range(n + 1, n + 5):
            h = fr[m]; hp = h[np.sign(h[:, 2]) == sg, :2] if len(h) else np.zeros((0, 2))
            if not len(hp): break
            dd = np.hypot(*(hp - c).T); kk = np.argmin(dd)
            if dd[kk] > 2.5: break
            c = hp[kk]; ch[m] = c
        for lag in LAGS:
            m = n + lag
            if m not in ch: continue
            a = angle_at(ch[m], m); out.append(("esc", q, lag, ctv(ch[m], a)) + sample(ch[m], a, m))
            if lag % 5 == 0:
                f = fr[m]; dd = tree.query(f[:, :2])[0]
                res = [jj for jj in np.flatnonzero((np.sign(f[:, 2]) == sg) & (dd < 4)) if np.hypot(*(f[jj, :2] - ch[m])) > 3 and abs(abs(f[jj, 0]) - abs(ch[m][0])) / A < 0.15]
                if res: jj = rng.choice(res); out.append(("ctl", q, lag, ctv(f[jj, :2], f[jj, 3])) + sample(f[jj, :2], f[jj, 3], m))
    for m in range(len(fr) // 2, len(fr), 9):                          # all +1/2 (depth map, bulk reference)
        f = fr[m]
        if not len(f): continue
        for fd in f[f[:, 2] > 0]: out.append(("all", -1, 0, ctv(fd[:2], fd[3])) + sample(fd[:2], fd[3], m))
    return out
if __name__ == "__main__":
    tag = sys.argv[1]; runs = [os.path.abspath(r) for r in sys.argv[2:]]
    if os.environ.get("REPLOT"):   # redraw from the saved samples
        a = np.load(os.path.join(H, f"vortdef_{tag}.npy"), allow_pickle=True)
    else:
      with Pool(min(12, len(runs))) as pool: res = pool.map(one, runs, chunksize=1)
      a = np.array(sum(res, []), dtype=object); np.save(os.path.join(H, f"vortdef_{tag}.npy"), a, allow_pickle=True)
      np.save(os.path.join(H, f"vortdef_{tag}_runs.npy"), np.array([(r, i) for r, o in zip(runs, res) for i in range(len(o))], dtype=object), allow_pickle=True)   # run of each row (for the region split)
    kind = a[:, 0].astype(str); lag = a[:, 2].astype(int); ct, d, wc, ww, wb, wl, wr, gam, sgm, lw, lat = (a[:, i].astype(float) for i in range(3, 14))
    wl, wr = wl * sgm, wr * sgm                                       # undo the circulation sign: raw lobe vorticity in the defect frame
    bulkm = (kind == "all") & (d > 8); O0 = np.nanmean(wl[bulkm] - wr[bulkm]) / 2         # free +1/2: left = +O0, right = -O0
    wallL = np.where(lw > 0, wl / O0, -wr / O0); bulkL = np.where(lw > 0, -wr / O0, wl / O0)  # lobe strengths in units of the free value (1 = free)
    wallL[lat < 0.5] = np.nan; bulkL[lat < 0.5] = np.nan                                 # only when the lobes straddle the wall normal (head roughly along the wall)
    asym = (wl + wr) / 2 / O0                                                             # net (non-dipolar) lobe vorticity in units of O0
    m_ = lambda x: np.nanmean(x); se = lambda x: np.nanstd(x) / np.sqrt(max(np.isfinite(x).sum(), 1))
    L = [f"vorticity signed by the cavity circulation (> 0 co-rotating); |wall-band circulation| {np.nanmean(gam):.3f}; units 1/t.u.; free +1/2 lobe strength O0 = {O0:+.4f} (bulk, depth > 8)",
         "lobes: wallL / bulkL = strength of the lobe on the wall / bulk side in units of the free value (only when |left.normal| > 0.5); net = (left+right)/2/O0",
         "ESCAPERS: lag(t.u.) depth cos_tau  core  wall-side  bulk-side  net-lobe  wallL  bulkL |  controls (same lags): core  wall-side  bulk-side  net-lobe  wallL  bulkL"]
    for lg in LAGS:
        e = (kind == "esc") & (lag == lg); c = (kind == "ctl") & (lag == lg)
        s = f"{2*lg:+4d}  {m_(d[e]):.2f}  {m_(ct[e]):+.2f}   {m_(wc[e]):+.4f}±{se(wc[e]):.4f}  {m_(ww[e]):+.4f}  {m_(wb[e]):+.4f}  {m_(asym[e]):+.2f}  {m_(wallL[e]):+.2f}  {m_(bulkL[e]):+.2f}"
        if c.any(): s += f"  |  {m_(wc[c]):+.4f}±{se(wc[c]):.4f}  {m_(ww[c]):+.4f}  {m_(wb[c]):+.4f}  {m_(asym[c]):+.2f}  {m_(wallL[c]):+.2f}  {m_(bulkL[c]):+.2f}"
        L.append(s)
    L.append("\nALL +1/2 by depth: depth  n  core  wall-side  bulk-side  net-lobe  |core|  wallL  bulkL")
    DB = [0, 1, 2, 3, 4, 5, 6, 8, 10, 99]; prof = []
    for i in range(len(DB) - 1):
        al = (kind == "all") & (d >= DB[i]) & (d < DB[i + 1])
        if al.sum() > 20:
            prof.append((0.5 * (DB[i] + min(DB[i + 1], 14)), m_(wc[al]), se(wc[al]), m_(ww[al]), m_(wb[al]), m_(asym[al]), m_(wallL[al]), m_(bulkL[al])))
            L.append(f"  {DB[i]}-{DB[i+1]}  {al.sum():5d}  {m_(wc[al]):+.4f}±{se(wc[al]):.4f}  {m_(ww[al]):+.4f}  {m_(wb[al]):+.4f}  {m_(asym[al]):+.2f}  {m_(np.abs(wc[al])):.4f}  {m_(wallL[al]):+.2f}  {m_(bulkL[al]):+.2f}")
    print("\n".join(L)); open(os.path.join(H, f"vortdef_{tag}.txt"), "w").write("\n".join(L) + "\n")
    # Figure (user, 2026-10-13: make the y quantities explicit). (a) sketch of the sampling geometry; (b) circulation-signed vorticity at three points along the wall normal;
    # (c) lobe strengths in the defect frame relative to a free +1/2; (d) the same by depth for all +1/2. Colour = where the vorticity is read, style = escaper / control.
    from matplotlib.lines import Line2D
    from matplotlib.patches import Circle, FancyArrowPatch, Rectangle
    fig = plt.figure(figsize=(19, 5.4)); gs = fig.add_gridspec(1, 4, width_ratios=[0.8, 1, 1, 1]); ax = [fig.add_subplot(gs[i]) for i in range(4)]; tl = 2 * np.array(LAGS)
    CW, CB, CC = "#c2452d", "#1f5fa8", "k"
    ax[0].remove(); sg0 = gs[0].subgridspec(2, 1, hspace=0.25); S = [fig.add_subplot(sg0[0]), fig.add_subplot(sg0[1])]
    p0 = np.array([0.0, 3.0]); th = np.deg2rad(30); hd = np.array([np.cos(th), np.sin(th)]); lv = np.array([-hd[1], hd[0]]); nv = np.array([0.0, 1.0])
    def seg(s_, a_, b_, **kw): s_.plot([a_[0], b_[0]], [a_[1], b_[1]], **kw)
    for j, s_ in enumerate(S):
        s_.set_aspect("equal"); s_.set_xlim(-6, 6); s_.set_ylim(-1.0, 7.0); s_.axis("off")
        s_.add_patch(Rectangle((-6, -1.0), 12, 1.0, color="0.85")); s_.plot([-6, 6], [0, 0], color="0.3", lw=1.2); s_.text(5.8, -0.6, "wall", ha="right", fontsize=7, color="0.3")
        seg(s_, p0 - 1.4 * hd, p0, color="0.2", lw=2); seg(s_, p0 - 1.1 * hd + 0.4 * lv, p0, color="0.2", lw=1); seg(s_, p0 - 1.1 * hd - 0.4 * lv, p0, color="0.2", lw=1)
        s_.add_patch(FancyArrowPatch(p0, p0 + 1.3 * hd, arrowstyle="-|>", mutation_scale=9, color="0.2"))
        if j == 0:
            pts = ((p0, CC, "core", (1.7, -0.9)), (p0 - 2.5 * nv, CW, "2.5 towards wall", (1.7, 0.0)), (p0 + 2.5 * nv, CB, "2.5 towards bulk", (1.7, 0.3))); mk = "s"; ttl = "(a) Along the wall normal, for (b)"
        else:
            pts = ((p0 - 2.5 * lv, CW, "wall-side lobe", (1.7, -0.3)), (p0 + 2.5 * lv, CB, "bulk-side lobe", (-1.0, 1.9))); mk = "o"; ttl = "Perpendicular to the head, for (c, d)"
            s_.text(*(p0 + 1.4 * hd + [0.1, 0.1]), "head", fontsize=6.5)
        for c_, col, lab, off in pts:
            s_.add_patch(Circle(c_, 1.5, fill=False, lw=1.1, color=col)); s_.plot(*c_, mk, ms=3.5, color=col); s_.text(c_[0] + off[0], c_[1] + off[1], lab, fontsize=6.5, color=col, va="center")
        s_.set_title(ttl, fontsize=8)
    for Y, col in ((wc, CC), (ww, CW), (wb, CB)):
        ax[1].plot(tl, [m_(Y[(kind == "esc") & (lag == l)]) for l in LAGS], "-o", ms=3, color=col)
        cl = sorted(set(lag[kind == "ctl"])); ax[1].plot(2 * np.array(cl), [m_(Y[(kind == "ctl") & (lag == l)]) for l in cl], "--", color=col, lw=1)
    ax[1].set_ylabel("vorticity $\\omega$ (1/t.u.)"); ax[1].set_xlabel("time relative to departure (t.u.)")
    ax[1].set_title("(b) Vorticity, $>0$ co-rotating with the cavity", fontsize=10)
    ax2 = ax[1].twinx(); ax2.plot(tl, [m_(d[(kind == "esc") & (lag == l)]) for l in LAGS], ":", color="0.5"); ax2.set_ylabel("escaper depth (cells)", color="0.5")
    for Y, col in ((wallL, CW), (bulkL, CB)):
        ax[2].plot(tl, [m_(Y[(kind == "esc") & (lag == l)]) for l in LAGS], "-o", ms=3, color=col)
        cl = sorted(set(lag[kind == "ctl"])); ax[2].plot(2 * np.array(cl), [m_(Y[(kind == "ctl") & (lag == l)]) for l in cl], "--", color=col, lw=1)
    ax[2].axhline(1, color="0.8", lw=0.6, ls="--"); ax[2].text(-40, 1.04, "free +1/2", fontsize=7, color="0.5")
    ax[2].set_xlabel("time relative to departure (t.u.)"); ax[2].set_ylabel("lobe vorticity / $\\Omega_0$"); ax[2].set_title(f"(c) Lobes relative to a free $+1/2$ ($\\Omega_0$ = {O0:.2f}/t.u.)", fontsize=10)
    P_ = np.array(prof)
    ax[3].plot(P_[:, 0], P_[:, 6], "-o", ms=3, color=CW); ax[3].plot(P_[:, 0], P_[:, 7], "-o", ms=3, color=CB); ax[3].axhline(1, color="0.8", lw=0.6, ls="--")
    ax[3].set_ylabel("lobe vorticity / $\\Omega_0$"); ax[3].set_xlabel("depth (cells from the wall)"); ax[3].set_title("(d) all +1/2 by depth", fontsize=10)
    a3 = ax[3].twinx(); a3.plot(P_[:, 0], P_[:, 1], "-s", ms=3, color=CC); a3.set_ylabel("core $\\omega$ (1/t.u.)"); a3.axhline(0, color="0.85", lw=0.6)
    sty = [Line2D([], [], color="0.3", ls="-", marker="o", ms=3, label="escaper"), Line2D([], [], color="0.3", ls="--", lw=1, label="control")]
    ax[1].legend(handles=[Line2D([], [], color=c, marker="s", lw=2, label=n) for c, n in ((CC, "at the core"), (CW, "2.5 cells towards the wall"), (CB, "2.5 cells towards the bulk"))] + sty
                 + [Line2D([], [], color="0.5", ls=":", label="escaper depth (right axis)")], fontsize=7, frameon=False, ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.17))
    ax[2].legend(handles=[Line2D([], [], color=c, marker="o", lw=2, label=n) for c, n in ((CW, "wall-side lobe"), (CB, "bulk-side lobe"))] + sty,
                 fontsize=7, frameon=False, ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.17))
    ax[3].legend(handles=[Line2D([], [], color=CW, marker="o", label="wall-side lobe / $\\Omega_0$"), Line2D([], [], color=CB, marker="o", label="bulk-side lobe / $\\Omega_0$"),
                          Line2D([], [], color=CC, marker="s", label="core $\\omega$ (right axis)")], fontsize=7, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=2)
    for a_ in ax[1:]: a_.axhline(0, color="0.8", lw=0.6)
    for a_ in ax[1:3]: a_.axvline(0, color="0.6", lw=0.8, ls=":")
    fig.tight_layout(); fig.savefig(os.path.join(H, f"fig_vortdef_{tag}.png"), dpi=200)
