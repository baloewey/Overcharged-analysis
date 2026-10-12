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
         "ESCAPERS: lag(t.u.) depth cos_tau  core  wall-side  bulk-side  net-lobe  wallL  bulkL |  residents (same lags): core  wall-side  bulk-side  net-lobe  wallL  bulkL"]
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
    fig, ax = plt.subplots(1, 3, figsize=(15, 4)); tl = 2 * np.array(LAGS)
    for Y, nm, col in ((wc, "core", "k"), (ww, "wall side (2.5 cells out)", "#c2452d"), (wb, "bulk side (2.5 cells in)", "#1f5fa8")):
        ax[0].plot(tl, [m_(Y[(kind == "esc") & (lag == l)]) for l in LAGS], "-o", ms=3, color=col, label=f"escaper: {nm}")
        cl = sorted(set(lag[kind == "ctl"])); ax[0].plot(2 * np.array(cl), [m_(Y[(kind == "ctl") & (lag == l)]) for l in cl], "--", color=col, lw=1, label=f"resident: {nm}")
    ax[0].set_ylabel("vorticity, co-rotating > 0 (1/t.u.)"); ax[0].set_xlabel("time relative to departure (t.u.)"); ax[0].set_title("(a) escaper vs resident", fontsize=10)
    ax2 = ax[0].twinx(); ax2.plot(tl, [m_(d[(kind == "esc") & (lag == l)]) for l in LAGS], ":", color="0.5"); ax2.set_ylabel("depth (dotted)", color="0.5")
    for Y, nm, col in ((wallL, "wall-side lobe", "#c2452d"), (bulkL, "bulk-side lobe", "#1f5fa8")):
        ax[1].plot(tl, [m_(Y[(kind == "esc") & (lag == l)]) for l in LAGS], "-o", ms=3, color=col, label=f"escaper: {nm}")
        cl = sorted(set(lag[kind == "ctl"])); ax[1].plot(2 * np.array(cl), [m_(Y[(kind == "ctl") & (lag == l)]) for l in cl], "--", color=col, lw=1, label=f"resident: {nm}")
    ax[1].axhline(1, color="0.8", lw=0.6, ls="--")
    ax[1].set_xlabel("time relative to departure (t.u.)"); ax[1].set_ylabel("lobe strength / free +1/2 value"); ax[1].set_title("(b) the two vortex lobes of the +1/2", fontsize=10)
    P_ = np.array(prof)
    for j_, nm, col in ((1, "core vorticity (1/t.u.)", "k"), (6, "wall-side lobe / free", "#c2452d"), (7, "bulk-side lobe / free", "#1f5fa8")): ax[2].plot(P_[:, 0], P_[:, j_] * (1 if j_ > 1 else 10), "-o", ms=3, color=col, label=nm if j_ > 1 else "core vorticity ×10 (1/t.u.)")
    ax[2].set_xlabel("depth (cells from the wall)"); ax[2].set_ylabel("value"); ax[2].set_title("(c) all +1/2 by depth", fontsize=10)
    for a_ in ax: a_.axhline(0, color="0.8", lw=0.6); a_.legend(fontsize=7, frameon=False)
    for a_ in ax[:2]: a_.axvline(0, color="0.6", lw=0.8, ls=":")
    fig.tight_layout(); fig.savefig(os.path.join(H, f"fig_vortdef_{tag}.png"), dpi=200)
