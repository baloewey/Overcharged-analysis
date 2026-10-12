#!/usr/bin/env python3
"""Zoomed movie of one depleting escape (fineflow run: flow + defects every 2 t.u.). Background = vorticity; arrows = flow; +1/2 comets (head = -angle), -1/2 tripods;
escaper highlighted with its trail, its head arrow and the local inward normal; panel text: time to departure, cos(tau), depth, inward flow at core and ring.
Usage: escape_movie.py RUNDIR [event_index] -> frames/ + escape_<run>_<ev>.mp4"""
import os, sys, subprocess
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from scipy.ndimage import gaussian_filter
from matplotlib.path import Path
H = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, H)
import energetics as En  # noqa
from selfflow import load_flow  # noqa
run = sys.argv[1].rstrip("/"); J, bc, k, t, fr, ev = En.events(run, max_events=999); Q = bc["Q"]; sg = 1 if k > 0 else -1
W, s, P = En.wall_polyline(bc); tree = cKDTree(W); Wc = W[:, 0] + 1j * W[:, 1]; Tn = np.roll(Wc, -1) - np.roll(Wc, 1); Nn = -1j * Tn / np.abs(Tn)
if np.mean((Nn * np.conj(Wc)).real) < 0: Nn = -Nn
rows, ft, idx = load_flow(run); X = rows[:, 1] + 0.5 - Q[0]; Y = rows[:, 2] + 0.5 - Q[1]
def flow(tt):
    ii = idx[ft[np.argmin(np.abs(ft - tt))]]; return X[ii], Y[ii], rows[ii, 4], rows[ii, 5]
def follow(e, back=30, fwd=12):
    n = e["n"]; ch = dict(e["chain"])
    for rng_, start in ((range(n + 1, n + fwd + 1), n), (range(min(ch) - 1, n - back - 1, -1), min(ch))):
        c = ch[start]
        for m in rng_:
            h = fr[m]; hp = h[np.sign(h[:, 2]) == sg, :2] if len(h) else np.zeros((0, 2))
            if not len(hp): break
            d = np.hypot(*(hp - c).T); kk = np.argmin(d)
            if d[kk] > 2.5: break
            c = hp[kk]; ch[m] = c
    return ch
def cost(p, m):
    f = fr[m]; i = int(np.argmin(np.hypot(*(f[:, :2] - p).T))); _, j = tree.query(p); return np.cos(np.angle(np.exp(1j * (f[i, 3] - np.angle(Nn[j]))))), f[i, 3]
O0 = 0.118   # free bulk +1/2 lobe vorticity (vortdef_ff12, sigma-1 smoothing as here)
print(len(ev), "depleting escapes")
if len(sys.argv) > 2: ie = int(sys.argv[2])
else:   # pick the clearest: long resident chain, clean exit, largest untilting
    sc = []
    for q, e in enumerate(ev):
        ch = follow(e); n = e["n"]
        if not all(m in ch for m in range(n - 30, n + 13)): sc.append(-9); continue
        sc.append(cost(ch[n], n)[0] - cost(ch[n - 25], n - 25)[0] + 0.05 * (tree.query(ch[n + 12])[0] > 8))
    ie = int(np.argmax(sc)); print("scores", np.round(sc, 2))
e = ev[ie]; ch = follow(e); n = e["n"]; site = ch[n]; _, j0 = tree.query(site); nin0 = -Nn[j0]
cen = site + 5 * np.array([nin0.real, nin0.imag]); L = 13
os.makedirs("frames", exist_ok=True); [os.remove("frames/" + f) for f in os.listdir("frames")]
ms = sorted(m for m in ch if n - 30 <= m <= n + 12); hist = []
for q, m in enumerate(ms):
    x, y, vx, vy = flow(t[m]); win = (abs(x - cen[0]) < L + 1) & (abs(y - cen[1]) < L + 1)
    xs, ys = np.unique(x[win]), np.unique(y[win]); G = np.full((len(ys), len(xs), 2), np.nan)
    G[np.searchsorted(ys, y[win]), np.searchsorted(xs, x[win])] = np.c_[vx[win], vy[win]]
    XX, YY = np.meshgrid(xs, ys); pts = np.c_[XX.ravel(), YY.ravel()]; inside = (Path(W).contains_points(pts) | (tree.query(pts)[0] < 1.5)).reshape(XX.shape) & np.isfinite(G[..., 0])
    Gs = np.zeros_like(G); wt = gaussian_filter(inside.astype(float), 1.0) + 1e-9
    for c_ in (0, 1): Gs[..., c_] = gaussian_filter(np.where(inside, G[..., c_], 0), 1.0) / wt
    Gs[~inside] = 0; w = np.gradient(Gs[..., 1], axis=1) - np.gradient(Gs[..., 0], axis=0); w[~inside] = np.nan; Gs[~inside] = np.nan
    p = ch[m]; d, j = tree.query(p); nin = -Nn[j]; r = np.hypot(x - p[0], y - p[1])
    vc = np.mean(vx[r < 2] * nin.real + vy[r < 2] * nin.imag); rg = (r > 4) & (r < 8); vr = np.mean(vx[rg] * nin.real + vy[rg] * nin.imag)
    ct, ang = cost(p, m)
    hd = -np.exp(1j * ang); lf = 1j * hd; lv = np.array([lf.real, lf.imag]); nv = np.array([nin.real, nin.imag]); XXw, YYw = np.meshgrid(xs, ys)
    def disc(c_, r_=1.5):
        mm_ = (np.hypot(XXw - c_[0], YYw - c_[1]) < r_) & np.isfinite(w); return w[mm_].mean() if mm_.any() else np.nan
    cL, cR = p + 2.5 * lv, p - 2.5 * lv; wl, wr = disc(cL), disc(cR); left_wall = np.dot(lv, nv) < 0      # free +1/2: left lobe = +O0, right = -O0
    wallS, bulkS = (wl / O0, -wr / O0) if left_wall else (-wr / O0, wl / O0); cW, cB = (cL, cR) if left_wall else (cR, cL)
    hist.append((t[m] - t[n], ct, vc, vr, wallS, bulkS))
    fig = plt.figure(figsize=(10.5, 7.2), dpi=110); ax = fig.add_axes([0.03, 0.08, 0.58, 0.86]); ax.set_aspect("equal")
    ax.pcolormesh(xs, ys, np.ma.masked_invalid(w), cmap="RdBu_r", vmin=-0.2, vmax=0.2, shading="nearest", alpha=0.6)
    ax.quiver(xs[::2], ys[::2], Gs[::2, ::2, 0], Gs[::2, ::2, 1], scale=4, width=0.004, color="0.15")
    for cc_, val_, col_, nm_ in ((cW, wallS, "#c2452d", "wall-side lobe"), (cB, bulkS, "#1f5fa8", "bulk-side lobe")):
        ax.add_patch(plt.Circle(cc_, 1.5, fill=False, ec=col_, lw=2.2, zorder=9))
        ax.text(cc_[0], cc_[1] + 1.9, f"{val_:+.1f}", color=col_, fontsize=10, fontweight="bold", ha="center", va="bottom", zorder=10,
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.75))
    for rr_, ls_ in ((2, "-"), (4, ":"), (8, ":")): ax.add_patch(plt.Circle(p, rr_, fill=False, ec="#e08a1e" if rr_ == 2 else "#1f5fa8", lw=1, ls=ls_, zorder=6))
    ax.plot(np.r_[W[:, 0], W[0, 0]], np.r_[W[:, 1], W[0, 1]], "k-", lw=3)
    for fd in fr[m]:
        if abs(fd[0] - cen[0]) > L or abs(fd[1] - cen[1]) > L: continue
        if fd[2] > 0:
            hx, hy = -np.cos(fd[3]), -np.sin(fd[3]); ax.plot(fd[0], fd[1], "o", ms=7, color="#2a7a3b", mec="k", zorder=5)
            ax.annotate("", (fd[0] + 1.6 * hx, fd[1] + 1.6 * hy), (fd[0], fd[1]), arrowprops=dict(arrowstyle="-|>", color="#2a7a3b", lw=2), zorder=5)
        else:
            ax.plot(fd[0], fd[1], "^", ms=8, color="#7a3b8f", mec="k", zorder=5)
    tr = np.array([ch[mm] for mm in ms if mm <= m]); ax.plot(tr[:, 0], tr[:, 1], "-", color="#e08a1e", lw=2, zorder=6)
    ax.plot(p[0], p[1], "o", ms=13, mfc="none", mec="#e08a1e", mew=3, zorder=7)
    hx, hy = -np.cos(ang), -np.sin(ang); ax.annotate("", (p[0] + 2.6 * hx, p[1] + 2.6 * hy), p, arrowprops=dict(arrowstyle="-|>", color="#e08a1e", lw=3), zorder=8)
    wp = W[j]; ax.annotate("", (wp[0] + 3 * nin.real, wp[1] + 3 * nin.imag), wp, arrowprops=dict(arrowstyle="-|>", color="k", lw=1.5, ls="--"), zorder=8)
    ax.set_xlim(cen[0] - L, cen[0] + L); ax.set_ylim(cen[1] - L, cen[1] + L); ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f"t − t_dep = {t[m] - t[n]:+.0f} t.u.    cos τ = {ct:+.2f}    depth = {d:.1f}", fontsize=11)
    H_ = np.array(hist); xl = (t[ms[0]] - t[n] - 2, t[ms[-1]] - t[n] + 2)
    a1 = fig.add_axes([0.70, 0.71, 0.28, 0.24]); a1.plot(H_[:, 0], H_[:, 1], "o", ms=3, color="#e08a1e", alpha=0.5); sm = np.convolve(H_[:, 1], np.ones(3) / 3, mode="same") if len(H_) > 2 else H_[:, 1]; a1.plot(H_[1:-1, 0], sm[1:-1], "-", color="#e08a1e", lw=2); a1.set_xlim(*xl); a1.set_ylim(-1.05, 1.05)
    a1.set_ylabel("cos τ  (1 = head into the bulk)", fontsize=9); a1.axvline(0, color="0.6", lw=0.8, ls=":")
    a2 = fig.add_axes([0.70, 0.40, 0.28, 0.24]); a2.plot(H_[:, 0], H_[:, 2], "-o", ms=3, color="#c2452d", label="at the defect (r<2)")
    a2.plot(H_[:, 0], H_[:, 3], "-s", ms=3, color="#1f5fa8", label="ring 4–8 cells"); a2.axhline(0, color="0.6", lw=0.8); a2.axvline(0, color="0.6", lw=0.8, ls=":")
    a2.set_xlim(*xl); a2.set_ylim(-0.4, 0.8); a2.set_ylabel("inward flow\n(cells/t.u.)", fontsize=9)
    a2.legend(fontsize=7, loc="upper left")
    a3 = fig.add_axes([0.70, 0.09, 0.28, 0.24])
    for jj_, col_, nm_ in ((4, "#c2452d", "wall-side lobe"), (5, "#1f5fa8", "bulk-side lobe")):
        a3.plot(H_[:, 0], H_[:, jj_], "o", ms=2.5, color=col_, alpha=0.4); smj = np.convolve(H_[:, jj_], np.ones(3) / 3, mode="same") if len(H_) > 2 else H_[:, jj_]
        a3.plot(H_[1:-1, 0], smj[1:-1], "-", color=col_, lw=2, label=nm_)
    a3.axhline(1, color="0.6", lw=0.8, ls="--"); a3.axhline(0, color="0.6", lw=0.8); a3.axvline(0, color="0.6", lw=0.8, ls=":")
    a3.set_xlim(*xl); a3.set_ylim(-1.5, 3.5); a3.set_xlabel("time relative to departure (t.u.)"); a3.set_ylabel("lobe strength\n(1 = free +1/2)", fontsize=9); a3.legend(fontsize=7, loc="upper left")
    fig.text(0.03, 0.005, "orange: escaping +1/2 (trail, head; circles r = 2, 4, 8)   red/blue rings: its wall-/bulk-side vortex lobes (strength vs a free +1/2)\n   green: other +1/2 (head)   purple: −1/2   dashed: inward wall normal   colour: vorticity (red = counter-clockwise)", fontsize=8)
    fig.savefig(f"frames/f{q:03d}.png"); plt.close(fig)
out = f"escape_{os.path.basename(run)}_ev{ie}.mp4"
subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", "4", "-i", "frames/f%03d.png", "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2", "-pix_fmt", "yuv420p", out], check=True)
print("wrote", out, "event", ie, "n", n, "t_dep", t[n])
