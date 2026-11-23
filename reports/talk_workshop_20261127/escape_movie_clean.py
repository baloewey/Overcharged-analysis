#!/usr/bin/env python3
"""Clean escape movie for the workshop talk (slide 9). Same event and tracking as reports/master/theory/energetics/movie/escape_movie.py
(fineflow eps2.0_s2, event 4: flow and defects every 2 t.u.), drawn for a general audience: wall, light flow arrows, +1/2 as a red dot with its tail
(as in Fig. 1d), -1/2 as a blue tripod, the escaping +1/2 in orange with its trail and an arrow along its self-propulsion (head) direction; big clock
relative to departure; one-line key. No side panels.
Usage: escape_movie_clean.py RUNDIR EVENT -> assets/S4_escape_clean.mp4 (+ poster)"""
import os, sys, subprocess
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from scipy.ndimage import gaussian_filter
from matplotlib.path import Path
E = "/home/baloewe/Physics/Overcharged/reports/master/theory/energetics"; sys.path.insert(0, E)
import energetics as En  # noqa
from selfflow import load_flow  # noqa
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "assets")
run = sys.argv[1].rstrip("/"); ie = int(sys.argv[2]); J, bc, k, t, fr, ev = En.events(run, max_events=999); Q = bc["Q"]; sg = 1 if k > 0 else -1
W, s, P = En.wall_polyline(bc); tree = cKDTree(W); Wc = W[:, 0] + 1j * W[:, 1]; Tn = np.roll(Wc, -1) - np.roll(Wc, 1); Nn = -1j * Tn / np.abs(Tn)
if np.mean((Nn * np.conj(Wc)).real) < 0: Nn = -Nn
rows, ft, idx = load_flow(run); X = rows[:, 1] + 0.5 - Q[0]; Y = rows[:, 2] + 0.5 - Q[1]
RED, BLUE, ORANGE = "#c0392b", "#1f5fa8", "#f08c00"


def flow(tt):
    ii = idx[ft[np.argmin(np.abs(ft - tt))]]; return X[ii], Y[ii], rows[ii, 4], rows[ii, 5]


def follow(e, back=30, fwd=12):                     # as escape_movie.py
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


e = ev[ie]; ch = follow(e); n = e["n"]; site = ch[n]; _, j0 = tree.query(site); nin0 = -Nn[j0]
cen = site + 6 * np.array([nin0.real, nin0.imag]); L = 15
fdir = os.path.join(OUT, "escape_frames"); os.makedirs(fdir, exist_ok=True); [os.remove(os.path.join(fdir, f)) for f in os.listdir(fdir)]
ms = sorted(m for m in ch if n - 30 <= m <= n + 12)
outside = ~Path(W).contains_points(np.c_[np.meshgrid(np.linspace(cen[0] - L, cen[0] + L, 300), np.linspace(cen[1] - L, cen[1] + L, 300))[0].ravel(),
                                          np.meshgrid(np.linspace(cen[0] - L, cen[0] + L, 300), np.linspace(cen[1] - L, cen[1] + L, 300))[1].ravel()]).reshape(300, 300)
for q, m in enumerate(ms):
    x, y, vx, vy = flow(t[m]); win = (abs(x - cen[0]) < L + 1) & (abs(y - cen[1]) < L + 1)
    xs, ys = np.unique(x[win]), np.unique(y[win]); G = np.full((len(ys), len(xs), 2), np.nan)
    G[np.searchsorted(ys, y[win]), np.searchsorted(xs, x[win])] = np.c_[vx[win], vy[win]]
    XX, YY = np.meshgrid(xs, ys); inside = Path(W).contains_points(np.c_[XX.ravel(), YY.ravel()]).reshape(XX.shape) & np.isfinite(G[..., 0])
    Gs = np.zeros_like(G); wt = gaussian_filter(inside.astype(float), 1.2) + 1e-9
    for c_ in (0, 1): Gs[..., c_] = gaussian_filter(np.where(inside, G[..., c_], 0), 1.2) / wt
    Gs[~inside] = np.nan
    fig = plt.figure(figsize=(8, 7.4), dpi=120); ax = fig.add_axes([0.02, 0.12, 0.96, 0.79]); ax.set_aspect("equal")
    ax.imshow(np.where(outside, 1.0, np.nan), extent=(cen[0] - L, cen[0] + L, cen[1] - L, cen[1] + L), origin="lower", cmap="Greys", vmin=0, vmax=4, zorder=0)
    ax.quiver(xs[::2], ys[::2], Gs[::2, ::2, 0], Gs[::2, ::2, 1], scale=7.0, width=0.003, color="0.72", zorder=1)
    ax.plot(np.r_[W[:, 0], W[0, 0]], np.r_[W[:, 1], W[0, 1]], "k-", lw=4, zorder=2)
    p = ch[m]; f = fr[m]; me = int(np.argmin(np.hypot(*(f[:, :2] - p).T)))
    for i, fd in enumerate(f):
        if abs(fd[0] - cen[0]) > L or abs(fd[1] - cen[1]) > L or i == me: continue
        if fd[2] > 0:
            ax.plot([fd[0], fd[0] + 1.6 * np.cos(fd[3])], [fd[1], fd[1] + 1.6 * np.sin(fd[3])], "-", color=RED, lw=4, solid_capstyle="round", zorder=4)
            ax.plot(fd[0], fd[1], "o", ms=10, color=RED, mec="white", mew=1.5, zorder=5)
        else:
            for a_ in (0, 2 * np.pi / 3, -2 * np.pi / 3):
                ax.plot([fd[0], fd[0] + 1.2 * np.cos(fd[3] + a_)], [fd[1], fd[1] + 1.2 * np.sin(fd[3] + a_)], "-", color=BLUE, lw=3.5, solid_capstyle="round", zorder=4)
    tr = np.array([ch[mm] for mm in ms if mm <= m])
    if len(tr) > 2: tr = np.r_[tr[:1], (tr[:-2] + tr[1:-1] + tr[2:]) / 3, tr[-1:]]      # 3-frame running mean: grid jitter of the tracked core
    ax.plot(tr[:, 0], tr[:, 1], "-", color=ORANGE, lw=2.5, alpha=0.7, zorder=6)
    ang = f[me, 3]; ax.plot([p[0], p[0] + 1.8 * np.cos(ang)], [p[1], p[1] + 1.8 * np.sin(ang)], "-", color=ORANGE, lw=5, solid_capstyle="round", zorder=7)
    ax.plot(p[0], p[1], "o", ms=14, color=ORANGE, mec="k", mew=1.5, zorder=8)
    ax.annotate("", (p[0] - 3.2 * np.cos(ang), p[1] - 3.2 * np.sin(ang)), p, arrowprops=dict(arrowstyle="-|>", color=ORANGE, lw=3, mutation_scale=22), zorder=9)
    ax.set_xlim(cen[0] - L, cen[0] + L); ax.set_ylim(cen[1] - L, cen[1] + L); ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(False)
    dt_ = t[m] - t[n]; lab = "leaves the wall" if dt_ > 0 else ("departure" if dt_ == 0 else "at the wall")
    fig.text(0.03, 0.94, f"t = {dt_:+.0f} t.u.", fontsize=24, fontweight="bold", color="k"); fig.text(0.97, 0.94, lab, fontsize=20, ha="right", color=ORANGE)
    fig.text(0.5, 0.015, "Orange: escaping +½ (arrow: its self-propulsion)    Grey arrows: flow\nRed: other +½ (dot + tail)    Blue: −½",
             fontsize=14, linespacing=1.4, ha="center", color="0.25")
    fig.savefig(os.path.join(fdir, f"f{q:03d}.png"), facecolor="white"); plt.close(fig)
for extra in range(6):                              # hold the last frame ~1 s
    os.link(os.path.join(fdir, f"f{len(ms) - 1:03d}.png"), os.path.join(fdir, f"f{len(ms) + extra:03d}.png"))
out = os.path.join(OUT, "S4_escape_clean.mp4")
subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", "6", "-i", os.path.join(fdir, "f%03d.png"), "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2",
                "-pix_fmt", "yuv420p", "-c:v", "libx264", out], check=True)
subprocess.run(["cp", os.path.join(fdir, "f000.png"), os.path.join(OUT, "S4_escape_clean_poster.png")], check=True)
print("wrote", out, "frames", len(ms), "event", ie, "t_dep", t[n])
