#!/usr/bin/env python3
"""Generic cavity movie for the SI (user, 2026-10-15). Background: flow vorticity (smoothed over 2 cells; colour scale from 7 frames across the window), flow arrows,
+1/2 red comets (dot = head), -1/2 blue tripods, the wall, and a running header: time, wall +1/2 count (within 4 cells), bulk net charge, tilt coherence of the layer.
Optional companion graph below: wall count and bulk charge against time over the whole window.
Usage: cavity_movie.py RUNDIR T0 T1 STEP NAME [title] -> movie_<NAME>.mp4 (20 fps)"""
import os, sys, json, subprocess, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from scipy.spatial import cKDTree
from matplotlib.path import Path
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, "/home/baloewe/Physics/Overcharged/runs/cavity"); sys.path.insert(0, os.path.join(H, "..", "energetics"))
from edge_criteria import read_with_angle  # noqa
from selfflow import load_flow  # noqa
import energetics as En  # noqa
run, T0, T1, STEP, NAME = sys.argv[1].rstrip("/"), float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), sys.argv[5]; TITLE = sys.argv[6] if len(sys.argv) > 6 else ""
J = json.load(open(run + "/input.json")); bc = J["BC"][0]; Q = np.array(bc["Q"][:2], float); k = bc.get("customAnch", 0)
W, s, P = En.wall_polyline(bc); tree = cKDTree(W); Wc = W[:, 0] + 1j * W[:, 1]; Tn = np.roll(Wc, -1) - np.roll(Wc, 1); Nn = -1j * Tn / np.abs(Tn)
if np.mean((Nn * np.conj(Wc)).real) < 0: Nn = -Nn
t, fr = read_with_angle(run, Q, 0.0); t = np.array(t)
rows, ft, idx = load_flow(run); X = rows[:, 1] + 0.5 - Q[0]; Y = rows[:, 2] + 0.5 - Q[1]
nx, ny = rows[:, 1].astype(int).max() + 1, rows[:, 2].astype(int).max() + 1; xs, ys = np.arange(nx) + 0.5 - Q[0], np.arange(ny) + 0.5 - Q[1]
XX, YY = np.meshgrid(xs, ys, indexing="ij"); inside = Path(W).contains_points(np.c_[XX.ravel(), YY.ravel()]).reshape(XX.shape)
sel = np.flatnonzero((t >= T0) & (t <= T1)); stats = []
for n in sel:
    f = fr[n]
    if not len(f): stats.append((t[n], np.nan, np.nan, np.nan)); continue
    d, j = tree.query(f[:, :2]); w = (f[:, 2] > 0) & (d < 4)
    tau = np.angle(np.exp(1j * (f[w, 3] - np.angle(Nn[j[w]])))); C = abs(np.mean(np.exp(1j * tau))) if w.any() else np.nan
    stats.append((t[n], w.sum(), 0.5 * np.sign(f[d >= 4, 2]).sum(), C))
S = np.array(stats, float); out = os.path.join(H, f"frames_{NAME}"); os.makedirs(out, exist_ok=True)
for f_ in os.listdir(out): os.remove(os.path.join(out, f_))
asp = np.ptp(W[:, 0]) / np.ptp(W[:, 1])
def vort(tt):
    ii = idx[ft[np.argmin(np.abs(ft - tt))]]; U = np.zeros((nx, ny)); V = np.zeros((nx, ny)); U[rows[ii, 1].astype(int), rows[ii, 2].astype(int)] = rows[ii, 4]; V[rows[ii, 1].astype(int), rows[ii, 2].astype(int)] = rows[ii, 5]
    w8 = gaussian_filter(inside.astype(float), 2.0) + 1e-9; Us = gaussian_filter(np.where(inside, U, 0), 2.0) / w8; Vs = gaussian_filter(np.where(inside, V, 0), 2.0) / w8
    om = np.gradient(Vs, axis=0) - np.gradient(Us, axis=1); om[~inside] = np.nan; return Us, Vs, om
VM = np.nanmax([np.nanpercentile(np.abs(vort(tq)[2]), 95) for tq in np.linspace(T0, T1, 7)])   # colour scale from frames across the window
for kk, (n, st) in enumerate(list(zip(sel, S))[::STEP]):
    f = fr[n]; tt = t[n]; Us, Vs, om = vort(tt)
    fig = plt.figure(figsize=(6.5 * max(asp, 1) ** 0.6 + 0.5, 8.6)); gs = fig.add_gridspec(2, 1, height_ratios=[3.2, 1.3], hspace=0.3)
    ax = fig.add_subplot(gs[0]); ax.pcolormesh(xs, ys, om.T, cmap="PuOr_r", vmin=-VM, vmax=VM, shading="nearest", alpha=0.85)
    m = inside[2::4, 2::4]; ax.quiver(XX[2::4, 2::4][m], YY[2::4, 2::4][m], Us[2::4, 2::4][m], Vs[2::4, 2::4][m], color="0.3", scale=12, width=0.003)
    ax.plot(np.r_[W[:, 0], W[0, 0]], np.r_[W[:, 1], W[0, 1]], "k-", lw=2)
    for x, y, q, a in f:
        if q > 0:
            h = a + np.pi; ax.plot([x - 2.5 * np.cos(h), x], [y - 2.5 * np.sin(h), y], color="#c2452d", lw=2.4); ax.plot(x, y, "o", ms=6, color="#c2452d", mec="w", mew=0.6)
        else:
            for j in range(3): b = a + 2 * np.pi * j / 3; ax.plot([x, x + 1.8 * np.cos(b)], [y, y + 1.8 * np.sin(b)], color="#1f5fa8", lw=2)
            ax.plot(x, y, "o", ms=5, color="#1f5fa8", mec="w", mew=0.6)
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_title((TITLE + "\n" if TITLE else "") + f"t = {tt:.0f} t.u.   wall $+1/2$: {int(st[1])} (of {2*abs(k)-2} sites)   bulk charge {st[2]:+.1f}   tilt coherence {st[3]:.2f}", fontsize=10)
    a2 = fig.add_subplot(gs[1]); a2.plot(S[:, 0], S[:, 1], color="#c2452d", lw=0.8, alpha=0.8, label="wall $+1/2$"); a2.axhline(2 * abs(k) - 2, color="#c2452d", ls=":", lw=0.8)
    a3 = a2.twinx(); a3.plot(S[:, 0], S[:, 2], color="#1f5fa8", lw=0.8, alpha=0.8); a3.set_ylabel("bulk charge", color="#1f5fa8")
    a2.axvline(tt, color="k", lw=1); a2.set_xlim(T0, T1); a2.set_xlabel("time (t.u.)"); a2.set_ylabel("wall $+1/2$", color="#c2452d")
    fig.savefig(os.path.join(out, f"f{kk:04d}.png"), dpi=100, bbox_inches="tight"); plt.close(fig)
mp4 = os.path.join(H, f"movie_{NAME}.mp4")
subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", "20", "-i", os.path.join(out, "f%04d.png"), "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", "-pix_fmt", "yuv420p", mp4], check=True)
print(mp4, len(sel[::STEP]), "frames")
