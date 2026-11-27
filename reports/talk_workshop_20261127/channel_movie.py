#!/usr/bin/env python3
"""Talk movie of an active channel with winding anchoring (user, 2026-11-26), same style as cavity_movie.py, fixed frame (no tight bbox -> no jitter).
Top: channel (flow vorticity smoothed over 2 cells, flow arrows, +1/2 red comets, -1/2 blue tripods), header with time, wall +1/2 count (within 4 cells of
either wall; 2m per wall if all at the walls) and bulk charge. Bottom: wall count and bulk charge against time over the window.
Usage: channel_movie.py RUNDIR T0 T1 STEP NAME [title] -> movies_fixed/movie_<NAME>.mp4 (20 fps)"""
import os, sys, json, subprocess, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
T = "/home/baloewe/Physics/Overcharged/reports/master/theory"; sys.path.insert(0, "/home/baloewe/Physics/Overcharged/runs/cavity"); sys.path.insert(0, T + "/energetics")
from edge_criteria import read_with_angle  # noqa
from selfflow import load_flow  # noqa
H = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(H, "movies_fixed")
run, T0, T1, STEP, NAME = sys.argv[1].rstrip("/"), float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), sys.argv[5]; TITLE = sys.argv[6] if len(sys.argv) > 6 else ""
J = json.load(open(run + "/input.json")); L, Hd = J["domain"][:2]; y0, y1 = J["BC"][2]["Q"][1], J["BC"][3]["Q"][1]; W = y1 - y0
mb, mt = J["BC"][2].get("customAnch", 0), -J["BC"][3].get("customAnch", 0); NW = 2 * abs(mb) + 2 * abs(mt)
t, fr = read_with_angle(run, np.array([0.0, 0.0]), 0.0); t = np.array(t)
rows, ft, idx = load_flow(run); ix = rows[:, 1].astype(int); iy = rows[:, 2].astype(int); xs, ys = np.arange(L) + 0.5, np.arange(Hd) + 0.5
sel = np.flatnonzero((t >= T0) & (t <= T1)); S = []
for n in sel:
    f = fr[n]
    if not len(f): S.append((t[n], np.nan, np.nan)); continue
    nearw = (np.abs(f[:, 1] - y0) < 4) | (np.abs(f[:, 1] - y1) < 4); w = nearw & (f[:, 2] > 0)
    S.append((t[n], w.sum(), 0.5 * np.sign(f[~nearw, 2]).sum()))
S = np.array(S, float); fdir = os.path.join(OUT, f"frames_{NAME}"); os.makedirs(fdir, exist_ok=True)
for f_ in os.listdir(fdir): os.remove(os.path.join(fdir, f_))


def vort(tt):
    ii = idx[ft[np.argmin(np.abs(ft - tt))]]; U = np.zeros((L, Hd)); V = np.zeros((L, Hd)); U[ix[ii], iy[ii]] = rows[ii, 4]; V[ix[ii], iy[ii]] = rows[ii, 5]
    Us, Vs = gaussian_filter(U, 2.0, mode="wrap"), gaussian_filter(V, 2.0, mode="wrap"); om = np.gradient(Vs, axis=0) - np.gradient(Us, axis=1)
    om[:, (ys < y0) | (ys > y1)] = np.nan; return Us, Vs, om


VM = np.nanmax([np.nanpercentile(np.abs(vort(tq)[2]), 95) for tq in np.linspace(T0, T1, 7)])
for kk, (n, st) in enumerate(list(zip(sel, S))[::STEP]):
    f = fr[n]; tt = t[n]; Us, Vs, om = vort(tt)
    fig = plt.figure(figsize=(14, 0.075 * W + 4.4))
    gs = fig.add_gridspec(2, 1, height_ratios=[0.075 * W + 0.6, 2.2], hspace=0.32, left=0.07, right=0.93, top=0.88, bottom=0.12)
    ax = fig.add_subplot(gs[0]); ax.pcolormesh(xs, ys, om.T, cmap="PuOr_r", vmin=-VM, vmax=VM, shading="nearest", alpha=0.85)
    ax.quiver(xs[2::5], ys[2::5], Us[2::5, 2::5].T, Vs[2::5, 2::5].T, color="0.3", scale=30, width=0.0012)
    ax.axhline(y0, color="k", lw=3); ax.axhline(y1, color="k", lw=3)
    for x, y, q, a in f:
        if q > 0:
            h = a + np.pi; ax.plot([x - 2.5 * np.cos(h), x], [y - 2.5 * np.sin(h), y], color="#c2452d", lw=2.4); ax.plot(x, y, "o", ms=5.5, color="#c2452d", mec="w", mew=0.6)
        else:
            for j in range(3): b = a + 2 * np.pi * j / 3; ax.plot([x, x + 1.8 * np.cos(b)], [y, y + 1.8 * np.sin(b)], color="#1f5fa8", lw=2)
            ax.plot(x, y, "o", ms=4.5, color="#1f5fa8", mec="w", mew=0.6)
    ax.set_xlim(0, L); ax.set_ylim(y0 - 1, y1 + 1); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title((TITLE + "\n" if TITLE else "") + f"t = {tt:.0f} t.u.   wall $+1/2$: {int(st[1])} (of {NW} if all at the walls)   bulk charge {st[2]:+.1f}", fontsize=11)
    a2 = fig.add_subplot(gs[1]); a2.plot(S[:, 0], S[:, 1], color="#c2452d", lw=0.8, alpha=0.8); a2.axhline(NW, color="#c2452d", ls=":", lw=0.8)
    a3 = a2.twinx(); a3.plot(S[:, 0], S[:, 2], color="#1f5fa8", lw=0.8, alpha=0.8); a3.set_ylabel("Bulk charge", color="#1f5fa8")
    a2.axvline(tt, color="k", lw=1); a2.set_xlim(T0, T1); a2.set_xlabel("Time (t.u.)"); a2.set_ylabel("Wall $+1/2$", color="#c2452d")
    fig.savefig(os.path.join(fdir, f"f{kk:04d}.png"), dpi=100, facecolor="white"); plt.close(fig)
mp4 = os.path.join(OUT, f"movie_{NAME}.mp4")
subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", "20", "-i", os.path.join(fdir, "f%04d.png"), "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2", "-pix_fmt", "yuv420p", mp4], check=True)
print(mp4, len(sel[::STEP]), "frames")
