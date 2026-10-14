#!/usr/bin/env python3
"""Movie of a winding channel (user, 2026-10-14: mixed channel, +6 wall facing a -6 wall). Background: flow vorticity (flow every 2 t.u., smoothed over 2 cells; colour scale fixed at the 95th percentile of the first frame);
arrows: flow (every 4th cell); +1/2 red comets (head = tail angle + pi), -1/2 blue tripods; top-left text: time, wall counts, bulk charge, cumulative mid-line crossings.
Usage: channel_movie.py RUNDIR [T_START T_LEN] -> frames_<name>/ + movie_<name>.mp4 (ffmpeg, 15 fps)"""
import os, sys, json, subprocess, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, "/home/baloewe/Physics/Overcharged/runs/cavity"); sys.path.insert(0, os.path.join(H, "..", "energetics"))
from edge_criteria import read_with_angle  # noqa
from selfflow import load_flow  # noqa
run = sys.argv[1].rstrip("/"); name = os.path.basename(run); T0 = float(sys.argv[2]) if len(sys.argv) > 2 else 3000.0; TL = float(sys.argv[3]) if len(sys.argv) > 3 else 600.0
J = json.load(open(run + "/input.json")); L, Hd = J["domain"][:2]; y0, y1 = J["BC"][2]["Q"][1], J["BC"][3]["Q"][1]; ym = 0.5 * (y0 + y1)
cb, ct = J["BC"][2].get("customAnch", 0), J["BC"][3].get("customAnch", 0)
t, fr = read_with_angle(run, np.array([0.0, 0.0]), 0.0); t = np.array(t)
rows, ft, idx = load_flow(run); ix = rows[:, 1].astype(int); iy = rows[:, 2].astype(int)
out = os.path.join(H, f"frames_{name}"); os.makedirs(out, exist_ok=True)
for f_ in os.listdir(out): os.remove(os.path.join(out, f_))
sel = np.flatnonzero((t >= T0) & (t < T0 + TL)); cross = 0.0; prev = None; k = 0
xs, ys = np.arange(L) + 0.5, np.arange(Hd) + 0.5
for n in sel:
    f = fr[n]; tt = t[n]; ii = idx[ft[np.argmin(np.abs(ft - tt))]]
    U = np.zeros((L, Hd)); V = np.zeros((L, Hd)); U[ix[ii], iy[ii]] = rows[ii, 4]; V[ix[ii], iy[ii]] = rows[ii, 5]
    Us, Vs = gaussian_filter(U, 2.0, mode="wrap"), gaussian_filter(V, 2.0, mode="wrap"); om = np.gradient(Vs, axis=0) - np.gradient(Us, axis=1)
    om[:, (ys < y0) | (ys > y1)] = np.nan
    if prev is not None and len(prev) and len(f):     # mid-line crossings (nearest same-sign link < 2.5 cells, periodic in x)
        for p in prev:
            s = np.flatnonzero(np.sign(f[:, 2]) == np.sign(p[2]))
            if not len(s): continue
            dx = (f[s, 0] - p[0] + L / 2) % L - L / 2; d = np.hypot(dx, f[s, 1] - p[1]); j = np.argmin(d)
            if d[j] < 2.5:
                if p[1] < ym <= f[s[j], 1]: cross += 0.5 * np.sign(p[2])
                elif f[s[j], 1] < ym <= p[1]: cross -= 0.5 * np.sign(p[2])
    prev = f
    fig, ax = plt.subplots(figsize=(14, 3.9))
    VM = globals().setdefault("VM", np.nanpercentile(np.abs(om), 95)); ax.pcolormesh(xs, ys, om.T, cmap="PuOr_r", vmin=-VM, vmax=VM, shading="nearest", alpha=0.85)
    ax.quiver(xs[2::5], ys[2::5], Us[2::5, 2::5].T, Vs[2::5, 2::5].T, color="0.3", scale=30, width=0.0012)
    ax.axhline(y0, color="k", lw=2); ax.axhline(y1, color="k", lw=2); ax.axhline(ym, color="0.5", lw=0.6, ls=":")
    if len(f):
        P_ = f[f[:, 2] > 0]; M_ = f[f[:, 2] < 0]
        for x, y, q, a in P_:
            h = a + np.pi; ax.plot([x - 2.5 * np.cos(h), x], [y - 2.5 * np.sin(h), y], color="#c2452d", lw=2.4); ax.plot(x, y, "o", ms=5.5, color="#c2452d", mec="w", mew=0.6)
        for x, y, q, a in M_:
            for j in range(3): b = a + 2 * np.pi * j / 3; ax.plot([x, x + 1.8 * np.cos(b)], [y, y + 1.8 * np.sin(b)], color="#1f5fa8", lw=2)
            ax.plot(x, y, "o", ms=4.5, color="#1f5fa8", mec="w", mew=0.6)
        d = np.minimum(f[:, 1] - y0, y1 - f[:, 1]); bot = f[:, 1] - y0 < 4; top = y1 - f[:, 1] < 4; bulk = d >= 4
        txt = (f"t = {tt:.0f} t.u.   bottom wall (charge {cb:+d}): {int(((f[:,2]>0)&bot).sum())} +1/2, {int(((f[:,2]<0)&bot).sum())} −1/2   |   "
               f"top wall (charge {-ct:+d}): {int(((f[:,2]>0)&top).sum())} +1/2, {int(((f[:,2]<0)&top).sum())} −1/2   |   bulk charge {0.5*np.sign(f[bulk,2]).sum():+.1f}   |   net mid-line crossings {cross:+.1f}")
        ax.set_title(txt, fontsize=9, loc="left")
    ax.set_xlim(0, L); ax.set_ylim(y0 - 1, y1 + 1); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    fig.tight_layout(); fig.savefig(os.path.join(out, f"f{k:04d}.png"), dpi=110); plt.close(fig); k += 1
mp4 = os.path.join(H, f"movie_{name}.mp4")
subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", "15", "-i", os.path.join(out, "f%04d.png"), "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2", "-pix_fmt", "yuv420p", mp4], check=True)
print(mp4, k, "frames")
