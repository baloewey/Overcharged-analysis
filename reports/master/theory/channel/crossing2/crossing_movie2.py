#!/usr/bin/env python3
"""Crossing movie with birth-to-death trajectories (user, 2026-10-22; replaces ../crossing_movie.py, whose nearest-neighbour count was biased by fast defects).
Tracks from tracker2.build (predictive Hungarian, gate 9 cells, 2-frame gaps). Top: +1/2 red comets, -1/2 blue tripods, the WHOLE trajectory of every
living defect since its birth (faint), dashed mid-line; a defect crossing the line is ringed for 10 frames (green: charge up, purple: charge down).
Bottom: tracking-free Q_up(t) - Q_up(0) (black, exact by topology) and the tracked count (crossings + pairs straddling the line; orange); marker = time.
Every frame (2 t.u.), FPS frames/s. Usage: crossing_movie2.py RUN [FPS=40] -> crossing2/movies/movie_crossing2_<name>.mp4"""
import os, sys, subprocess, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from tracker2 import build
H = os.path.join(os.path.dirname(os.path.abspath(__file__)), "movies")
def main(run, FPS=40):
    B = build(run, gate=9.0); t, L, ym, y0, y1, fr, J = B["t"], B["L"], B["ym"], B["y0"], B["y1"], B["fr"], B["J"]; Hd = J["domain"][1]
    mb, mt = J["BC"][2].get("customAnch", 0), J["BC"][3].get("customAnch", 0); name = os.path.basename(run.rstrip("/"))
    q = B["Qup"] - np.nanmean(B["Qup"][:25]); cum = B["cum"]
    byn = {}
    for k, tr in enumerate(B["tracks"]):
        for i, n in enumerate(tr["n"]): byn.setdefault(n, []).append((k, i))
    xev = {}
    for k, tr in enumerate(B["tracks"]):
        y = np.array(tr["xy"])[:, 1]; s = np.sign(y - ym)
        for i in np.flatnonzero(s[1:] != s[:-1]): xev.setdefault(tr["n"][i + 1], []).append((k, (1 if s[i + 1] > 0 else -1) * tr["sg"]))
    out = os.path.join(H, f"frames_{name}"); os.makedirs(out, exist_ok=True)
    for f_ in os.listdir(out): os.remove(os.path.join(out, f_))
    hgt = 8.1 * (y1 - y0 + 2) / L; fig = plt.figure(figsize=(9.0, hgt + 3.1))
    g = fig.add_gridspec(2, 1, height_ratios=[hgt, 2.0], hspace=1.8 / (hgt + 2.0), left=0.09, right=0.98, top=1 - 0.35 / (hgt + 3.1), bottom=0.55 / (hgt + 3.1))
    ax = fig.add_subplot(g[0]); ay = fig.add_subplot(g[1])
    ay.plot(t, q, "k-", lw=1.2, label=r"Tracking-free $Q_\mathrm{up}(t)-Q_\mathrm{up}(0)$"); ay.plot(t, cum, "-", color="#e6550d", lw=1.0, label="Tracked (crossings + straddling pairs)")
    ay.axhline(0, color="0.8", lw=0.6); ay.set_xlabel("Time (t.u.)"); ay.set_ylabel("Charge carried up\nacross the mid-line"); ay.legend(loc="upper left", fontsize=8, ncol=2); ay.set_xlim(t[0], t[-1])
    lo, hi = min(np.nanmin(q), cum.min()) - 1, max(np.nanmax(q), cum.max()) + 1; ay.set_ylim(min(lo, -6), max(hi, 6))
    mk, = ay.plot([], [], "o", color="k", ms=6); vl = ay.axvline(t[0], color="0.5", lw=0.8); ring = []
    bl = "+" if mb > 0 else "−"; tl_ = "+" if mt < 0 else "−"
    for n in range(len(fr)):
        ring = [r for r in ring if n - r[0] < 10] + [(n, k, s) for k, s in xev.get(n, [])]
        ax.clear(); ax.axhspan(-1, y0, color="0.75"); ax.axhspan(y1, Hd + 1, color="0.75"); ax.axhline(ym, color="0.4", ls="--", lw=0.8); pos = {}
        for k, i in byn.get(n, []):
            tr = B["tracks"][k]; xy = np.array(tr["xy"][:i + 1]); col = "#c0392b" if tr["sg"] > 0 else "#2e6fb5"; pos[k] = xy[-1]
            if len(xy) > 1:
                brk = np.flatnonzero(np.abs(np.diff(xy[:, 0])) > L / 2) + 1
                for seg in np.split(xy, brk): ax.plot(seg[:, 0], seg[:, 1], "-", color=col, lw=0.8, alpha=0.35)
        f = fr[n]
        for x, y, qq, a in f:
            if qq > 0:
                h_ = a + np.pi; ax.plot([x - 3 * np.cos(h_), x], [y - 3 * np.sin(h_), y], color="#c0392b", lw=1.6); ax.plot(x, y, "o", ms=3.2, color="#c0392b", mec="w", mew=0.4)
            else:
                for j in range(3): b_ = a + 2 * np.pi * j / 3; ax.plot([x, x + 2.2 * np.cos(b_)], [y, y + 2.2 * np.sin(b_)], color="#2e6fb5", lw=1.3)
        for n0, k, s in ring:
            if k in pos: ax.plot(*pos[k], "o", ms=13, mfc="none", mec="#1a9850" if s > 0 else "#7b3294", mew=1.8)
        ax.set_xlim(0, L); ax.set_ylim(y0 - 1, y1 + 1); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f"{name}: bottom {bl}{abs(mb)}, top {tl_}{abs(mt)}  |  t = {t[n]:.0f}  |  charge carried up: exact {q[n]:+.1f}, tracked {cum[n]:+.1f}", fontsize=9)
        mk.set_data([t[n]], [q[n]]); vl.set_xdata([t[n], t[n]]); fig.savefig(os.path.join(out, f"f{n:05d}.png"), dpi=110)
    plt.close(fig); mp4 = os.path.join(H, f"movie_crossing2_{name}.mp4")
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(FPS), "-i", os.path.join(out, "f%05d.png"), "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2", "-pix_fmt", "yuv420p", mp4], check=True)
    print(mp4, len(fr), "frames; exact", q[-1], "tracked", cum[-1], B["parts"], flush=True)
if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 40)
