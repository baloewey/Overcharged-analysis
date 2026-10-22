#!/usr/bin/env python3
"""Defects-only movie of a mixed-winding channel (user, 2026-10-22): defects, their trajectories, mid-line crossings, and the charge carried across.
Top: the channel (walls grey), +1/2 red comets (dot = head), -1/2 blue tripods, fading trajectories (last TAIL frames), dashed mid-line; a defect that
crosses the mid-line is ringed for 10 frames (green: carries net charge from the bottom wall towards the top, i.e. +1/2 up or -1/2 down; purple: the
reverse). Bottom: cumulative net charge carried across the mid-line upwards (black), and its parts: +1/2 crossings (red, up minus down) and -1/2
crossings (blue, down minus up, i.e. their contribution to the upward charge), against time; marker = movie time. Wall charges taken from the anchoring
(customAnch; bottom +m means a +1/2 layer). Every defect frame (2 t.u.), FPS frames per second.
Usage: crossing_movie.py RUN [FPS=40] [TAIL=25] -> movie_crossing_<name>.mp4 (in reports/master/theory/channel/crossing_movies)"""
import os, sys, json, subprocess, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, "/home/baloewe/Physics/Overcharged/runs/cavity")
from edge_criteria import read_with_angle  # noqa
H = os.path.join(os.path.dirname(os.path.abspath(__file__)), "crossing_movies")
def track(run):
    J = json.load(open(run + "/input.json")); L = J["domain"][0]; y0, y1 = J["BC"][2]["Q"][1], J["BC"][3]["Q"][1]; ym = 0.5 * (y0 + y1)
    t, fr = read_with_angle(run, np.array([0.0, 0.0]), 0.0); ids, nid, prev, previd = [], 0, None, None
    cum = np.zeros((len(fr), 3)); ev = []                       # columns: net upward charge, +1/2 part, -1/2 part
    for n, f in enumerate(fr):
        idn = np.full(len(f), -1); c = cum[n - 1].copy() if n else np.zeros(3)
        if prev is not None and len(prev) and len(f):
            used = set()
            for a, p in enumerate(prev):
                s = np.flatnonzero(np.sign(f[:, 2]) == np.sign(p[2]))
                if not len(s): continue
                dx = (f[s, 0] - p[0] + L / 2) % L - L / 2; d = np.hypot(dx, f[s, 1] - p[1]); j = int(np.argmin(d)); b = s[j]
                if d[j] < 2.5 and b not in used:
                    used.add(b); idn[b] = previd[a]; q = 0.5 * np.sign(p[2])
                    if p[1] < ym <= f[b, 1]: c[0] += q; c[1 if q > 0 else 2] += q; ev.append((n, b, +1 if q > 0 else -1))
                    elif f[b, 1] < ym <= p[1]: c[0] -= q; c[1 if q > 0 else 2] -= q; ev.append((n, b, -1 if q > 0 else +1))
        new = idn < 0; idn[new] = nid + np.arange(new.sum()); nid += new.sum(); ids.append(idn); cum[n] = c; prev, previd = f, idn
    return J, np.array(t), fr, ids, cum, ev
def main(run, FPS=40, TAIL=25):
    J, t, fr, ids, cum, ev = track(run); L, Hd = J["domain"][:2]; y0, y1 = J["BC"][2]["Q"][1], J["BC"][3]["Q"][1]; ym = 0.5 * (y0 + y1)
    mb, mt = J["BC"][2].get("customAnch", 0), J["BC"][3].get("customAnch", 0)
    name = os.path.basename(run.rstrip("/")); out = os.path.join(H, f"frames_{name}"); os.makedirs(out, exist_ok=True)
    for f_ in os.listdir(out): os.remove(os.path.join(out, f_))
    hist = {}; ring = []                                         # id -> list of (n, x, y)
    evn = {}
    for n, b, s in ev: evn.setdefault(n, []).append((b, s))
    hgt = 8.1 * (y1 - y0 + 2) / L                                  # channel panel height (in) at equal aspect over 8.1 in of width
    fig = plt.figure(figsize=(9.0, hgt + 3.1)); g = fig.add_gridspec(2, 1, height_ratios=[hgt, 2.0], hspace=1.8 / (hgt + 2.0), left=0.09, right=0.98,
                                                                   top=1 - 0.35 / (hgt + 3.1), bottom=0.55 / (hgt + 3.1))
    ax = fig.add_subplot(g[0]); ay = fig.add_subplot(g[1])
    ay.plot(t, cum[:, 0], "k-", lw=1.4, label="Net charge carried up"); ay.plot(t, cum[:, 1], "-", color="#c0392b", lw=1.0, label=r"Carried by $+\frac{1}{2}$")
    ay.plot(t, cum[:, 2], "-", color="#2e6fb5", lw=1.0, label=r"Carried by $-\frac{1}{2}$ ($-\frac{1}{2}$ moving up counts $-\frac{1}{2}$)"); ay.axhline(0, color="0.8", lw=0.6)
    ay.set_xlabel("Time (t.u.)"); ay.set_ylabel("Charge across\nthe mid-line"); ay.legend(loc="upper left", fontsize=8, ncol=3); ay.set_xlim(t[0], t[-1])
    mk, = ay.plot([], [], "o", color="k", ms=6); vl = ay.axvline(t[0], color="0.5", lw=0.8)
    bl = "+" if mb > 0 else "−"; tl_ = "+" if mt < 0 else "−"   # top normal points down: customAnch sign flips the layer
    for n in range(len(fr)):
        f = fr[n]
        for (x, y, q, a), i_ in zip(f, ids[n]): hist.setdefault(int(i_), []).append((n, x, y))
        for b, s in evn.get(n, []): ring.append((n, int(ids[n][b]), s))
        ring = [r for r in ring if n - r[0] < 10]
        ax.clear(); ax.axhspan(-1, y0, color="0.75"); ax.axhspan(y1, Hd + 1, color="0.75"); ax.axhline(ym, color="0.4", ls="--", lw=0.8)
        live = {int(i_) for i_ in ids[n]}
        for i_ in live:
            h = [p for p in hist[i_] if n - p[0] < TAIL]
            if len(h) < 2: continue
            xs = np.array([p[1] for p in h]); ys = np.array([p[2] for p in h]); col = "#c0392b" if np.sign(f[list(ids[n]).index(i_), 2]) > 0 else "#2e6fb5"
            brk = np.flatnonzero(np.abs(np.diff(xs)) > L / 2) + 1
            for k_, (sx, sy) in enumerate(zip(np.split(xs, brk), np.split(ys, brk))): ax.plot(sx, sy, "-", color=col, lw=1.0, alpha=0.45)
        for x, y, q, a in f:
            if q > 0:
                h_ = a + np.pi; ax.plot([x - 3 * np.cos(h_), x], [y - 3 * np.sin(h_), y], color="#c0392b", lw=1.6); ax.plot(x, y, "o", ms=3.2, color="#c0392b", mec="w", mew=0.4)
            else:
                for j in range(3): b_ = a + 2 * np.pi * j / 3; ax.plot([x, x + 2.2 * np.cos(b_)], [y, y + 2.2 * np.sin(b_)], color="#2e6fb5", lw=1.3)
        pos = {int(i_): (x, y) for (x, y, q, a), i_ in zip(f, ids[n])}
        for n0, i_, s in ring:
            if i_ in pos: ax.plot(*pos[i_], "o", ms=13, mfc="none", mec="#1a9850" if s > 0 else "#7b3294", mew=1.8)
        ax.set_xlim(0, L); ax.set_ylim(y0 - 1, y1 + 1); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f"{name}:  bottom wall {bl}{abs(mb)} (+1/2 layer)   top wall {tl_}{abs(mt)}   |   t = {t[n]:.0f} t.u.   net charge carried up = {cum[n, 0]:+.1f}", fontsize=9)
        mk.set_data([t[n]], [cum[n, 0]]); vl.set_xdata([t[n], t[n]])
        fig.savefig(os.path.join(out, f"f{n:05d}.png"), dpi=110)
    plt.close(fig); mp4 = os.path.join(H, f"movie_crossing_{name}.mp4")
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(FPS), "-i", os.path.join(out, "f%05d.png"), "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2", "-pix_fmt", "yuv420p", mp4], check=True)
    print(mp4, len(fr), "frames; net charge carried up", cum[-1, 0], "(+1/2 part", cum[-1, 1], ", -1/2 part", cum[-1, 2], ")", flush=True)
if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 40, int(sys.argv[3]) if len(sys.argv) > 3 else 25)
