#!/usr/bin/env python3
"""Movie of the topological current in a mixed-winding channel, with companion graph (user, 2026-10-15).
Top: the channel (flow vorticity smoothed over 2 cells, flow arrows, +1/2 red comets, -1/2 blue tripods, wall charges labelled).
Bottom: accumulated charge transferred across the mid-line from the + wall towards the - wall, per unit length (net crossings / L), against time, for this run
and a symmetric control (+6/+6) of the same width, with the straight line of the mean current; a marker follows the movie time.
Usage: current_movie.py MIXED_RUN CONTROL_RUN [T_START T_END STEP] -> movie_current_<name>.mp4"""
import os, sys, json, subprocess, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, "/home/baloewe/Physics/Overcharged/runs/cavity"); sys.path.insert(0, os.path.join(H, "..", "energetics"))
from edge_criteria import read_with_angle  # noqa
from selfflow import load_flow  # noqa


def crossings(run):
    """cumulative net charge crossing the mid-line upwards (per unit length) at every defect frame"""
    J = json.load(open(run + "/input.json")); L = J["domain"][0]; y0, y1 = J["BC"][2]["Q"][1], J["BC"][3]["Q"][1]; ym = 0.5 * (y0 + y1)
    t, fr = read_with_angle(run, np.array([0.0, 0.0]), 0.0); c = 0.0; out = []; prev = None
    for f in fr:
        if prev is not None and len(prev) and len(f):
            for p in prev:
                s = np.flatnonzero(np.sign(f[:, 2]) == np.sign(p[2]))
                if not len(s): continue
                dx = (f[s, 0] - p[0] + L / 2) % L - L / 2; d = np.hypot(dx, f[s, 1] - p[1]); j = np.argmin(d)
                if d[j] < 2.5:
                    if p[1] < ym <= f[s[j], 1]: c += 0.5 * np.sign(p[2])
                    elif f[s[j], 1] < ym <= p[1]: c -= 0.5 * np.sign(p[2])
        out.append(c / L); prev = f
    return np.array(t), np.array(out), t, fr, J


if __name__ == "__main__":
    run, ctl = sys.argv[1].rstrip("/"), sys.argv[2].rstrip("/"); T0 = float(sys.argv[3]) if len(sys.argv) > 3 else 1000; T1 = float(sys.argv[4]) if len(sys.argv) > 4 else 4000
    STEP = int(sys.argv[5]) if len(sys.argv) > 5 else 3
    tm, cm, t, fr, J = crossings(run); tc, cc, _, _, Jc = crossings(ctl)
    L, Hd = J["domain"][:2]; y0, y1 = J["BC"][2]["Q"][1], J["BC"][3]["Q"][1]; ym = 0.5 * (y0 + y1); W = y1 - y0
    cb, ct = J["BC"][2].get("customAnch", 0), -J["BC"][3].get("customAnch", 0)
    rows, ft, idx = load_flow(run); ix = rows[:, 1].astype(int); iy = rows[:, 2].astype(int); xs, ys = np.arange(L) + 0.5, np.arange(Hd) + 0.5
    m0 = tm >= T0; cm0 = cm - np.interp(T0, tm, cm); cc0 = cc - np.interp(T0, tc, cc)
    sl = np.polyfit(tm[m0 & (tm <= T1)], cm0[m0 & (tm <= T1)], 1)[0]; slc = np.polyfit(tc[(tc >= T0) & (tc <= T1)], cc0[(tc >= T0) & (tc <= T1)], 1)[0]
    name = os.path.basename(run); out = os.path.join(H, f"frames_current_{name}"); os.makedirs(out, exist_ok=True)
    for f_ in os.listdir(out): os.remove(os.path.join(out, f_))
    sel = np.flatnonzero((t >= T0) & (t <= T1))[::STEP]; VM = None
    for k, n in enumerate(sel):
        f = fr[n]; tt = t[n]; ii = idx[ft[np.argmin(np.abs(ft - tt))]]
        U = np.zeros((L, Hd)); V = np.zeros((L, Hd)); U[ix[ii], iy[ii]] = rows[ii, 4]; V[ix[ii], iy[ii]] = rows[ii, 5]
        Us, Vs = gaussian_filter(U, 2.0, mode="wrap"), gaussian_filter(V, 2.0, mode="wrap"); om = np.gradient(Vs, axis=0) - np.gradient(Us, axis=1); om[:, (ys < y0) | (ys > y1)] = np.nan
        VM = VM or np.nanpercentile(np.abs(om), 95)
        fig = plt.figure(figsize=(14, 0.075 * W + 4.2)); gs = fig.add_gridspec(2, 1, height_ratios=[0.075 * W + 0.6, 2.4], hspace=0.28)
        ax = fig.add_subplot(gs[0]); ax.pcolormesh(xs, ys, om.T, cmap="PuOr_r", vmin=-VM, vmax=VM, shading="nearest", alpha=0.85)
        ax.quiver(xs[2::5], ys[2::5], Us[2::5, 2::5].T, Vs[2::5, 2::5].T, color="0.3", scale=30, width=0.0012)
        ax.axhline(y0, color="#c2452d", lw=3); ax.axhline(y1, color="#1f5fa8", lw=3); ax.axhline(ym, color="0.5", lw=0.6, ls=":")
        for x, y, q, a in f:
            if q > 0:
                h = a + np.pi; ax.plot([x - 2.5 * np.cos(h), x], [y - 2.5 * np.sin(h), y], color="#c2452d", lw=2.4); ax.plot(x, y, "o", ms=5.5, color="#c2452d", mec="w", mew=0.6)
            else:
                for j in range(3): b = a + 2 * np.pi * j / 3; ax.plot([x, x + 1.8 * np.cos(b)], [y, y + 1.8 * np.sin(b)], color="#1f5fa8", lw=2)
                ax.plot(x, y, "o", ms=4.5, color="#1f5fa8", mec="w", mew=0.6)
        ax.set_xlim(0, L); ax.set_ylim(y0 - 1, y1 + 1); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        ax.text(L - 1, y1 + 1.5, f"top wall: charge {ct:+d}", color="#1f5fa8", fontsize=9, va="bottom", ha="right")
        ax.text(L - 1, y0 - 1.5, f"bottom wall: charge {cb:+d}", color="#c2452d", fontsize=9, va="top", ha="right")
        ax.set_title(f"t = {tt:.0f} t.u.   (channel L = {L}, W = {W:.0f}; $+1/2$ red, $-1/2$ blue; background: vorticity)", fontsize=9.5, loc="left", pad=14)
        a2 = fig.add_subplot(gs[1])
        a2.plot(tc[tc >= T0], cc0[tc >= T0] * 1000, color="0.55", lw=1.2, label=f"symmetric control (+6/+6): {slc*1e3:+.3f}")
        a2.plot(tm[m0 & (tm <= tt)], cm0[m0 & (tm <= tt)] * 1000, color="#c2452d", lw=1.8, label=f"mixed (+6/$-$6): {sl*1e3:+.3f}")
        a2.plot(tm[m0], (tm[m0] - T0) * sl * 1000, color="#c2452d", lw=0.8, ls="--", alpha=0.6)
        a2.plot([tt], [np.interp(tt, tm, cm0) * 1000], "o", color="#c2452d", ms=6)
        a2.set_xlim(T0, T1); a2.set_xlabel("time (t.u.)"); a2.set_ylabel("charge transferred\nper unit length ($10^{-3}$)")
        a2.legend(fontsize=8.5, frameon=False, loc="upper left", title="Mean current (per unit length per 1000 t.u.)", title_fontsize=8.5)
        a2.axhline(0, color="0.8", lw=0.6); [a2.spines[s_].set_visible(False) for s_ in ("top", "right")]
        fig.savefig(os.path.join(out, f"f{k:04d}.png"), dpi=100, bbox_inches="tight"); plt.close(fig)
    mp4 = os.path.join(H, f"movie_current_{name}.mp4")
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", "20", "-i", os.path.join(out, "f%04d.png"), "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2,scale=trunc(iw/2)*2:trunc(ih/2)*2", "-pix_fmt", "yuv420p", mp4], check=True)
    print(mp4, len(sel), "frames; mean currents", sl * 1e3, slc * 1e3)
