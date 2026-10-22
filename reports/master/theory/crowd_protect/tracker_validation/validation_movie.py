#!/usr/bin/env python3
"""Tracker checking movies (user, 2026-10-22: 'movies that allow us to visually check if the trackers are working; draw the trajectories; not too fast').
Same frames, two panels: OLD = nearest-neighbour linker (mutual nearest neighbour within 3.5 cells, 2-frame gaps, no prediction; as exchange_count.py)
and NEW = predictive linker v2 (trackers.track_new: cost min(distance to prediction, distance to last position), Hungarian, gate 9, 2-frame gaps).
Each track has its own colour; its trajectory since birth is drawn (older parts fainter). +1/2 dots, -1/2 triangles. Green ring = track born this frame,
black cross = track ended. Red segment = a link longer than 4 cells (the risky ones). Orange star = an escape by the shape_runs definition (start in the
wall band, deeper than 6 cells 20 t.u. later) according to that linker, shown at the start frame and held for 20 frames. Grey band: 4-cell wall band;
dashed: 6-cell escape depth (cavities) or the mid-line (channels). Usage: validation_movie.py RUN T_START (negative = auto: 40 frames before the first escape NEW finds and OLD misses) N_FRAMES [FPS=6] -> movies/val_<name>.mp4"""
import json, os, subprocess, sys, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.path import Path
from trackers import track_new, track_nn
sys.path.insert(0, "../exchange"); import exchange_count as X  # noqa
H = os.path.join(os.path.dirname(os.path.abspath(__file__)), "movies"); os.makedirs(H, exist_ok=True)
PAL = plt.cm.tab20(np.linspace(0, 1, 20))
def escapes(T, dist, nf):
    out = {}
    for k, tr in enumerate(T):
        pos = dict(zip(tr["n"], tr["xy"]))
        for n, xy in pos.items():
            if tr["sg"] < 0 or n % 5 or n + 10 >= nf or dist(np.array([xy]))[0] >= 4: continue
            e = pos.get(n + 10)
            if e is not None and dist(np.array([e]))[0] > 6: out.setdefault(n, []).append((k, xy))
    return out
def main(run, t0, nfr, fps=6):
    kind, Q, dist, Lw, info = X.geometry(run); t, fr = X.read_with_angle(run, Q, 0.0); t = np.array(t); J = json.load(open(run + "/input.json"))
    L = J["domain"][0] if kind == "channel" else None
    if t0 < 0:                                                            # auto: start 40 frames before the first escape found by NEW but not by OLD
        En, Eo = escapes(track_new(fr, 9.0, L=L), dist, len(fr)), escapes(track_nn(fr, 3.5, L=L), dist, len(fr)); first = None
        for n0 in sorted(En):
            for k, xy in En[n0]:
                if not any(abs(n1 - n0) <= 3 and np.hypot(*(np.array(xy) - np.array(xy2))) < 4 for n1, v in Eo.items() for _, xy2 in v): first = n0; break
            if first is not None: break
        t0 = t[max((first or 0) - 40, 0)]; print("auto window starts at t =", t0, flush=True)
    i0 = int(np.argmin(np.abs(t - t0))); i1 = min(i0 + nfr, len(fr)); fr = fr[i0:i1]; t = t[i0:i1]; nf = len(fr)
    TR = {"OLD: nearest neighbour, 3.5 cells, no prediction": track_nn(fr, 3.5, L=L), "NEW: predictive + Hungarian (v2), gate 9": track_new(fr, 9.0, L=L)}
    ESC = {k: escapes(T, dist, nf) for k, T in TR.items()}; name = run.rstrip("/").split("runs/")[1].replace("/", "_")
    out = os.path.join(H, f"frames_{name}"); os.makedirs(out, exist_ok=True)
    for f_ in os.listdir(out): os.remove(os.path.join(out, f_))
    if kind == "channel":
        y0, y1 = J["BC"][2]["Q"][1], J["BC"][3]["Q"][1]; fig, axs = plt.subplots(2, 1, figsize=(15, 2 * (2.2 + 12.5 * (y1 - y0) / L)))
    else:
        bc = J["BC"][0]; W, s, P = X.wall_polyline(bc); asp_ = np.ptp(W[:, 1]) / np.ptp(W[:, 0]); fig, axs = plt.subplots(1, 2, figsize=(15, 7.3 * asp_ + 1.4))
        xs = np.linspace(W[:, 0].min() - 1, W[:, 0].max() + 1, 260); ys = np.linspace(W[:, 1].min() - 1, W[:, 1].max() + 1, 260); XX, YY = np.meshgrid(xs, ys)
        ins = Path(W).contains_points(np.c_[XX.ravel(), YY.ravel()]).reshape(XX.shape); DD = dist(np.c_[XX.ravel(), YY.ravel()]).reshape(XX.shape)
    fig.subplots_adjust(left=0.02, right=0.98, top=0.86 if kind != "channel" else 0.9, bottom=0.03, wspace=0.04, hspace=0.25)
    for n in range(nf):
        for ax, (lab, T) in zip(axs, TR.items()):
            ax.clear()
            if kind == "channel":
                ax.axhspan(y0 - 1, y0, color="0.6"); ax.axhspan(y1, y1 + 1, color="0.6"); ax.axhspan(y0, y0 + 4, color="0.92"); ax.axhspan(y1 - 4, y1, color="0.92")
                ax.axhline(0.5 * (y0 + y1), color="0.4", ls="--", lw=0.8); ax.set_xlim(0, L); ax.set_ylim(y0 - 1, y1 + 1)
            else:
                ax.contourf(XX, YY, np.where(ins, (DD < 4).astype(float), np.nan), levels=[-0.5, 0.5, 1.5], colors=["white", "0.9"])
                ax.contour(XX, YY, np.where(ins, DD, np.nan), levels=[6], colors="0.5", linestyles="--", linewidths=0.8)
                ax.plot(np.r_[W[:, 0], W[0, 0]], np.r_[W[:, 1], W[0, 1]], "k-", lw=1.2)
            for k, tr in enumerate(T):
                nn = np.array(tr["n"]); m = nn <= n
                if not m.any() or nn[0] > n or (nn[-1] < n - 2): continue
                xy = np.array(tr["xy"])[m]; c = PAL[k % 20]
                if len(xy) > 1:
                    seg = np.split(xy, np.flatnonzero(np.abs(np.diff(xy[:, 0])) > (L / 2 if L else 1e9)) + 1)
                    for sg_ in seg:
                        for a in range(1, len(sg_)):
                            age = n - nn[m][-len(sg_) + a] if len(sg_) == len(xy) else 0; al = max(0.15, 1 - (n - nn[m][a]) / 60) if len(seg) == 1 else 0.6
                            st = np.hypot(*(sg_[a] - sg_[a - 1])); ax.plot(sg_[a - 1:a + 1, 0], sg_[a - 1:a + 1, 1], "-", color="red" if st > 4 else c, lw=2.2 if st > 4 else 1.1, alpha=1 if st > 4 else al)
                if nn[m][-1] == n:
                    ax.plot(*xy[-1], "o" if tr["sg"] > 0 else "^", ms=6, color=c, mec="k", mew=0.5)
                    if nn[0] == n and n > 0: ax.plot(*xy[-1], "o", ms=13, mfc="none", mec="#1a9850", mew=2)
                if nn[-1] == n and n < nf - 3: ax.plot(*xy[-1], "x", ms=11, color="k", mew=2)
            for n0 in range(max(0, n - 20), n + 1):
                for k, xy in ESC[lab].get(n0, []): ax.plot(*xy, "*", ms=16, color="#ff7f00", mec="k", mew=0.6)
            ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
            nE = sum(len(v) for kk, v in ESC[lab].items() if kk <= n)
            ax.set_title(f"{lab}\nlive tracks {sum(1 for tr in T if tr['n'][0] <= n <= tr['n'][-1])}, escapes so far {nE}", fontsize=10)
        fig.suptitle(f"{name}   t = {t[n]:.0f} t.u.   (frame {n + 1}/{nf}, 2 t.u. per frame; {fps} frames/s)", fontsize=11)
        fig.savefig(os.path.join(out, f"f{n:04d}.png"), dpi=95)
    plt.close(fig); mp4 = os.path.join(H, f"val_{name}.mp4")
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(fps), "-i", os.path.join(out, "f%04d.png"), "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2", "-pix_fmt", "yuv420p", mp4], check=True)
    print(mp4, nf, "frames", {k: sum(len(v) for v in e.values()) for k, e in ESC.items()}, flush=True)
if __name__ == "__main__":
    main(sys.argv[1], float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]) if len(sys.argv) > 4 else 6)
