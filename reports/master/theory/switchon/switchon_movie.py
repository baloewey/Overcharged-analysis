#!/usr/bin/env python3
"""Movies of the layer forming after activity is switched on (runs/switchon; 2026-10-22). Frames every 10 t.u. (flow output) from the switch to
+1500 t.u. Background: flow speed |v| (cell velocity, 3x3 smoothed), arrows: flow (every 4th cell); wall +1/2 / bulk +1/2 as comets (dot + tail along
the radial arm), -1/2 as triangles; ticks outside the wall: homeotropic sites (grey) and the lock-selected sites for the final chirality (red,
shifted by tau/10 with tau = +-63.5 deg). Header: time since the switch, wall +1/2 count, circulation.
Usage: switchon_movie.py single RUN | montage   -> movie/<RUN>.mp4 or movie/montage.mp4"""
import sys, os, gzip, json, subprocess, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from multiprocessing import Pool
from scipy.ndimage import uniform_filter
sys.path.insert(0, "/home/baloewe/Physics/Overcharged/runs/cavity"); from edge_criteria import read_with_angle
RD = "/home/baloewe/Physics/Overcharged/runs/switchon"; RUNS = ["align45", "alignX", "alignY_s1", "alignY_s2", "rand_s1", "rand_s2", "rand_s3", "rand_s4"]
T0, T1, DT = 3000.0, 4500.0, 10.0
VMAX, QSC = 0.8, 15
def load(run):
    d = RD + "/" + run; bc = json.load(open(d + "/input.json"))["BC"][0]; R = bc["R"]; c = np.array(bc["Q"][:2]); L = 67
    F = {}; cur, t = None, None
    with gzip.open(d + "/flowfield.dat.gz", "rt") as f:
        for line in f:
            p = line.split()
            if len(p) != 7 or not p[0][0].isdigit(): continue
            tt = float(p[0])
            if tt > T1 + 1: break
            if tt != t:
                if cur is not None: F[t] = cur
                t = tt; cur = np.zeros((L, L, 2))
            cur[int(p[1]), int(p[2])] = (float(p[4]), float(p[5]))
    F[t] = cur
    td, fr = read_with_angle(d, c, 0.0); td = np.array(td)
    circ_sign = np.sign(np.mean([np.sum((-v[..., 0] * (np.arange(L)[None, :] + .5 - c[1]) + v[..., 1] * (np.arange(L)[:, None] + .5 - c[0]))) for tt, v in F.items() if tt > T1 - 300]))
    return dict(R=R, c=c, F=F, td=td, fr=fr, sg=circ_sign)
def draw(a, D, t, small=False):
    R, c = D["R"], D["c"]; v = D["F"][t]; L = v.shape[0]; xs = np.arange(L) + 0.5 - c[0]; X, Y = np.meshgrid(xs, np.arange(L) + 0.5 - c[1], indexing="ij")
    inside = np.hypot(X, Y) < R; sp = np.hypot(uniform_filter(v[..., 0], 3), uniform_filter(v[..., 1], 3)); sp[~inside] = np.nan
    a.pcolormesh(X, Y, sp, cmap="Blues", vmin=0, vmax=VMAX, shading="nearest")
    s = slice(1, None, 4) if not small else slice(2, None, 6); vs = np.where(inside[..., None], v, 0); ux, uy = uniform_filter(vs[..., 0], 3), uniform_filter(vs[..., 1], 3); m = (np.hypot(X, Y) < R - 2)[s, s]
    a.quiver(X[s, s][m], Y[s, s][m], ux[s, s][m], uy[s, s][m], color="0.3", scale=QSC, width=0.004 if not small else 0.006)
    th = np.linspace(0, 2 * np.pi, 200); a.plot(R * np.cos(th), R * np.sin(th), "k", lw=1.5)
    for j in range(10):
        for ph, col in ((2 * np.pi * j / 10, "0.55"), (2 * np.pi * j / 10 + D["sg"] * np.radians(63.5) / 10, "C3")):
            a.plot([(R + 0.6) * np.cos(ph), (R + 2.4) * np.cos(ph)], [(R + 0.6) * np.sin(ph), (R + 2.4) * np.sin(ph)], color=col, lw=2 if not small else 1.5)
    i = int(np.argmin(np.abs(D["td"] - t))); f = D["fr"][i]; nw = 0
    for x, y, q, ang in f:
        wall = R - np.hypot(x, y) < 4
        if q > 0:
            nw += wall; col = "C3" if wall else "C1"; ms = 6 if not small else 4
            a.plot(x, y, "o", color=col, ms=ms, mec="k", mew=0.5); a.plot([x, x + 3.5 * np.cos(ang)], [y, y + 3.5 * np.sin(ang)], color=col, lw=2)
        else:
            a.plot(x, y, "^", color="C2", ms=6 if not small else 4, mec="k", mew=0.5)
    num = np.nansum((-v[..., 0] * Y + v[..., 1] * X) / np.maximum(np.hypot(X, Y), 0.5) * inside); den = np.nansum(np.hypot(v[..., 0], v[..., 1]) * inside)
    a.set_xlim(-R - 3, R + 3); a.set_ylim(-R - 3, R + 3); a.set_aspect("equal"); a.axis("off")
    return nw, num / max(den, 1e-12)
def frame_single(args):
    run, k, t = args; D = DATA[run]; fig, a = plt.subplots(figsize=(7, 7.6)); nw, cc = draw(a, D, t)
    a.set_title(f"{run.replace('_', ' ')}: {t - T0:5.0f} t.u. after activity is switched on\nwall +1/2: {nw}   circulation: {cc:+.2f}", fontsize=14)
    fig.text(0.5, 0.02, "Red comets: wall +1/2; orange: bulk +1/2; green: -1/2.\nTicks: homeotropic sites (grey), lock-selected sites (red). Blue: flow speed; arrows: flow.", ha="center", fontsize=10)
    fig.savefig(f"movie/frames_{run}/f{k:04d}.png", dpi=110); plt.close(fig)
def frame_montage(args):
    k, t = args; fig, ax = plt.subplots(2, 4, figsize=(16, 8.6))
    for a, run in zip(ax.ravel(), RUNS):
        nw, cc = draw(a, DATA[run], t, small=True); a.set_title(f"{run.replace('_', ' ')}: wall {nw}, circ {cc:+.2f}", fontsize=11)
    fig.suptitle(f"Passive layer (12 wall +1/2, registry set by the initial director) -> active (0.07): {t - T0:5.0f} t.u. after the switch", fontsize=15)
    fig.text(0.5, 0.01, "Red comets: wall +1/2; orange: bulk +1/2; green: -1/2. Ticks outside the wall: homeotropic sites (grey) and lock-selected sites for the final chirality (red). Blue: flow speed; arrows: flow.", ha="center", fontsize=10)
    fig.tight_layout(rect=(0, 0.03, 1, 0.95)); fig.savefig(f"movie/frames_montage/f{k:04d}.png", dpi=90); plt.close(fig)
if __name__ == "__main__":
    mode = sys.argv[1]; runs = [sys.argv[2]] if mode == "single" else RUNS
    with Pool(len(runs)) as p: DATA = dict(zip(runs, p.map(load, runs)))
    times = sorted(t for t in DATA[runs[0]]["F"] if T0 <= t <= T1)
    tag = runs[0] if mode == "single" else "montage"; os.makedirs(f"movie/frames_{tag}", exist_ok=True)
    with Pool(8) as p:
        if mode == "single": p.map(frame_single, [(runs[0], k, t) for k, t in enumerate(times)])
        else: p.map(frame_montage, list(enumerate(times)))
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", "10", "-i", f"movie/frames_{tag}/f%04d.png", "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2",
                    "-c:v", "libx264", "-pix_fmt", "yuv420p", f"movie/{tag}.mp4"], check=True)
    print("done", tag, len(times), "frames")
