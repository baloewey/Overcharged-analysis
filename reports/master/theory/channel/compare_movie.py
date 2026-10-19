#!/usr/bin/env python3
"""Forced vs unforced winding channels, side by side (user, 2026-10-19). Three stacked panels (same time window, L 192, W 24): unforced active channel, and
pressure-driven (body force g along x) at g = 0.001 and 0.003. Background: flow vorticity (smoothed 2 cells; one colour scale for all panels, 95th percentile of the
unforced first frame); arrows: flow; +1/2 red comets (head = tail angle + pi), -1/2 blue tripods. Panel title: case, live wall +1/2 / -1/2 counts, bulk charge and the
tilt coherence C of the bottom +1/2 layer. Usage: compare_movie.py sym|mixed [T0 TL] -> movie_forced_vs_unforced_<set>.mp4 (15 fps, one frame per 2 t.u.)"""
import os, sys, json, gzip, subprocess, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from multiprocessing import Pool
sys.path.insert(0, "/home/baloewe/Physics/Overcharged/runs/cavity"); from edge_criteria import read_with_angle  # noqa
H = os.path.dirname(os.path.abspath(__file__)); R_ = "/home/baloewe/Physics/Overcharged/runs"
SETS = {"sym": [("channel_current/W24_sym_s5", "Unforced"), ("channel_pressure/sym_g+0.001_s2", "Pressure-driven, g = 0.001"), ("channel_pressure/sym_g+0.003_s1", "Pressure-driven, g = 0.003")],
        "mixed": [("channel_current/W24_mixed_s3", "Unforced"), ("channel_pressure/mixed_g+0.001_s1", "Pressure-driven, g = 0.001"), ("channel_pressure/mixed_g+0.003_s1", "Pressure-driven, g = 0.003")]}
which = sys.argv[1]; T0 = float(sys.argv[2]) if len(sys.argv) > 2 else 4000.0; TL = float(sys.argv[3]) if len(sys.argv) > 3 else 400.0
def load(r):
    run = f"{R_}/{r}"; J = json.load(open(run + "/input.json")); L, Hd = J["domain"][:2]; y0, y1 = J["BC"][2]["Q"][1], J["BC"][3]["Q"][1]
    F = {}
    with gzip.open(run + "/flowfield.dat.gz", "rt") as fh:
        for l in fh:
            p = l.split()
            if len(p) != 7 or not p[0][0].isdigit(): continue
            tt = float(p[0])
            if tt < T0 - 2: continue
            if tt > T0 + TL + 2: break
            F.setdefault(tt, []).append([float(v) for v in p[1:6]])
    F = {k: np.array(v) for k, v in F.items()}
    t, fr = read_with_angle(run, np.array([0.0, 0.0]), 0.0); t = np.array(t); sel = np.flatnonzero((t >= T0) & (t < T0 + TL))
    return dict(L=L, Hd=Hd, y0=y0, y1=y1, cb=J["BC"][2].get("customAnch", 0), ct=-J["BC"][3].get("customAnch", 0), g=(J.get("grav") or [0])[0], F=F, t=t[sel], fr=[fr[i] for i in sel])
with Pool(3) as p: D = p.map(load, [r for r, _ in SETS[which]])
labels = [l for _, l in SETS[which]]
def field(d, tt):
    ft = np.array(sorted(d["F"])); a = d["F"][ft[np.argmin(np.abs(ft - tt))]]; U = np.zeros((d["L"], d["Hd"])); V = np.zeros_like(U)
    U[a[:, 0].astype(int), a[:, 1].astype(int)] = a[:, 3]; V[a[:, 0].astype(int), a[:, 1].astype(int)] = a[:, 4]
    Us, Vs = gaussian_filter(U, 2.0, mode="wrap"), gaussian_filter(V, 2.0, mode="wrap"); om = np.gradient(Vs, axis=0) - np.gradient(Us, axis=1); return Us, Vs, om
d0 = D[0]; _, _, om0 = field(d0, d0["t"][0]); ys = np.arange(d0["Hd"]) + 0.5; VM = np.nanpercentile(np.abs(om0[:, (ys > d0["y0"]) & (ys < d0["y1"])]), 95)
out = os.path.join(H, f"frames_forced_{which}"); os.makedirs(out, exist_ok=True)
for f_ in os.listdir(out): os.remove(os.path.join(out, f_))
n = min(len(d["t"]) for d in D)
for k in range(n):
    fig, axs = plt.subplots(3, 1, figsize=(14, 7.0)); fig.subplots_adjust(left=0.01, right=0.92, top=0.86, bottom=0.02, hspace=0.42)
    for ax, d, lab in zip(axs, D, labels):
        tt = d["t"][k]; f = d["fr"][k]; Us, Vs, om = field(d, tt); xs, ys = np.arange(d["L"]) + 0.5, np.arange(d["Hd"]) + 0.5; y0, y1 = d["y0"], d["y1"]
        om[:, (ys < y0) | (ys > y1)] = np.nan
        im = ax.pcolormesh(xs, ys, om.T, cmap="PuOr_r", vmin=-VM, vmax=VM, shading="nearest", alpha=0.85)
        ax.quiver(xs[2::6], ys[2::4], Us[2::6, 2::4].T, Vs[2::6, 2::4].T, color="0.25", scale=25, width=0.0012)
        ax.axhline(y0, color="k", lw=2); ax.axhline(y1, color="k", lw=2)
        C = np.nan
        if len(f):
            for x, y, q, a_ in f[f[:, 2] > 0]:
                h = a_ + np.pi; ax.plot([x - 2.5 * np.cos(h), x], [y - 2.5 * np.sin(h), y], color="#c2452d", lw=2.2); ax.plot(x, y, "o", ms=5, color="#c2452d", mec="w", mew=0.6)
            for x, y, q, a_ in f[f[:, 2] < 0]:
                for j in range(3): b_ = a_ + 2 * np.pi * j / 3; ax.plot([x, x + 1.8 * np.cos(b_)], [y, y + 1.8 * np.sin(b_)], color="#1f5fa8", lw=1.8)
                ax.plot(x, y, "o", ms=4, color="#1f5fa8", mec="w", mew=0.6)
            bot = f[:, 1] - y0 < 4; top = y1 - f[:, 1] < 4; bulk = ~bot & ~top; w = bot & (f[:, 2] > 0)
            if w.sum() >= 3: C = abs(np.mean(np.exp(1j * (f[w, 3] + np.pi / 2))))
            ax.set_title(f"{lab}:  bottom {int(w.sum())} +1/2 (of {2*abs(d['cb'])}),  top {int(((f[:,2]>0)&top).sum())} +1/2 / {int(((f[:,2]<0)&top).sum())} −1/2,"
                         f"  bulk charge {0.5*np.sign(f[bulk,2]).sum():+.1f},  bottom-layer tilt coherence C = {C:.2f}", fontsize=9.5, loc="left")
        ax.set_xlim(0, d["L"]); ax.set_ylim(y0 - 1, y1 + 1); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    fig.suptitle(f"{'Symmetric (+6 / +6)' if which == 'sym' else 'Mixed (+6 / −6)'} winding channels, W = 24, L = 192,  t = {D[0]['t'][k]:.0f} t.u.\n"
                 f"Colour: flow vorticity (same scale in all panels); arrows: flow; red comets: +1/2; blue tripods: −1/2. Forced flow runs left to right.", fontsize=10.5, x=0.01, ha="left", y=0.985)
    cax = fig.add_axes([0.935, 0.25, 0.012, 0.5]); cb = fig.colorbar(im, cax=cax); cb.set_label("Vorticity (1/t.u.)", fontsize=9); cb.ax.tick_params(labelsize=8)
    fig.savefig(os.path.join(out, f"f{k:04d}.png"), dpi=100); plt.close(fig)
mp4 = os.path.join(H, f"movie_forced_vs_unforced_{which}.mp4")
subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", "15", "-i", os.path.join(out, "f%04d.png"), "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2", "-pix_fmt", "yuv420p", mp4], check=True)
print(mp4, n, "frames")
