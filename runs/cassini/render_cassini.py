#!/usr/bin/env python3
"""
Render Cassini (or circle) cavity runs: snapshots and an MP4 movie.

Per frame (one per director/flow output):
  background  vorticity w = d_x v_y - d_y v_x of the cell flow field, smoothed over 1 cell; blue = clockwise
              (w < 0), orange = counter-clockwise (w > 0); colour scale fixed per run (98th percentile of |w|)
  lines       director field (one segment per cell, length ~ S)
  defects     +1/2: black dot with a tail along its orientation; -1/2: white triangle (nearest defect frame)
  title       time, and the mean azimuthal velocity about the centroid of the left (x < Q_x) and right
              (x > Q_x) halves (the mean vorticity of a no-slip region is ~0 by Stokes, so it can't give the sense)
Works on finished (gzipped) and running runs (a truncated last frame is dropped).

Usage:
  render_cassini.py RUN [RUN ...] [--frames all|last|N] [--outdir DIR] [--fps 8]
Writes DIR/<run>_last.png and, unless --frames last, DIR/<run>.mp4.
"""
import argparse
import gzip
import json
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np
from scipy.ndimage import gaussian_filter

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "cavity"))
from defect_audit import read_frames  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402

VORT = LinearSegmentedColormap.from_list("vort", ["#184f95", "#6da7ec", "#f4f5f7", "#f4a27f", "#b84a1d"])


def opendat(path):
    return open(path) if os.path.exists(path) else gzip.open(path + ".gz", "rt")


def read_grid_frames(path, ncol, Lx, Ly):
    """-> times, array (nframes, Lx, Ly, ncol-4) of the per-cell values."""
    times, frames, cur, t0 = [], [], {}, None
    with opendat(path) as f:
        for line in f:
            p = line.split()
            if len(p) != ncol:
                continue
            try:
                v = [float(x) for x in p]
            except ValueError:
                continue
            if t0 is None:
                t0 = v[0]
            if v[0] != t0:
                frames.append(cur); times.append(t0); cur = {}; t0 = v[0]
            cur[(int(v[1]), int(v[2]))] = v[4:]
    if cur and len(cur) == Lx * Ly:
        frames.append(cur); times.append(t0)
    out = np.zeros((len(frames), Lx, Ly, ncol - 4))
    for n, fr in enumerate(frames):
        for (i, j), val in fr.items():
            if i < Lx and j < Ly:
                out[n, i, j] = val
    return np.array(times), out


def read_defects_angle(run):
    """-> times, list of (n,4) arrays [x, y, charge, angle] (defects.dat has the +1/2 orientation)."""
    times, frames, cur, nexp = [], [], None, 0
    with opendat(os.path.join(run, "defects.dat")) as f:
        for line in f:
            p = line.split()
            if len(p) == 2:
                try:
                    t, nn = float(p[0]), int(p[1])
                except ValueError:
                    continue
                if cur is not None and len(cur) == nexp:
                    frames.append(np.array(cur).reshape(-1, 4)); times.append(tcur)
                cur, nexp, tcur = [], nn, t
            elif len(p) == 4 and cur is not None:
                try:
                    cur.append([float(x) for x in p])
                except ValueError:
                    pass
    if cur is not None and len(cur) == nexp:
        frames.append(np.array(cur).reshape(-1, 4)); times.append(tcur)
    return np.array(times), frames


def geom(run):
    j = json.load(open(os.path.join(run, "input.json")))
    bc = j["BC"][0]
    a = bc.get("cassiniA", 0.0) if bc.get("shape", 0) == 1 else 0.0
    return dict(Lx=j["domain"][0], Ly=j["domain"][1], Q=np.array(bc["Q"][:2], float), b=bc["R"], a=a,
                k=bc.get("customAnch", 0), act=j["species"][0].get("act", 0.0),
                dt_def=j.get("defectsOut", 20) * 0.1, warm=j.get("warmUp", 0) * 0.1)


def inside(g, x, y):
    X, Y = x - g["Q"][0], y - g["Q"][1]
    return np.hypot(X - g["a"], Y) * np.hypot(X + g["a"], Y) < g["b"] ** 2


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--frames", default="all")
    ap.add_argument("--outdir", default="movies")
    ap.add_argument("--fps", type=int, default=8)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    for run in a.runs:
        g = geom(run)
        Lx, Ly = g["Lx"], g["Ly"]
        td, D = read_grid_frames(os.path.join(run, "directorfield.dat"), 8, Lx, Ly)
        tf, F = read_grid_frames(os.path.join(run, "flowfield.dat"), 7, Lx, Ly)
        n = min(len(td), len(tf))
        if n == 0:
            print(run, ": no complete frames yet"); continue
        tdef, frd = read_defects_angle(run)
        xc, yc = np.meshgrid(np.arange(Lx) + 0.5, np.arange(Ly) + 0.5, indexing="ij")
        mask = inside(g, xc, yc)
        W = []
        for m in range(n):
            vx, vy = F[m, :, :, 0], F[m, :, :, 1]
            w = np.gradient(vy, axis=0) - np.gradient(vx, axis=1)
            W.append(np.where(mask, gaussian_filter(np.where(mask, w, 0.0), 1.0), np.nan))
        W = np.array(W)
        wmax = np.nanpercentile(np.abs(W[n // 5:]), 98) if n > 5 else np.nanmax(np.abs(W))
        th = np.linspace(0, 2 * np.pi, 800)
        r2 = g["a"] ** 2 * np.cos(2 * th) + np.sqrt(g["b"] ** 4 - g["a"] ** 4 * np.sin(2 * th) ** 2)
        ox, oy = g["Q"][0] + np.sqrt(r2) * np.cos(th), g["Q"][1] + np.sqrt(r2) * np.sin(th)
        left = mask & (xc < g["Q"][0]); right = mask & (xc > g["Q"][0])
        cl = np.array([xc[left].mean(), yc[left].mean()]); cr = np.array([xc[right].mean(), yc[right].mean()])

        def vphi(m, sel, c):
            X, Y = xc[sel] - c[0], yc[sel] - c[1]
            vx, vy = F[m, :, :, 0][sel], F[m, :, :, 1][sel]
            r = np.hypot(X, Y)
            ok = r > 1
            return np.mean((X[ok] * vy[ok] - Y[ok] * vx[ok]) / r[ok])
        name = os.path.basename(os.path.normpath(run))
        sel = range(n) if a.frames == "all" else ([n - 1] if a.frames == "last" else range(0, n, max(1, n // int(a.frames))))
        tmp = tempfile.mkdtemp()
        w_in = 9.0
        h_in = w_in * Ly / Lx + 1.1
        for c, m in enumerate(sel):
            fig, ax = plt.subplots(figsize=(w_in, h_in))
            ax.imshow(W[m].T, origin="lower", extent=[0, Lx, 0, Ly], cmap=VORT, vmin=-wmax, vmax=wmax, interpolation="bilinear")
            nx, ny, S = D[m, :, :, 0], D[m, :, :, 1], D[m, :, :, 3]
            mm = mask & (S > 0.05)
            hl = 0.45 * np.clip(S[mm], 0.3, 1.0)
            segs = np.stack([np.stack([xc[mm] - hl * nx[mm], yc[mm] - hl * ny[mm]], -1),
                             np.stack([xc[mm] + hl * nx[mm], yc[mm] + hl * ny[mm]], -1)], 1)
            ax.add_collection(LineCollection(segs, colors="#2b2f36", linewidths=0.55, alpha=0.75))
            ax.plot(ox, oy, color="#2b2f36", lw=1.4)
            if len(tdef):
                i = int(np.argmin(np.abs(tdef - td[m])))
                dfr = frd[i]
                if len(dfr):
                    P = dfr[dfr[:, 2] > 0]; M = dfr[dfr[:, 2] < 0]
                    ax.scatter(P[:, 0] + 0.5, P[:, 1] + 0.5, s=26, color="black", zorder=5, lw=0)
                    for (x0, y0, ag) in zip(P[:, 0] + 0.5, P[:, 1] + 0.5, P[:, 3]):
                        ax.plot([x0, x0 - 2.2 * np.cos(ag)], [y0, y0 - 2.2 * np.sin(ag)], color="black", lw=1.6, zorder=5)
                    ax.scatter(M[:, 0] + 0.5, M[:, 1] + 0.5, s=34, marker="^", facecolor="white", edgecolor="black", lw=0.9, zorder=6)
            wl, wr = vphi(m, left, cl), vphi(m, right, cr)
            ba = g["b"] / g["a"] if g["a"] else float("inf")
            ax.set_title(f"{name}   k = {g['k']}, b/a = {ba:.2f}, ζ = {g['act']}   t = {td[m]:.0f}\n"
                         f"mean azimuthal velocity about each lobe centre: left {wl:+.4f}, right {wr:+.4f} (+ = counter-clockwise)\n"
                         f"background: vorticity, orange = counter-clockwise, blue = clockwise", fontsize=9)
            ax.set_xlim(0, Lx); ax.set_ylim(0, Ly); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_visible(False)
            fig.tight_layout()
            fig.savefig(os.path.join(tmp, f"f{c:04d}.png"), dpi=130)
            if m == n - 1:
                fig.savefig(os.path.join(a.outdir, f"{name}_last.png"), dpi=150)
            plt.close(fig)
        if a.frames != "last" and shutil.which("ffmpeg"):
            out = os.path.join(a.outdir, f"{name}.mp4")
            subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(a.fps), "-i", os.path.join(tmp, "f%04d.png"),
                            "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2", "-c:v", "libx264", "-pix_fmt", "yuv420p", out], check=True)
            print("written", out)
        shutil.rmtree(tmp)
        print("written", os.path.join(a.outdir, f"{name}_last.png"), f"({n} frames)")


if __name__ == "__main__":
    main()
