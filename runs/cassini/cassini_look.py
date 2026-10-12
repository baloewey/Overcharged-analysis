#!/usr/bin/env python3
"""
Quick look at Cassini (or circle) cavity runs: wall alignment, charge audit, accessible volume, figure.

Per run directory (input.json, synopsis.dat, defects.dat[.gz], directorfield.dat[.gz]):
  alignment   <cos 2(phi_dir - phi_imposed)> over cells whose centre lies within 1.5 cells inside the
              wall, all director frames after the first; phi_imposed = k * 2 pi s/P (eMode 3: s from the
              polar angle; eMode 0: k * normal angle). 1 = perfect, 0 = random.
  audit       fraction of frames with net charge = k and with the 2|k| imposed defects present
              (defect_audit.py), number of +1/2 and -1/2 per frame
  volume      accessible volume from synopsis.dat vs the numerical oval area
  figure      oval outline, last director frame, time-averaged defect positions (+1/2 blue, -1/2 orange)

Usage: cassini_look.py RUN [RUN ...] [--out PREFIX]
"""
import argparse
import gzip
import json
import math
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "cavity"))
from defect_audit import audit_frames, read_frames  # noqa: E402


def opendat(path):
    if os.path.exists(path):
        return open(path)
    return gzip.open(path + ".gz", "rt")


def read_director(run):
    """-> list of (n,5) arrays [x, y, nx, ny, S] per frame (cell centres)."""
    frames, cur, t0 = [], [], None
    with opendat(os.path.join(run, "directorfield.dat")) as f:
        for line in f:
            p = line.split()
            if len(p) != 8:
                continue
            try:
                t, x, y, z, nx, ny, nz, S = map(float, p)
            except ValueError:
                continue
            if t0 is None:
                t0 = t
            if t != t0:
                frames.append(np.array(cur)); cur = []; t0 = t
            cur.append((x + 0.5, y + 0.5, nx, ny, S))
    if cur:
        frames.append(np.array(cur))
    return frames


def geom(run):
    j = json.load(open(os.path.join(run, "input.json")))
    bc = j["BC"][0]
    return dict(Q=np.array(bc["Q"][:2], float), b=bc["R"], a=bc.get("cassiniA", 0.0) if bc.get("shape", 0) == 1 else 0.0,
                k=bc.get("customAnch", 0), emode=bc.get("eMode", 0), cassini=bc.get("shape", 0) == 1,
                act=j["species"][0].get("act", 0.0), j=j)


def W(g, x, y):
    X, Y = x - g["Q"][0], y - g["Q"][1]
    a = g["a"]
    return np.hypot(X - a, Y) * np.hypot(X + a, Y) - g["b"] ** 2


def grad(g, x, y):
    X, Y = x - g["Q"][0], y - g["Q"][1]
    a = g["a"]
    d1, d2 = np.hypot(X - a, Y), np.hypot(X + a, Y)
    return d2 / d1 * (X - a) + d1 / d2 * (X + a), (d2 / d1 + d1 / d2) * Y


def arclength_fraction(g, theta):
    a, b = g["a"], g["b"]
    th = np.linspace(0, 2 * np.pi, 20001)
    r2 = a * a * np.cos(2 * th) + np.sqrt(b ** 4 - a ** 4 * np.sin(2 * th) ** 2)
    x, y = np.sqrt(r2) * np.cos(th), np.sqrt(r2) * np.sin(th)
    s = np.concatenate([[0], np.cumsum(np.hypot(np.diff(x), np.diff(y)))])
    return np.interp(theta, th, s / s[-1])


def imposed_angle(g, x, y):
    if g["emode"] == 3:
        th = np.mod(np.arctan2(y - g["Q"][1], x - g["Q"][0]), 2 * np.pi)
        return g["k"] * 2 * np.pi * arclength_fraction(g, th)
    gx, gy = grad(g, x, y)
    return g["k"] * np.mod(np.arctan2(gy, gx), 2 * np.pi)


def wall_distance(g, x, y):
    gx, gy = grad(g, x, y)
    return -W(g, x, y) / np.hypot(gx, gy)       # first-order distance inside the wall (>0 inside)


def analyse(run):
    g = geom(run)
    out = dict(run=run, k=g["k"], emode=g["emode"], a=g["a"], b=g["b"], act=g["act"])
    D = read_director(run)
    al = []
    for f in D[1:]:
        x, y, nx, ny, S = f.T
        dist = wall_distance(g, x, y)
        m = (dist > 0) & (dist < 1.5) & (S > 0)
        if m.any():
            psi = imposed_angle(g, x[m], y[m])
            al.append(np.mean(np.cos(2 * (np.arctan2(ny[m], nx[m]) - psi))))
    out["align"] = float(np.mean(al)) if al else np.nan
    t, fr = read_frames(os.path.join(run, "defects.dat"))
    fr = fr[int(0.2 * len(fr)):]
    a_ = audit_frames(fr, g["k"])
    out["frac_net_k"] = float(np.mean(a_["ok_net"]))
    out["frac_imposed"] = float(np.mean(a_["ok_imposed"]))
    out["n_plus"] = float(np.mean([np.sum(f[:, 2] > 0) for f in fr]))
    out["n_minus"] = float(np.mean([np.sum(f[:, 2] < 0) for f in fr]))
    syn = open(os.path.join(run, "synopsis.dat")).read()
    m = re.search(r"Accessible volume:\s*([0-9.]+)", syn)
    out["vol"] = float(m.group(1)) if m else np.nan
    th = np.linspace(0, 2 * np.pi, 20000, endpoint=False)
    r2 = g["a"] ** 2 * np.cos(2 * th) + np.sqrt(g["b"] ** 4 - g["a"] ** 4 * np.sin(2 * th) ** 2)
    out["area"] = float(0.5 * r2.mean() * 2 * np.pi)
    out["_g"], out["_D"], out["_fr"] = g, D, fr
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--out", default="cassini_look")
    a = ap.parse_args()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    R = [analyse(r) for r in a.runs]
    print(f"{'run':<22} {'k':>3} {'eMode':>5} {'b/a':>6} {'align':>6} {'net=k':>6} {'imposed':>7} {'n+':>6} {'n-':>5} {'vol':>8} {'area':>8}")
    for r in R:
        ba = r["b"] / r["a"] if r["a"] else float("inf")
        print(f"{r['run']:<22} {r['k']:>3} {r['emode']:>5} {ba:6.2f} {r['align']:6.3f} {r['frac_net_k']:6.3f} {r['frac_imposed']:7.3f} "
              f"{r['n_plus']:6.2f} {r['n_minus']:5.2f} {r['vol']:8.1f} {r['area']:8.1f}")
    n = len(R)
    fig, axs = plt.subplots(1, n, figsize=(3.6 * n, 3.4), squeeze=False)
    for ax, r in zip(axs[0], R):
        g = r["g"] if "g" in r else r["_g"]
        th = np.linspace(0, 2 * np.pi, 800)
        r2 = g["a"] ** 2 * np.cos(2 * th) + np.sqrt(g["b"] ** 4 - g["a"] ** 4 * np.sin(2 * th) ** 2)
        ax.plot(g["Q"][0] + np.sqrt(r2) * np.cos(th), g["Q"][1] + np.sqrt(r2) * np.sin(th), color="#555a62", lw=1)
        f = r["_D"][-1]
        x, y, nx, ny, S = f.T
        m = S > 0.1
        h = 0.4
        ax.plot(np.vstack([x[m] - h * nx[m], x[m] + h * nx[m]]), np.vstack([y[m] - h * ny[m], y[m] + h * ny[m]]),
                color="#9aa0a8", lw=0.6)
        allf = np.vstack([q for q in r["_fr"] if len(q)])
        for sg, col in ((1, "#2a78d6"), (-1, "#eb6834")):
            s = allf[np.sign(allf[:, 2]) == sg]
            ax.scatter(s[:, 0], s[:, 1], s=1, color=col, alpha=0.05, lw=0)
        ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        ba = g["b"] / g["a"] if g["a"] else float("inf")
        ax.set_title(f"{os.path.basename(os.path.normpath(r['run']))}\nk={g['k']}, b/a={ba:.2f}, eMode {g['emode']}, align {r['align']:.2f}", fontsize=8)
    fig.tight_layout(); fig.savefig(a.out + ".png", dpi=150)
    print("written", a.out + ".png")


if __name__ == "__main__":
    main()
