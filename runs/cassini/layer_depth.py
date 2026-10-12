#!/usr/bin/env python3
"""
Layer depth from the accumulated charge Q(c) (user's idea, 2026-10-01).

c = sqrt(d1 d2) (circle: c = r). c* = first upward crossing of Q(c) = 1 on the way to Q = k (scanning inwards from
the wall: the outermost c where Q drops below 1). d* = arclength-weighted mean wall distance along the level
curve c = c* (circle: R - c*). Plateau width w = (mean wall distance at the inner end) - (at the outer end) of the
outermost interval around c* where |Q - 1| < --tol. Reported in cells, / l_a, and / size (R or half-length).

Usage: layer_depth.py RUN [RUN ...] [--labels ...] [--la 22.21] [--tol 0.15]
"""
import argparse, json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from lobe_profiles import profile, geometry, cval  # noqa: E402


def wall_curve(a, b, n=4000):
    th = np.linspace(0, 2 * np.pi, n, endpoint=False)
    r2 = a * a * np.cos(2 * th) + np.sqrt(b ** 4 - a ** 4 * np.sin(2 * th) ** 2)
    return np.column_stack([np.sqrt(r2) * np.cos(th), np.sqrt(r2) * np.sin(th)])


def level_mean_distance(a, b, c, wall, h=0.15):
    """mean wall distance over the level set c (points with |c(x) - c| < h/2 on a fine grid ~ arclength weight)."""
    X = np.arange(-np.hypot(a, b) - 1, np.hypot(a, b) + 1, h)
    xx, yy = np.meshgrid(X, X, indexing="ij")
    cc = cval(a, xx, yy)
    m = np.abs(cc - c) < 0.5 * h * max(1.0, np.median(np.abs(np.gradient(cc)[0][np.abs(cc - c) < 2 * h])) / h)
    if not m.any():
        return np.nan
    P = np.column_stack([xx[m], yy[m]])
    d = np.empty(len(P))
    for i in range(0, len(P), 2000):
        blk = P[i:i + 2000]
        d[i:i + 2000] = np.sqrt(((blk[:, None, :] - wall[None, :, :]) ** 2).sum(-1)).min(1)
    return d.mean()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+"); ap.add_argument("--labels", nargs="+")
    ap.add_argument("--la", type=float, default=22.21); ap.add_argument("--tol", type=float, default=0.15)
    a = ap.parse_args()
    labels = a.labels or [os.path.basename(os.path.normpath(r)) for r in a.runs]
    print(f"{'run':<34} {'size':>5} {'c*/b':>5} {'d* cells':>8} {'d*/la':>6} {'d*/size':>7} {'plateau w':>9} {'w/la':>5}")
    rows = []
    for run, lab in zip(a.runs, labels):
        A, B = geometry(run)
        p = profile(run, 200, 0.2)
        x = p["edges"][1:]; Q = p["Q"]; k = p["k"]
        # scan inwards from the wall: outermost index where Q < 1
        idx = np.where(Q < 1.0)[0]
        if not len(idx):
            print(f"{lab:<34} Q never below 1"); continue
        i = idx[-1]
        if i + 1 >= len(Q):
            continue
        xs = x[i] + (1 - Q[i]) / (Q[i + 1] - Q[i]) * (x[i + 1] - x[i])        # linear interpolation of the crossing
        cstar = xs * B
        W = wall_curve(A, B)
        size = np.hypot(A, B)
        dstar = level_mean_distance(A, B, cstar, W)
        # plateau: contiguous band around the crossing with |Q-1| < tol
        ok = np.abs(Q - 1) < a.tol
        lo = hi = i
        while lo > 0 and ok[lo - 1]: lo -= 1
        while hi + 1 < len(Q) and ok[hi + 1]: hi += 1
        if ok[i] or ok[min(i + 1, len(Q) - 1)]:
            d_in = level_mean_distance(A, B, x[lo] * B, W); d_out = level_mean_distance(A, B, x[hi] * B, W)
            w = d_in - d_out
        else:
            w = 0.0
        rows.append((lab, size, xs, dstar, w))
        print(f"{lab:<34} {size:5.1f} {xs:5.3f} {dstar:8.2f} {dstar / a.la:6.3f} {dstar / size:7.3f} {w:9.2f} {w / a.la:5.2f}", flush=True)


if __name__ == "__main__":
    main()
