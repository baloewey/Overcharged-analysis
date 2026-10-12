#!/usr/bin/env python3
"""
Area occupied by the "layer" (user's idea, 2026-10-01): region outside the level c* where Q(c) crosses +1.
f* = A(c > c*) / A_total;  w_eff = A(c > c*) / P (equivalent uniform strip width; P = wall perimeter), also / l_a.
c* as in layer_depth.py (outermost upward crossing of Q = +1). Usage: layer_area.py RUN [RUN ...] [--labels ...]
"""
import argparse, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from lobe_profiles import profile, geometry, cval  # noqa: E402
from layer_depth import wall_curve  # noqa: E402


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("runs", nargs="+"); ap.add_argument("--labels", nargs="+")
    ap.add_argument("--la", type=float, default=22.21)
    a = ap.parse_args()
    labels = a.labels or [os.path.normpath(r) for r in a.runs]
    print(f"{'run':<34} {'k':>3} {'c*/b':>5} {'f*':>6} {'w_eff':>6} {'w_eff/la':>8}")
    for run, lab in zip(a.runs, labels):
        A, B = geometry(run)
        p = profile(run, 200, 0.2)
        x = p["edges"][1:]; Q = p["Q"]
        idx = np.where(Q < 1.0)[0]; i = idx[-1]
        xs = x[i] + (1 - Q[i]) / (Q[i + 1] - Q[i]) * (x[i + 1] - x[i])
        h = 0.1; X = np.arange(-np.hypot(A, B) - 1, np.hypot(A, B) + 1, h)
        xx, yy = np.meshgrid(X, X, indexing="ij"); c = cval(A, xx, yy)
        tot = np.sum(c < B) * h * h; out = np.sum((c < B) & (c > xs * B)) * h * h
        W = wall_curve(A, B); P = np.sum(np.hypot(*np.diff(np.vstack([W, W[:1]]), axis=0).T))
        print(f"{lab:<34} {p['k']:>3} {xs:5.3f} {out / tot:6.3f} {out / P:6.2f} {out / P / a.la:8.3f}", flush=True)


if __name__ == "__main__":
    main()
