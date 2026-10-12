#!/usr/bin/env python3
"""
Bulk vs edge charge in Cassini (or circle) cavities, shape-independent.

Per run (first --tfrac of defect frames dropped):
  site rule    pinned sites of the time-averaged +1/2 map and capture radius as in runs/cavity/edge_sites.py
               (classify(); nothing circle-specific); edge = +1/2 within the capture radius of a site.
               Reported: number of sites, mean edge charge Q_edge, bulk charge q_b = k - Q_edge (all
               non-edge defects), P(q_b = +1), and the bulk charge in the left / right half (x < Q_x, x > Q_x).
  Q(s)         enclosed charge: net charge of defects farther than s from the wall, with s the first-order
               distance to the Cassini wall, -W/|grad W| (W = d1 d2 - b^2). Plateau value = bulk charge
               without any cut.
Usage: cassini_charge.py RUN [RUN ...]
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "cavity"))
from defect_spatial import load  # noqa: E402
from edge_sites import classify  # noqa: E402


def wall_distance(bc, x, y):
    a = bc.get("cassiniA", 0.0) if bc.get("shape", 0) == 1 else 0.0
    b = bc["R"]
    X, Y = x, y                       # positions relative to the centre (load() subtracts Q)
    d1, d2 = np.hypot(X - a, Y), np.hypot(X + a, Y)
    W = d1 * d2 - b * b
    gx = d2 / d1 * (X - a) + d1 / d2 * (X + a)
    gy = (d2 / d1 + d1 / d2) * Y
    return -W / np.hypot(gx, gy)


def main():
    runs = sys.argv[1:]
    sgrid = [1, 2, 3, 4, 5, 6, 8, 10]
    print(f"{'run':<20} {'b/a':>5} {'k':>2} {'sites':>5} {'Q_edge':>6} {'q_b':>5} {'P(+1)':>5}  {'q_b left':>8} {'right':>6}  "
          + "  Q(s>" + ") Q(s>".join(str(s) for s in sgrid) + ")")
    for run in runs:
        bc = json.load(open(os.path.join(run, "input.json")))["BC"][0]
        d = load(run, 0.2)
        k = d["k"]
        c = classify(d, 1.0, 0.02, 0.0, 4.2, 0.4)
        sites = c["sites"]
        qL, qR, Qs = [], [], []
        for f in d["frames"]:
            p = f[:, 2] > 0
            edge = np.zeros(len(f), bool)
            if len(sites) and p.any():
                dd = np.hypot(*(f[p, None, :2] - sites[None, :, :]).transpose(2, 0, 1)).min(1)
                edge[np.where(p)[0]] = dd < c["dsite"]
            bulk = ~edge
            qL.append(f[bulk & (f[:, 0] < 0), 2].sum()); qR.append(f[bulk & (f[:, 0] > 0), 2].sum())
            s = wall_distance(bc, f[:, 0], f[:, 1])
            Qs.append([f[s > s0, 2].sum() for s0 in sgrid])
        a = bc.get("cassiniA", 0.0) if bc.get("shape", 0) == 1 else 0.0
        ba = bc["R"] / a if a else float("inf")
        Qs = np.mean(Qs, 0)
        print(f"{os.path.basename(os.path.normpath(run)):<20} {ba:5.2f} {k:2d} {c['nsite']:5d} {k - c['qb']:6.2f} {c['qb']:5.2f} {c['P1']:5.2f}  "
              f"{np.mean(qL):8.2f} {np.mean(qR):6.2f}  " + " ".join(f"{q:7.2f}" for q in Qs))


if __name__ == "__main__":
    main()
