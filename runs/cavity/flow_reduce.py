#!/usr/bin/env python3
"""
Per-frame flow reduction of circle runs (runs on the server, where the flow fields are).

flowfield.dat columns: t QX QY QZ VX VY VZ (all cells of the box). Cell centre = (QX + 0.5, QY + 0.5) - BC Q.
Cells with r < R - 1 are used (cells cut by the wall dropped). Per frame:
  vth   spatial mean of the azimuthal velocity v_theta = (x vy - y vx)/r      (sign = circulation sense, + = ccw)
  vrms  sqrt of the spatial mean of |v|^2
  ens   spatial mean of omega^2 / 2, omega = d vy/dx - d vx/dy (centred differences, both neighbours inside)
  wm    spatial mean of omega
Writes OUT/<run path with / -> __>.npz with t, vth, vrms, ens, wm.
Usage: flow_reduce.py RUNLIST ROOT OUT [--workers 40]
"""
import argparse, gzip, json, os, sys
from multiprocessing import Pool
import numpy as np


def one(args):
    root, run, out = args
    d = os.path.join(root, run)
    try:
        bc = json.load(open(os.path.join(d, "input.json")))["BC"][0]
        R = bc["R"]; cx, cy = bc["Q"][0], bc["Q"][1]
        with gzip.open(os.path.join(d, "flowfield.dat.gz"), "rt") as f:
            rows = []
            for line in f:
                p = line.split()
                if len(p) == 7 and p[0][0] in "0123456789":
                    rows.append(p)
        A = np.array(rows, dtype=float)
        ts = np.unique(A[:, 0]); nx = int(A[:, 1].max()) + 1; ny = int(A[:, 2].max()) + 1
        A = A[np.lexsort((A[:, 2], A[:, 1], A[:, 0]))]
        V = A[:, 4:6].reshape(len(ts), nx, ny, 2)
        X, Y = np.meshgrid(np.arange(nx) + 0.5 - cx, np.arange(ny) + 0.5 - cy, indexing="ij")
        r = np.hypot(X, Y); m = r < R - 1
        mi = m.copy(); mi[1:-1, 1:-1] &= m[2:, 1:-1] & m[:-2, 1:-1] & m[1:-1, 2:] & m[1:-1, :-2]
        mi[0, :] = mi[-1, :] = mi[:, 0] = mi[:, -1] = False
        vx, vy = V[..., 0], V[..., 1]
        vth = ((X * vy - Y * vx) / np.maximum(r, 1e-9))[:, m].mean(1)
        vrms = np.sqrt((vx ** 2 + vy ** 2)[:, m].mean(1))
        w = np.zeros_like(vx)
        w[:, 1:-1, 1:-1] = 0.5 * (vy[:, 2:, 1:-1] - vy[:, :-2, 1:-1]) - 0.5 * (vx[:, 1:-1, 2:] - vx[:, 1:-1, :-2])
        ens = 0.5 * (w ** 2)[:, mi].mean(1); wm = w[:, mi].mean(1)
        np.savez(os.path.join(out, run.replace("/", "__") + ".npz"), t=ts, vth=vth, vrms=vrms, ens=ens, wm=wm, R=R, ncell=m.sum())
        return None
    except Exception as ex:
        return f"{run}: {ex}"


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("runlist"); ap.add_argument("root"); ap.add_argument("out")
    ap.add_argument("--workers", type=int, default=40); a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    runs = [l.strip() for l in open(a.runlist) if l.strip()]
    with Pool(a.workers) as pool:
        for i, e in enumerate(pool.imap_unordered(one, [(a.root, r, a.out) for r in runs])):
            if e: print("FAILED", e, flush=True)
            if i % 100 == 0: print(i, flush=True)
    print("done")


if __name__ == "__main__":
    main()
