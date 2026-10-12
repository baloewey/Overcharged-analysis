#!/usr/bin/env python3
"""
Cut-free edge/bulk criteria for +1/2 defects, tested against the pinned-site rule.

Per +1/2 defect occurrence (first --tfrac of frames dropped):
  s        distance from the wall (circle, or Cassini via -W/|grad W|)
  D        mobility: net displacement |r(t + tau/2) - r(t - tau/2)| along its linked track (frame-to-frame
           optimal assignment of same-sign defects, max jump --rlink); occurrences whose track does not
           cover the whole window are skipped
  Delta    orientation relative to the local wall normal: angle column of defects.dat minus the angle of the
           outward wall normal at the defect position, wrapped to (-pi, pi]
  site     pinned-site label from edge_sites.classify (reference where it works)
Output: text summary (distributions split by site label, best single thresholds and their agreement with the
site rule) and a figure PREFIX.png per run set.

Usage: edge_criteria.py RUN [RUN ...] [--tau 20] [--labels ...] [--out PREFIX]
"""
import argparse
import gzip
import json
import os
import sys

import numpy as np
from scipy.optimize import linear_sum_assignment

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from defect_spatial import load  # noqa: E402
from edge_sites import classify  # noqa: E402


def opendat(path):
    return open(path) if os.path.exists(path) else gzip.open(path + ".gz", "rt")


def read_with_angle(run, Q, tfrac):
    times, frames, cur, nexp, tcur = [], [], None, 0, None
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
    n0 = int(tfrac * len(frames))
    out = []
    for f in frames[n0:]:
        g = f.copy()
        g[:, 0] -= Q[0]; g[:, 1] -= Q[1]
        out.append(g)
    return np.array(times[n0:]), out


def wall(bc, x, y):
    """-> distance from the wall (>0 inside), outward-normal angle."""
    a = bc.get("cassiniA", 0.0) if bc.get("shape", 0) == 1 else 0.0
    b = bc["R"]
    d1, d2 = np.hypot(x - a, y), np.hypot(x + a, y)
    W = d1 * d2 - b * b
    gx = d2 / np.maximum(d1, 1e-9) * (x - a) + d1 / np.maximum(d2, 1e-9) * (x + a)
    gy = (d2 / np.maximum(d1, 1e-9) + d1 / np.maximum(d2, 1e-9)) * y
    return -W / np.hypot(gx, gy), np.arctan2(gy, gx)


def link(frames, rlink):
    """Track ids for +1/2 defects: list of arrays (one per frame, -1 for -1/2)."""
    ids, nxt = [], 0
    prev = None
    for f in frames:
        lab = -np.ones(len(f), int)
        p = np.where(f[:, 2] > 0)[0]
        if prev is not None and len(p) and len(prev[1]):
            A, ia = prev[0], prev[1]
            D = np.hypot(*(A[:, None, :] - f[p, None, :2].transpose(1, 0, 2)).transpose(2, 0, 1))
            r, c = linear_sum_assignment(np.where(D < rlink, D, 1e6))
            for i, j in zip(r, c):
                if D[i, j] < rlink:
                    lab[p[j]] = ia[i]
        for j in p:
            if lab[j] < 0:
                lab[j] = nxt; nxt += 1
        ids.append(lab)
        prev = (f[p, :2], lab[p])
    return ids


def analyse(run, tau, rlink, tfrac, flip=False):
    j = json.load(open(os.path.join(run, "input.json")))
    bc = j["BC"][0]
    Q = np.array(bc["Q"][:2], float)
    t, fr = read_with_angle(run, Q, tfrac)
    if flip:                                   # mirror runs: follow the wall's -1/2 as if they were +1/2
        for f in fr:
            f[:, 2] *= -1
    dt = t[1] - t[0]
    h = max(1, int(round(0.5 * tau / dt)))
    ids = link(fr, rlink)
    # position lookup per track
    pos = {}
    for n, (f, lab) in enumerate(zip(fr, ids)):
        for i in np.where(lab >= 0)[0]:
            pos[(lab[i], n)] = f[i, :2]
    d = load(run, tfrac)
    c = classify(d, 1.0, 0.02, 0.0, 4.2, 0.4)
    sites = c["sites"]
    S, D, Dl, site = [], [], [], []
    for n, (f, lab) in enumerate(zip(fr, ids)):
        if n - h < 0 or n + h >= len(fr):
            continue
        p = np.where(lab >= 0)[0]
        if not len(p):
            continue
        s, nang = wall(bc, f[p, 0], f[p, 1])
        dl = np.angle(np.exp(1j * (f[p, 3] - nang)))
        if len(sites):
            ds = np.hypot(*(f[p, None, :2] - sites[None, :, :]).transpose(2, 0, 1)).min(1)
            st = ds < c["dsite"]
        else:
            st = np.zeros(len(p), bool)
        for i, si, di, sti in zip(p, s, dl, st):
            a0, a1 = pos.get((lab[i], n - h)), pos.get((lab[i], n + h))
            if a0 is None or a1 is None:
                continue
            S.append(si); D.append(np.hypot(*(a1 - a0))); Dl.append(di); site.append(sti)
    return dict(run=run, k=bc.get("customAnch"), R=bc["R"], S=np.array(S), D=np.array(D), Dl=np.array(Dl),
                site=np.array(site), nsite=c["nsite"], tau=2 * h * dt, nframes=max(len(fr) - 2 * h, 1),
                nplus_frame=float(np.mean([np.sum(f[:, 2] > 0) for f in fr])))


def best_threshold(x, y, lo, hi, n=200):
    """threshold on x (x < thr -> edge) maximising agreement with boolean y."""
    thr = np.linspace(lo, hi, n)
    acc = [np.mean((x < t_) == y) for t_ in thr]
    i = int(np.argmax(acc))
    return thr[i], acc[i]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--labels", nargs="+")
    ap.add_argument("--tau", type=float, default=20.0)
    ap.add_argument("--rlink", type=float, default=3.0)
    ap.add_argument("--tfrac", type=float, default=0.2)
    ap.add_argument("--out", default="edge_criteria")
    a = ap.parse_args()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    labels = a.labels or [os.path.normpath(r) for r in a.runs]
    res = [analyse(r, a.tau, a.rlink, a.tfrac) for r in a.runs]
    print(f"tau = {res[0]['tau']:.0f} time units; edge reference = pinned-site rule")
    print(f"{'run':<26} {'k':>3} {'sites':>5} {'n':>7} {'D edge':>7} {'D bulk':>7} {'D thr':>6} {'agree':>6} | "
          f"{'|Dl| edge':>9} {'|Dl| bulk':>9} {'Dl thr':>6} {'agree':>6} | {'both agree':>10}")
    for r, lab in zip(res, labels):
        e, b = r["site"], ~r["site"]
        if e.sum() and b.sum():
            dthr, dacc = best_threshold(r["D"], e, 0.5, 30)
            athr, aacc = best_threshold(np.abs(r["Dl"]), e, 0.05, np.pi)
            both = np.mean(((r["D"] < dthr) & (np.abs(r["Dl"]) < athr)) == e)
            print(f"{lab:<26} {r['k']:>3} {r['nsite']:>5} {len(r['D']):>7} {np.median(r['D'][e]):7.2f} {np.median(r['D'][b]):7.2f} "
                  f"{dthr:6.2f} {dacc:6.3f} | {np.median(np.abs(r['Dl'][e])):9.2f} {np.median(np.abs(r['Dl'][b])):9.2f} {athr:6.2f} {aacc:6.3f} | {both:10.3f}")
        else:
            print(f"{lab:<26} {r['k']:>3} {r['nsite']:>5} {len(r['D']):>7}  (no site reference)")
    n = len(res)
    fig, axs = plt.subplots(3, n, figsize=(3.2 * n, 8.4), squeeze=False)
    for col, (r, lab) in enumerate(zip(res, labels)):
        e = r["site"]
        ax = axs[0, col]
        bins = np.linspace(0, 40, 81)
        ax.hist(r["D"][e], bins, color="#2a78d6", alpha=0.75, density=True, label="site rule: edge")
        ax.hist(r["D"][~e], bins, color="#eb6834", alpha=0.6, density=True, label="site rule: bulk")
        ax.set_yscale("log"); ax.set_xlabel(f"displacement over {r['tau']:.0f} time units"); ax.set_title(lab, fontsize=9)
        ax = axs[1, col]
        ax.scatter(r["S"], r["D"], s=1, alpha=0.05, color="#2b2f36", lw=0)
        ax.set_xlim(0, 15); ax.set_ylim(0, 40); ax.set_xlabel("distance from wall s"); ax.set_ylabel("displacement")
        ax = axs[2, col]
        ax.scatter(r["S"], r["Dl"], s=1, alpha=0.05, color="#2b2f36", lw=0)
        ax.set_xlim(0, 15); ax.set_ylim(-np.pi, np.pi); ax.set_xlabel("distance from wall s"); ax.set_ylabel("orientation − wall normal")
    axs[0, 0].legend(fontsize=7, frameon=False)
    fig.tight_layout(); fig.savefig(a.out + ".png", dpi=140)
    np.savez_compressed(a.out + ".npz", **{f"{i}_{k}": v for i, r in enumerate(res) for k, v in r.items() if isinstance(v, np.ndarray)})
    print("written", a.out + ".png")


if __name__ == "__main__":
    main()
