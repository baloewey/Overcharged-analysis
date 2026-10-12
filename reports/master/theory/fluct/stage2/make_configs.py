#!/usr/bin/env python3
"""Stage 2: imperfect wall lattices -> configurations for the Stokes feedback matrix (stage2_modes.wl).
Each configuration: ellipse semi-axes A, B and sites [x, y, tx, ty] (position of a wall +1/2 and the ccw wall tangent at the nearest wall point).
Synthetic (area 2500, N = 10 evenly spaced in arclength at depth 1.5, as in stokes_tilt_modes.wl):
  'disorder' : Gaussian tangential displacement sigma_t in {1, 2, 4} cells and normal sigma_n = 0.5 (depth clipped to [0.8, 4]);
               sigma_t = 1 is the passive level (stage 1), 2 and 4 stand for active wandering;
  'vacancy'  : 1 or 2 sites removed at random, plus sigma_t = 1;
  20 samples each, eps in {1.0001, 1.2, 1.5, 2.0}.
Measured: W3 active ellipses (seed 1), one defect snapshot at the middle of each of the 40 flow intervals used in stokes_frames
  (wall = +1/2 within 4 cells of the wall), so that the predicted top mode can be compared with the measured flow of that interval.
Writes cfg_<tag>.json, one file per eps (run in parallel)."""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = "/home/baloewe/Physics/Overcharged/runs"
sys.path.insert(0, os.path.join(RUNS, "cavity"))
from edge_criteria import read_with_angle  # noqa
rng = np.random.default_rng(20261006)
NS, NN, DEPTH = 20, 10, 1.5


def wall(A, B, n=8000):
    th = np.linspace(0, 2 * np.pi, n, endpoint=False)
    W = np.column_stack([A * np.cos(th), B * np.sin(th)])
    seg = np.hypot(*np.diff(np.vstack([W, W[:1]]), axis=0).T)
    s = np.concatenate([[0], np.cumsum(seg)[:-1]])
    T = np.column_stack([-A * np.sin(th), B * np.cos(th)]); T /= np.linalg.norm(T, axis=1)[:, None]
    Nn = np.column_stack([T[:, 1], -T[:, 0]])                      # outward normal
    return W, s, seg.sum(), T, Nn


def at(Wd, sa, depth):
    W, s, P, T, Nn = Wd
    i = np.searchsorted(s, sa % P) % len(s)
    p = W[i] - depth * Nn[i]
    return [float(p[0]), float(p[1]), float(T[i, 0]), float(T[i, 1])]


def nearest(Wd, xy):
    W, s, P, T, Nn = Wd
    i = ((W - xy) ** 2).sum(1).argmin()
    return i, np.hypot(*(W[i] - xy))


out = {}
for eps, tag, rel in ((1.0001, "e10", None), (1.2, "e12", "cassini/W3_ellipse/eps1.2_k06/s01"),
                      (1.5, "e15", "cassini/W3_ellipse/eps1.5_k06/s01"), (2.0, "e20", "cassini/W3_ellipse/eps2.0_k06/s01")):
    A = np.sqrt(2500 * eps / np.pi); B = A / eps
    Wd = wall(A, B); P = Wd[2]
    base = (np.arange(NN) + 0.5) * P / NN
    cfgs = [dict(id="perfect", kind="perfect", A=A, B=B, sites=[at(Wd, sb, DEPTH) for sb in base])]
    for st in (1.0, 2.0, 4.0):
        for k in range(NS):
            sa = base + rng.normal(0, st, NN); dp = np.clip(DEPTH + rng.normal(0, 0.5, NN), 0.8, 4)
            cfgs.append(dict(id=f"dis{st:g}_{k}", kind="disorder", sigma=st, A=A, B=B, sites=[at(Wd, a, d) for a, d in zip(sa, dp)]))
    for nv in (1, 2):
        for k in range(NS):
            keep = np.sort(rng.choice(NN, NN - nv, replace=False))
            sa = base[keep] + rng.normal(0, 1.0, NN - nv); dp = np.clip(DEPTH + rng.normal(0, 0.5, NN - nv), 0.8, 4)
            cfgs.append(dict(id=f"vac{nv}_{k}", kind="vacancy", vac=nv, A=A, B=B, sites=[at(Wd, a, d) for a, d in zip(sa, dp)]))
    if rel:
        run = os.path.join(RUNS, rel); bc = json.load(open(run + "/input.json"))["BC"][0]
        fr_json = json.load(open(os.path.join(HERE, "..", "..", "stokes_frames", tag + ".json")))
        Am, Bm = fr_json["A"], fr_json["B"]; Wm = wall(Am, Bm)
        t, fr = read_with_angle(run, np.array(bc["Q"][:2], float), 0.0); t = np.array(t)
        for n, F in enumerate(fr_json["frames"]):
            j = int(np.abs(t - (F["t"] - 50)).argmin()); f = fr[j]
            sites = []
            for xx, yy, q, ang in f:
                if q <= 0: continue
                i, d = nearest(Wm, np.array([xx, yy]))
                if d < 4: sites.append([float(xx), float(yy), float(Wm[3][i, 0]), float(Wm[3][i, 1])])
            if len(sites) >= 2:
                cfgs.append(dict(id=f"meas_{n}", kind="measured", interval=n, t=float(t[j]), A=Am, B=Bm, sites=sites))
    json.dump(cfgs, open(os.path.join(HERE, f"cfg_{tag}.json"), "w"))
    print(tag, len(cfgs), "configs,", sum(len(c["sites"]) for c in cfgs), "solves; measured wall counts:",
          [len(c["sites"]) for c in cfgs if c["kind"] == "measured"][:40])
