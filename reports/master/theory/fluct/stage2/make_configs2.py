#!/usr/bin/env python3
"""Stage 2, overnight extension (2026-10-06): larger disorder ensembles + measured snapshots of the second seed (s02).
Synthetic: eps in {1.0001, 1.2, 1.5, 2.0, 2.5}; tangential disorder sigma_t in {0.5, 1, 1.5, 2, 3, 4, 6} (normal sigma_n = 0.5), and
1/2/3 vacancies with sigma_t = 1; 60 new samples each (seed differs from make_configs.py).
Measured: W3 ellipses seed s02, one snapshot per flow interval (40 intervals, as for s01); the measured flow of each interval is saved as
s02_<tag>_vel.npz (same picking as stokes_frames/prep_shape.py).
All configurations are shuffled into NCHUNK files cfg2_cNN.json (each config carries its tag) for parallel solving."""
import gzip, json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = "/home/baloewe/Physics/Overcharged/runs"
sys.path.insert(0, os.path.join(RUNS, "cavity"))
from edge_criteria import read_with_angle  # noqa
from wall_fourier import wall_polyline  # noqa
ns = {"__file__": os.path.join(HERE, "make_configs.py")}; exec(open(os.path.join(HERE, "make_configs.py")).read().split("out = {}")[0], ns)
wall, at, nearest = ns["wall"], ns["at"], ns["nearest"]
rng = np.random.default_rng(777); NS, NN, DEPTH, NCHUNK, NF = 60, 10, 1.5, 8, 40
cfgs = []
for eps, tag in ((1.0001, "e10"), (1.2, "e12"), (1.5, "e15"), (2.0, "e20"), (2.5, "e25")):
    A = np.sqrt(2500 * eps / np.pi); B = A / eps; Wd = wall(A, B); P = Wd[2]
    base = (np.arange(NN) + 0.5) * P / NN
    if tag == "e25":
        cfgs.append(dict(id="perfect", tag=tag, kind="perfect", A=A, B=B, sites=[at(Wd, sb, DEPTH) for sb in base]))
    for st in (0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0):
        for k in range(NS):
            sa = base + rng.normal(0, st, NN); dp = np.clip(DEPTH + rng.normal(0, 0.5, NN), 0.8, 4)
            cfgs.append(dict(id=f"big_dis{st:g}_{k}", tag=tag, kind="disorder", sigma=st, A=A, B=B, sites=[at(Wd, a, d) for a, d in zip(sa, dp)]))
    for nv in (1, 2, 3):
        for k in range(NS):
            keep = np.sort(rng.choice(NN, NN - nv, replace=False))
            sa = base[keep] + rng.normal(0, 1.0, NN - nv); dp = np.clip(DEPTH + rng.normal(0, 0.5, NN - nv), 0.8, 4)
            cfgs.append(dict(id=f"big_vac{nv}_{k}", tag=tag, kind="vacancy", vac=nv, A=A, B=B, sites=[at(Wd, a, d) for a, d in zip(sa, dp)]))
for tag, rel in (("e12", "cassini/W3_ellipse/eps1.2_k06/s02"), ("e15", "cassini/W3_ellipse/eps1.5_k06/s02"), ("e20", "cassini/W3_ellipse/eps2.0_k06/s02")):
    run = os.path.join(RUNS, rel); bc = json.load(open(run + "/input.json"))["BC"][0]; Q = bc["Q"]
    W, s, P = wall_polyline(bc); A, B = np.abs(W[:, 0]).max(), np.abs(W[:, 1]).max(); Wm = wall(A, B)
    t, fr = read_with_angle(run, np.array(Q[:2], float), 0.0); t = np.array(t)
    rows = [l.split() for l in gzip.open(run + "/flowfield.dat.gz", "rt") if len(l.split()) == 7 and l.split()[0][0].isdigit()]
    Ar = np.array(rows, float); ts = np.unique(Ar[:, 0]); nc = len(Ar) // len(ts); Ar = Ar[:len(ts) * nc].reshape(len(ts), nc, 7)
    x = Ar[0, :, 1] + 0.5 - Q[0]; y = Ar[0, :, 2] + 0.5 - Q[1]
    ins = (x / A) ** 2 + (y / B) ** 2 < 1
    pick = np.linspace(len(ts) // 5, len(ts) - 1, NF).astype(int)
    np.savez(os.path.join(HERE, f"s02_{tag}_vel.npz"), vx=Ar[pick][:, ins, 4], vy=Ar[pick][:, ins, 5], x=x[ins], y=y[ins])
    for n, i in enumerate(pick):
        j = int(np.abs(t - (ts[i] - 50)).argmin()); f = fr[j]; sites = []
        for xx, yy, q, ang in f:
            if q <= 0: continue
            ii, d = nearest(Wm, np.array([xx, yy]))
            if d < 4: sites.append([float(xx), float(yy), float(Wm[3][ii, 0]), float(Wm[3][ii, 1])])
        if len(sites) >= 2:
            cfgs.append(dict(id=f"s02_meas_{n}", tag=tag, kind="measured", seed="s02", interval=n, t=float(t[j]), A=A, B=B, sites=sites))
order = rng.permutation(len(cfgs))
for c in range(NCHUNK):
    json.dump([cfgs[i] for i in order[c::NCHUNK]], open(os.path.join(HERE, f"cfg_2_c{c:02d}.json"), "w"))
print(len(cfgs), "configs,", sum(len(c["sites"]) for c in cfgs), "solves in", NCHUNK, "chunks")
