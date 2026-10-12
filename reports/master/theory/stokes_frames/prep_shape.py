#!/usr/bin/env python3
"""Interval-averaged defect forces and measured cell velocities for ellipse runs (W3, +k; flow every 100 time units).
Wall = within 4 cells of the wall polyline; every 10th defect frame in each flow interval (tf-100, tf], weight 1/n."""
import gzip, json, os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); RUNS = os.path.abspath(os.path.join(HERE, "..", "..", "..", "..", "runs"))
sys.path.insert(0, os.path.join(RUNS, "cavity"))
from edge_criteria import read_with_angle
from wall_fourier import wall_polyline, project
NF = 40
for tag, rel in (("e12", "cassini/W3_ellipse/eps1.2_k06/s01"), ("e15", "cassini/W3_ellipse/eps1.5_k06/s01"), ("e20", "cassini/W3_ellipse/eps2.0_k06/s01")):
    run = os.path.join(RUNS, rel); bc = json.load(open(run + "/input.json"))["BC"][0]; Q = bc["Q"]
    W, s, P = wall_polyline(bc); A, B = np.abs(W[:, 0]).max(), np.abs(W[:, 1]).max()
    t, fr = read_with_angle(run, np.array(Q[:2], float), 0.0)
    rows = [l.split() for l in gzip.open(run + "/flowfield.dat.gz", "rt") if len(l.split()) == 7 and l.split()[0][0].isdigit()]
    Ar = np.array(rows, float); ts = np.unique(Ar[:, 0]); nc = len(Ar) // len(ts); Ar = Ar[:len(ts) * nc].reshape(len(ts), nc, 7)
    x = Ar[0, :, 1] + 0.5 - Q[0]; y = Ar[0, :, 2] + 0.5 - Q[1]
    ins = (x / A) ** 2 + (y / B) ** 2 < 1
    dc, _ = project(W, s, np.column_stack([x, y])); ins &= dc > 1.0
    pick = np.linspace(len(ts) // 5, len(ts) - 1, NF).astype(int)
    frames = []
    for i in pick:
        tf = ts[i]; sel = [j for j in range(len(t)) if tf - 100 < t[j] <= tf][::10]; w = 1.0 / max(len(sel), 1)
        wp, bp, mm = [], [], []
        for j in sel:
            f = fr[j]; d, _ = project(W, s, f[:, :2])
            for (xx, yy, q, ang), dd in zip(f, d):
                if q > 0: (wp if dd < 4 else bp).append([xx, yy, np.cos(ang), np.sin(ang), w])
                else: mm.append([xx, yy, ang, w])
        frames.append(dict(t=float(tf), wallp=wp, bulkp=bp, minus=mm))
    json.dump(dict(A=float(A), B=float(B), cells=np.column_stack([x[ins], y[ins]]).tolist(), frames=frames), open(f"{HERE}/{tag}.json", "w"))
    np.savez(f"{HERE}/{tag}_vel.npz", vx=Ar[pick][:, ins, 4], vy=Ar[pick][:, ins, 5], x=x[ins], y=y[ins])
    print(tag, "A", round(A, 2), "B", round(B, 2), len(frames), "intervals; per defect frame: wall +1/2", round(np.mean([len(f["wallp"]) for f in frames]) / 5, 2),
          "bulk +1/2", round(np.mean([len(f["bulkp"]) for f in frames]) / 5, 2), "-1/2", round(np.mean([len(f["minus"]) for f in frames]) / 5, 2))
