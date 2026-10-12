#!/usr/bin/env python3
"""Prepare interval-averaged defect force configurations and measured cell velocities (fixed-radius scan, seed 2) for the Stokes fit.
Wall/bulk: +1/2 within 4 cells of the wall = wall. Every 5th defect frame inside each flow interval (tf-50, tf], weight 1/5.
Writes stokes_frames/<run>.json (forces, cells) and <run>_vel.npz (measured velocities, cell centres)."""
import gzip, json, os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..", "..", "..", "runs", "cavity")))
from edge_criteria import read_with_angle
NF = 40
for k in (2, 6, 10):
    run = os.path.abspath(os.path.join(HERE, "..", "..", "..", "..", "runs", "cavity", f"Rscan_R31.83_seeds/s02/k{k:02d}"))
    bc = json.load(open(run + "/input.json"))["BC"][0]; R = bc["R"]; Q = bc["Q"]
    t, fr = read_with_angle(run, np.array(Q[:2], float), 0.0)
    rows = [l.split() for l in gzip.open(run + "/flowfield.dat.gz", "rt") if len(l.split()) == 7 and l.split()[0][0].isdigit()]
    A = np.array(rows, float); ts = np.unique(A[:, 0]); nc = len(A) // len(ts); A = A[:len(ts) * nc].reshape(len(ts), nc, 7)
    x = A[0, :, 1] + 0.5 - Q[0]; y = A[0, :, 2] + 0.5 - Q[1]; ins = np.hypot(x, y) < R - 1
    pick = np.linspace(len(ts) // 5, len(ts) - 1, NF).astype(int)
    frames = []
    for i in pick:
        tf = ts[i]; sel = [j for j in range(len(t)) if tf - 50 < t[j] <= tf][::5]
        w = 1.0 / max(len(sel), 1); wp, bp, mm = [], [], []
        for j in sel:
            for xx, yy, q, ang in fr[j]:
                s = R - np.hypot(xx, yy)
                if q > 0: (wp if s < 4 else bp).append([xx, yy, np.cos(ang), np.sin(ang), w])
                else: mm.append([xx, yy, ang, w])
        frames.append(dict(t=float(tf), wallp=wp, bulkp=bp, minus=mm))
    json.dump(dict(R=R, k=k, cells=np.column_stack([x[ins], y[ins]]).tolist(), frames=frames), open(f"{HERE}/k{k:02d}.json", "w"))
    np.savez(f"{HERE}/k{k:02d}_vel.npz", vx=A[pick][:, ins, 4], vy=A[pick][:, ins, 5], x=x[ins], y=y[ins])
    print(k, len(frames), "frames;", np.mean([len(f["wallp"]) for f in frames]) / 5, "wall +1/2,", np.mean([len(f["bulkp"]) for f in frames]) / 5, "bulk +1/2,",
          np.mean([len(f["minus"]) for f in frames]) / 5, "-1/2 per defect frame (approx)")
