#!/usr/bin/env python3
"""Spontaneous vs forced net flow in winding channels (user, 2026-10-19): why does an imposed pressure flow re-screen the walls while the spontaneous
directed flow of narrow channels does not? Per run (second half, W 20-24): channel-mean <v_x>, rms speed; WALL SHEAR at each wall s(t) = slope of the
x-averaged v_x(y) over 1-4 cells from the wall (sign: along the mean flow), its time mean, rms and fraction of frames with the mean's sign; bulk/wall flow
ratio; wall counts and tilt coherence of the bottom layer (frames C > 0.9). -> flow_compare.txt"""
import os, sys, json, gzip, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, "/home/baloewe/Physics/Overcharged/runs/cavity")
from edge_criteria import read_with_angle  # noqa
from multiprocessing import Pool
R_ = "/home/baloewe/Physics/Overcharged/runs"
RUNS = ["channel_current/W24_mixed_s3", "channel_current/W24_mixed_s4", "channel_current/W24_sym_s5", "channel_current/W24_sym_s6", "channel_curve/W20_sym_s11", "channel_curve/W20_mixed_s12",
        "channel_pressure/mixed_g+0.001_s1", "channel_pressure/mixed_g+0.003_s1", "channel_pressure/sym_g+0.001_s1", "channel_pressure/sym_g+0.001_s2", "channel_pressure/sym_g+0.003_s1", "channel_pressure/mixed_passive_g+0.003"]
def opn(fn): return gzip.open(fn + ".gz", "rt") if os.path.exists(fn + ".gz") else open(fn)
def one(r):
    run = f"{R_}/{r}"; J = json.load(open(run + "/input.json")); BC = J["BC"]; y0, y1 = BC[2]["Q"][1], BC[3]["Q"][1]; W = y1 - y0; g = (J.get("grav") or [0])[0]
    rows = np.loadtxt((l for l in opn(run + "/flowfield.dat") if len(l.split()) == 7 and l.split()[0][0].isdigit()))
    ts = np.unique(rows[:, 0]); ts = ts[ts > ts.max() / 2]; M, SB, ST, VR = [], [], [], []
    for t in ts:
        r_ = rows[rows[:, 0] == t]; yb = (r_[:, 2] - y0).astype(int); ok = (yb >= 0) & (yb < W)
        prof = np.bincount(yb[ok], r_[ok, 4], minlength=W) / np.maximum(np.bincount(yb[ok], minlength=W), 1)
        M.append(prof.mean()); VR.append(np.sqrt(np.mean(r_[ok, 4] ** 2 + r_[ok, 5] ** 2)))
        yy = np.arange(W) + 0.5; SB.append(np.polyfit(yy[1:4], prof[1:4], 1)[0]); ST.append(-np.polyfit(yy[-4:-1], prof[-4:-1], 1)[0])
    M, SB, ST = np.array(M), np.array(SB), np.array(ST); sg = np.sign(M.mean()) or 1
    t, fr = read_with_angle(run, np.array([0., 0.]), 0.5); C, nb = [], []
    for f in fr:
        if not len(f): continue
        w = (f[:, 2] > 0) & (f[:, 1] - y0 < 4); nb.append(w.sum())
        if w.sum() >= 3: C.append(abs(np.mean(np.exp(1j * (f[w, 3] + np.pi / 2)))))
    return dict(run=r, g=g, W=W, mean=M.mean(), rms_mean=M.std(), vrms=np.mean(VR), sb=np.mean(SB * sg), sb_sd=np.std(SB), sb_frac=np.mean(np.sign(SB) == sg),
                st=np.mean(ST * sg), nb=np.mean(nb), C=np.mean(C) if C else np.nan, C90=np.mean(np.array(C) > 0.9) if C else np.nan)
if __name__ == "__main__":
    with Pool(6) as p: res = p.map(one, RUNS)
    L = ["SPONTANEOUS vs FORCED net flow (second half). wall shear s = d<v_x>/dy over 1-4 cells from the wall, signed along the mean flow (time mean, frame sd, fraction of frames with the mean's sign).",
         "run | g | W | <v_x> | sd of <v_x>(t) | v_rms | <v_x>/v_rms | wall shear bottom: mean, sd, frac same sign | top mean | bottom +1/2 (of 12) | coherence C, frac C>0.9"]
    for r in res:
        L.append(f"{r['run']:36s} | {r['g']:+.3f} | {r['W']:.0f} | {r['mean']:+.3f} | {r['rms_mean']:.3f} | {r['vrms']:.3f} | {abs(r['mean'])/r['vrms']:.2f} | {r['sb']:+.3f} {r['sb_sd']:.3f} {r['sb_frac']:.2f} | {r['st']:+.3f} | {r['nb']:5.2f} | {r['C']:.2f} {r['C90']:.2f}")
    open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "flow_compare.txt"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
