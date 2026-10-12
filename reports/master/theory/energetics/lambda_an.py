#!/usr/bin/env python3
"""Causal test of the tilt protection (outline plan B; runs/lambda, reference lambda 2 = fineflow/fineflow2): per run (second half): wall +1/2 count (< 4 cells), bulk charge,
resident tilt <cos tau> and <|sin tau|>, depleting-escape rate per wall defect per 1000 t.u. (energetics.events), flow rms (cavity) and wall-band tangential rms.
Usage: lambda_an.py TAG RUNDIR [...] -> lambda_<TAG>.txt (per run + grouped by (lambda, eps))"""
import os, sys, json, gzip, numpy as np
from multiprocessing import Pool
from scipy.spatial import cKDTree
from matplotlib.path import Path
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
import energetics as En  # noqa
from selfflow import load_flow  # noqa
def one(run):
    J, bc, k, t, fr, ev = En.events(run, max_events=10 ** 6); lam = J["species"][0].get("tumble", 2.0); eps = float(bc.get("aInv", [1])[0] or 1); sg = 1 if k > 0 else -1
    W, s, P = En.wall_polyline(bc); tree = cKDTree(W); Wc = W[:, 0] + 1j * W[:, 1]; Tn = np.roll(Wc, -1) - np.roll(Wc, 1); Tn /= np.abs(Tn); Nn = -1j * Tn
    if np.mean((Nn * np.conj(Wc)).real) < 0: Nn = -Nn; Tn = -Tn
    h = len(fr) // 2; nw, qb, cs, ss = [], [], [], []
    for f in fr[h:]:
        if not len(f): continue
        d, j = tree.query(f[:, :2]); w = (np.sign(f[:, 2]) == sg) & (d < 4); nw.append(w.sum()); qb.append(0.5 * np.sign(f[~w & (d >= 4), 2]).sum() if (~w).any() else 0)
        ang = np.angle(np.exp(1j * (f[w, 3] - np.angle(Nn[j[w]])))); cs += list(np.cos(ang)); ss += list(np.abs(np.sin(ang)))
    T = t[-1] - t[h]; ne = sum(1 for e in ev if t[e["n"]] >= t[h]); rate = 1000 * ne / (np.mean(nw) * T)
    rows, ft, idx = load_flow(run); Q = bc["Q"]; X = rows[:, 1] + 0.5 - Q[0]; Y = rows[:, 2] + 0.5 - Q[1]; ins = Path(W).contains_points(np.c_[X, Y]); dW, jW = tree.query(np.c_[X, Y])
    vr, vt = [], []
    for tt in ft[len(ft) // 2:]:
        ii = idx[tt]; m = ins[ii]; vr.append(np.sqrt(np.mean(rows[ii, 4][m] ** 2 + rows[ii, 5][m] ** 2)))
        b = m & (dW[ii] > 1) & (dW[ii] < 5); vt.append(np.mean(rows[ii, 4][b] * Tn[jW[ii][b]].real + rows[ii, 5][b] * Tn[jW[ii][b]].imag))
    return (os.path.basename(os.path.dirname(run)) + "/" + os.path.basename(run), lam, eps, k, np.mean(nw), np.mean(qb), np.mean(cs), np.mean(ss), ne, rate, np.mean(vr), np.sqrt(np.mean(np.square(vt))))
if __name__ == "__main__":
    tag = sys.argv[1]; runs = [os.path.abspath(r) for r in sys.argv[2:]]
    with Pool(min(16, len(runs))) as pool: res = pool.map(one, runs, chunksize=1)
    L = ["run  lambda  eps  k  wall+1/2  bulk_charge  <cos tau>  <|sin tau|>  n_depl  rate/1000tu  v_rms  wall_vt_rms"]
    for r in res: L.append(f"{r[0]:32s} {r[1]:4.1f} {r[2]:3.1f} {r[3]:2d}  {r[4]:5.2f}  {r[5]:+5.2f}  {r[6]:+.3f}  {r[7]:.3f}  {r[8]:3d}  {r[9]:.3f}  {r[10]:.4f}  {r[11]:.4f}")
    L.append("\nGROUPED (mean ± s.e. over runs): lambda eps | n_runs | wall+1/2 | bulk charge | <cos tau> | rate | v_rms | wall_vt_rms")
    for key in sorted(set((r[1], r[2]) for r in res)):
        S = np.array([r[4:] for r in res if (r[1], r[2]) == key], float); m = S.mean(0); e = S.std(0) / np.sqrt(max(len(S) - 1, 1))
        L.append(f"  λ {key[0]:3.1f} ε {key[1]:3.1f} | {len(S)} | {m[0]:5.2f}±{e[0]:.2f} | {m[1]:+.2f}±{e[1]:.2f} | {m[2]:+.3f}±{e[2]:.3f} | {m[5]:.3f}±{e[5]:.3f} ({int(S[:, 4].sum())} ev) | {m[6]:.4f} | {m[7]:.4f}")
    print("\n".join(L)); open(os.path.join(H, f"lambda_{tag}.txt"), "w").write("\n".join(L) + "\n")
