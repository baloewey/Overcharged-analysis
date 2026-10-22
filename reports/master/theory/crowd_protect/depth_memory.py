#!/usr/bin/env python3
"""Test 1a (2026-10-22): (i) how long do depth and tilt fluctuations of resident wall +1/2 last? (ii) does depth or tilt AVERAGED over the 20 t.u. before
predict a lagged escape? Shape sweep, all frames (2 t.u.), +1/2 tracked frame to frame (mutual nearest neighbours < 2 cells).
(i) autocorrelation of depth and cos tilt (deviation from the track mean) along resident stretches (depth < 4 throughout), lags 0-60 t.u.
(ii) starts every 10 t.u.: tracks in the wall (depth < 4) for the previous 20 t.u.; window means of depth and cos tilt; outcome lagged escape (depth < 4 at
+40 t.u., > 6 at +60 t.u., same track) and prompt escape (> 6 at +20 t.u.). -> depth_memory.npz (per-run arrays)"""
import glob, json, os, sys
import numpy as np
from multiprocessing import Pool
from scipy.spatial import cKDTree
ROOT = "/home/baloewe/Physics/Overcharged/runs"; sys.path.insert(0, ROOT + "/cavity")
from edge_criteria import read_with_angle  # noqa
from wall_fourier import wall_polyline  # noqa
SW = ROOT + "/shapes_server"; H = os.path.dirname(os.path.abspath(__file__)); LAGS = np.arange(0, 31)
def one(rel):
    run = os.path.join(SW, rel)
    try:
        bc = json.load(open(run + "/input.json"))["BC"][0]; Q = bc["Q"]; W, s, P = wall_polyline(bc); tree = cKDTree(W)
        Wc = W[:, 0] + 1j * W[:, 1]; T = np.roll(Wc, -1) - np.roll(Wc, 1); T /= np.abs(T); Nn = -1j * T
        if np.mean((Nn * np.conj(Wc - Wc.mean())).real) < 0: Nn = -Nn
        t, fr = read_with_angle(run, np.array(Q[:2], float), 0.5); nf = len(fr)
        # tracks: id per defect per frame
        X, D, C, ID = [], [], [], []; nid = 0
        for n, f in enumerate(fr):
            p = f[f[:, 2] > 0] if len(f) else np.zeros((0, 4)); d, j = tree.query(p[:, :2]) if len(p) else (np.zeros(0), np.zeros(0, int))
            c = np.cos(np.angle(np.exp(1j * (p[:, 3] - np.angle(Nn[j]))))) if len(p) else np.zeros(0)
            ids = np.full(len(p), -1)
            if n and len(p) and len(X[-1]):
                t1 = cKDTree(X[-1]); dd, jj = t1.query(p[:, :2]); t2 = cKDTree(p[:, :2]); _, kk = t2.query(X[-1])
                ok = (dd < 2) & (kk[jj] == np.arange(len(p))); ids[ok] = ID[-1][jj[ok]]
            new = ids < 0; ids[new] = nid + np.arange(new.sum()); nid += new.sum()
            X.append(p[:, :2]); D.append(d); C.append(c); ID.append(ids)
        dep = np.full((nid, nf), np.nan); cos = np.full((nid, nf), np.nan)
        for n in range(nf): dep[ID[n], n] = D[n]; cos[ID[n], n] = C[n]
        # (i) autocorrelation along resident stretches of >= 60 frames
        acd = np.zeros(len(LAGS)); acc_ = np.zeros(len(LAGS)); na = np.zeros(len(LAGS))
        for i in range(nid):
            r = np.isfinite(dep[i]) & (dep[i] < 4)
            if r.sum() < 60: continue
            idx = np.flatnonzero(r); brk = np.flatnonzero(np.diff(idx) > 1); segs = np.split(idx, brk + 1)
            for sg in segs:
                if len(sg) < 60: continue
                a = dep[i, sg] - dep[i, sg].mean(); b = cos[i, sg] - cos[i, sg].mean()
                for q, L in enumerate(LAGS):
                    acd[q] += np.sum(a[:len(a) - L] * a[L:]); acc_[q] += np.sum(b[:len(b) - L] * b[L:]); na[q] += len(a) - L
        # (ii) window means and outcomes
        S = []
        for n in range(10, nf - 30, 5):
            for i in ID[n]:
                w = dep[i, n - 10:n + 1]
                if not (np.all(np.isfinite(w)) and np.all(w < 4)): continue
                late = np.isfinite(dep[i, n + 20]) and dep[i, n + 20] < 4 and np.isfinite(dep[i, n + 30]) and dep[i, n + 30] > 6
                prompt = np.isfinite(dep[i, n + 10]) and dep[i, n + 10] > 6
                S.append((w.mean(), np.nanmean(cos[i, n - 10:n + 1]), dep[i, n], cos[i, n], late, prompt))
        return rel, (acd / np.maximum(na, 1), acc_ / np.maximum(na, 1), np.array(S)), None
    except Exception as ex:
        return rel, None, str(ex)
if __name__ == "__main__":
    runs = sorted(os.path.relpath(os.path.dirname(f), SW) for f in glob.glob(SW + "/*/*/*/*/s*/DONE"))
    with Pool(int(os.environ.get("NPROC", 8))) as p: res = p.map(one, runs, chunksize=2)
    ok = [(r, a) for r, a, e in res if a is not None]; print("errors", sum(1 for r in res if r[2]), [r[2] for r in res if r[2]][:3])
    np.savez_compressed(os.path.join(H, "depth_memory.npz"), runs=np.array([r for r, a in ok]), acd=np.array([a[0] for r, a in ok]), acc=np.array([a[1] for r, a in ok]),
                        S=np.array([a[2] for r, a in ok], dtype=object), LAGS=LAGS * 2.0)
