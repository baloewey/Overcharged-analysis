#!/usr/bin/env python3
"""Kink statistics in the shape sweep (user, 2026-10-22). All frames (2 t.u.); wall +1/2 (depth < 4) tracked frame to frame (mutual nearest neighbours < 2 cells).
Per frame: majority tilt sign of the wall layer; a defect is WRONG if its tilt sign differs. Per run: fraction of wrong defects; episodes = maximal stretches
of consecutive frames a tracked resident is wrong (durations, t.u.; episodes cut by track loss are kept and flagged); episode rate per defect per 1000 t.u.;
tilt histogram (signed relative to the majority sign, 5-deg bins); per-frame C, majority fraction, number of sign changes around the ring.
-> kink_life.npz"""
import glob, json, os, sys, numpy as np
from multiprocessing import Pool
from scipy.spatial import cKDTree
ROOT = "/home/baloewe/Physics/Overcharged/runs"; sys.path.insert(0, ROOT + "/cavity")
from edge_criteria import read_with_angle  # noqa
from wall_fourier import wall_polyline  # noqa
SW = ROOT + "/shapes_server"; H = os.path.dirname(os.path.abspath(__file__)); HB = np.arange(-180, 181, 5)
def one(rel):
    run = os.path.join(SW, rel)
    try:
        bc = json.load(open(run + "/input.json"))["BC"][0]; Q = bc["Q"]; W, s, P = wall_polyline(bc); tree = cKDTree(W); s = np.asarray(s)
        Wc = W[:, 0] + 1j * W[:, 1]; T = np.roll(Wc, -1) - np.roll(Wc, 1); T /= np.abs(T); Nn = -1j * T
        if np.mean((Nn * np.conj(Wc - Wc.mean())).real) < 0: Nn = -Nn
        t, fr = read_with_angle(run, np.array(Q[:2], float), 0.5); prevX = None; prevID = None; nid = 0
        wrong = {}; hist = np.zeros(len(HB) - 1); fstat = []; ntw = 0
        for n, f in enumerate(fr):
            p = f[f[:, 2] > 0] if len(f) else np.zeros((0, 4))
            ids = np.full(len(p), -1)
            if prevX is not None and len(p) and len(prevX):
                dd, jj = cKDTree(prevX).query(p[:, :2]); _, kk = cKDTree(p[:, :2]).query(prevX); ok = (dd < 2) & (kk[jj] == np.arange(len(p))); ids[ok] = prevID[jj[ok]]
            new = ids < 0; ids[new] = nid + np.arange(new.sum()); nid += new.sum(); prevX, prevID = p[:, :2], ids
            if not len(p): continue
            d, j = tree.query(p[:, :2]); w = d < 4
            if w.sum() < 3: continue
            tl = np.angle(np.exp(1j * (p[w, 3] - np.angle(Nn[j[w]])))); sg = np.sign(tl); maj = np.sign(np.sum(sg)) or 1.0
            hist += np.histogram(np.degrees(tl) * maj, HB)[0]; o = np.argsort(s[j[w]]); so = sg[o]
            fstat.append((abs(np.mean(np.exp(1j * tl))), np.mean(sg == maj), np.sum(so != np.roll(so, -1)), w.sum()))
            for i_, wr in zip(ids[w], sg != maj): wrong.setdefault(int(i_), []).append((n, bool(wr)))
            ntw += w.sum()
        eps = []
        for i_, seq in wrong.items():
            cur = 0; last = None
            for n, wr in seq:
                if wr and last is not None and n == last + 1 and cur: cur += 1
                elif wr: 
                    if cur: eps.append(cur)
                    cur = 1
                else:
                    if cur: eps.append(cur)
                    cur = 0
                last = n
            if cur: eps.append(cur)
        dt = (t[1] - t[0]) if len(t) > 1 else 2.0
        return rel, dict(eps=np.array(eps) * dt, rate=len(eps) / max(ntw * dt, 1) * 1000, hist=hist, fstat=np.array(fstat), dt=dt), None
    except Exception as ex:
        return rel, None, str(ex)
if __name__ == "__main__":
    runs = sorted(os.path.relpath(os.path.dirname(f), SW) for f in glob.glob(SW + "/*/*/*/*/s*/DONE"))
    with Pool(int(os.environ.get("NPROC", 8))) as p: res = p.map(one, runs, chunksize=2)
    ok = [(r, a) for r, a, e in res if a is not None]; print("errors", sum(1 for r in res if r[2]), [r[2] for r in res if r[2]][:3], flush=True)
    np.savez_compressed(os.path.join(H, "kink_life.npz"), runs=np.array([r for r, a in ok]), eps=np.array([a["eps"] for r, a in ok], dtype=object),
                        rate=np.array([a["rate"] for r, a in ok]), hist=np.array([a["hist"] for r, a in ok]), fstat=np.array([a["fstat"] for r, a in ok], dtype=object),
                        dt=np.array([a["dt"] for r, a in ok]), HB=HB)
