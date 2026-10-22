#!/usr/bin/env python3
"""Per-event escape data with the predictive tracker (2026-10-22; replaces the outcomes of escape_events.py, whose 10 t.u. nearest-neighbour chain
(< 6 cells per 10 t.u.) missed fast escapes, most in crowded layers). Starts: every wall +1/2 (d < 4) at frames every 10 t.u. (as before). Features as
escape_events.py (tilt from the outward normal, arc gaps to the two wall neighbours, depth, |curvature| R_eq, wall count). Outcomes from the start's track
(crowd_protect/exchange/exchange_count2.tracks, all frames, gate 9, 2-frame gaps): esc = deeper than 6 cells at +20 t.u.; late = in the wall (< 4) at +40 t.u.
and deeper than 6 at +60 t.u. -> events2.npz (same columns as events.npz) + runs.csv order unchanged."""
import csv, glob, json, os, sys, numpy as np
from multiprocessing import Pool
from scipy.spatial import cKDTree
sys.path.insert(0, "/home/baloewe/Physics/Overcharged/reports/master/theory/crowd_protect/exchange"); import exchange_count2 as X2  # noqa
ROOT = "/home/baloewe/Physics/Overcharged/runs"; sys.path.insert(0, ROOT + "/cavity")
from edge_criteria import read_with_angle  # noqa
from wall_fourier import wall_polyline  # noqa
H = os.path.dirname(os.path.abspath(__file__)); SW = ROOT + "/shapes_server"
def one(rel):
    run = os.path.join(SW, rel)
    try:
        bc = json.load(open(run + "/input.json"))["BC"][0]; Q = bc["Q"]; W, s, P = wall_polyline(bc); n_ = len(W); Wc = W[:, 0] + 1j * W[:, 1]
        T = np.roll(Wc, -1) - np.roll(Wc, 1); T /= np.abs(T); Nn = -1j * T; kap = np.gradient(np.unwrap(np.angle(Nn)), P / n_)
        if np.mean((Nn * np.conj(Wc - Wc.mean())).real) < 0: Nn = -Nn
        Req = np.sqrt(abs(0.5 * np.sum(W[:, 0] * np.roll(W[:, 1], -1) - np.roll(W[:, 0], -1) * W[:, 1])) / np.pi); tree = cKDTree(W); sarc = np.asarray(s)
        t, fr = read_with_angle(run, np.array(Q[:2], float), 0.5); nf = len(fr)
        TR = [tr for tr in X2.tracks(fr, 9.0) if tr["sg"] > 0]; at = {}
        for k, tr in enumerate(TR):
            for n, xy in zip(tr["n"], tr["xy"]): at[(n, round(xy[0], 3), round(xy[1], 3))] = k
        dep = [dict(zip(tr["n"], tree.query(np.array(tr["xy"]))[0])) for tr in TR]
        out = []
        for a in range(0, nf - 30, 5):
            f = fr[a]
            if not len(f): continue
            p = f[f[:, 2] > 0]; d, j = tree.query(p[:, :2]); w = d < 4
            if w.sum() < 2: continue
            p0 = p[w]; i0 = j[w]; sa = sarc[i0]; o = np.argsort(sa); ss = sa[o]; g = np.diff(np.r_[ss, ss[0] + P]); gl = np.empty(len(ss)); gr = np.empty(len(ss)); gr[o] = g; gl[o] = np.roll(g, 1)
            tl = np.degrees(np.angle(np.exp(1j * (p0[:, 3] - np.angle(Nn[i0])))))
            for q, (x, y) in enumerate(p0[:, :2]):
                k = at.get((a, round(x, 3), round(y, 3)))
                if k is None: continue
                dd = dep[k]; d10, d20, d30 = dd.get(a + 10), dd.get(a + 20), dd.get(a + 30)
                esc = d10 is not None and d10 > 6; late = d20 is not None and d30 is not None and d20 < 4 and d30 > 6
                out.append((a, tl[q], gl[q], gr[q], d[w][q], abs(kap[i0[q]]) * Req, esc, w.sum(), late))
        return rel, np.array(out) if out else np.zeros((0, 9)), None
    except Exception as ex:
        return rel, None, str(ex)
if __name__ == "__main__":
    meta = list(csv.DictReader(open(os.path.join(H, "runs.csv")))); runs = [m["run"] for m in meta]
    with Pool(int(os.environ.get("NPROC", 8))) as pool: res = pool.map(one, runs, chunksize=2)
    E = np.vstack([np.column_stack([np.full(len(ev), j), ev]) for j, (rel, ev, err) in enumerate(res) if ev is not None and len(ev)])
    np.savez_compressed(os.path.join(H, "events2.npz"), E=E, cols=np.array(["run", "frame", "tilt", "gap_l", "gap_r", "depth", "kapR", "esc", "nwall", "late"]))
    print(f"{len(E)} events, {int(E[:, 7].sum())} escapes, {int(E[:, 9].sum())} lagged; errors {sum(1 for r in res if r[2])}", flush=True)
