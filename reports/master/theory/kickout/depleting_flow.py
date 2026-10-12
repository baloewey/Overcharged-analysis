#!/usr/bin/env python3
"""DEPLETING escapes (user, 2026-10-11): escapes that take a resident wall +1/2 out of a layer at or below 2k-2 and leave a persistent vacancy, vs other escapes.
Shape sweep, frames every 2 t.u., second half. Wall count n_w(t) = +1/2 within 4 cells of the wall, smoothed over +-5 frames (20 t.u.) to remove flicker.
Escape events as in channels.py (leaves the band, deeper than 6 cells within 20 t.u.). An escape at frame n is DEPLETING if
  (i) smoothed n_w at n-5 <= 2k-2 (layer not crowded), (ii) RESIDENT: the escaper's backward chain stays in the band (depth < 4) for the previous 20 frames (40 t.u.),
  (iii) no SWAP (no other bulk +1/2 enters the band within 4 cells of the site within 20 t.u.), (iv) smoothed n_w over n+5 .. n+25 (10-50 t.u. after) stays below its value at n-5.
Precursors at frame n (escaper vs other escapes vs random resident wall defects in the same frame, same region): depth, cos(tau), nearest other bulk +1/2 distance and heading,
layer tilt coherence, and the flow at the site from the nearest flowfield frame (|v_n| towards/away from the wall at depth 2-6 within 4 cells along the wall, and |v_t|).
Rates of depleting escapes per wall defect per 1000 t.u., circle-like vs elongated, by R_eq. -> depleting.txt"""
import csv, glob, gzip, json, os, sys
import numpy as np
from multiprocessing import Pool
from scipy.spatial import cKDTree
ROOT = "/home/baloewe/Physics/Overcharged/runs"; sys.path.insert(0, ROOT + "/cavity")
from edge_criteria import read_with_angle  # noqa
from wall_fourier import wall_polyline  # noqa
SW = ROOT + "/shapes_server"; H = os.path.dirname(os.path.abspath(__file__))
def one(rel):
    try:
        run = rel if rel.startswith('/') else os.path.join(SW, rel); J = json.load(open(run + "/input.json")); bc = J["BC"][0]; Q = bc["Q"]; k = bc["customAnch"]; N0 = 2 * k - 2
        W, s, P = wall_polyline(bc); tree = cKDTree(W); A = np.abs(W[:, 0]).max(); Wc = W[:, 0] + 1j * W[:, 1]; Tn = np.roll(Wc, -1) - np.roll(Wc, 1); Nn = -1j * Tn / np.abs(Tn)
        if np.mean((Nn * np.conj(Wc)).real) < 0: Nn = -Nn
        t, fr = read_with_angle(run, np.array(Q[:2], float), 0.5); t = np.array(t); D = [tree.query(f[:, :2])[0] if len(f) else np.zeros(0) for f in fr]
        nw = np.array([((f[:, 2] > 0) & (D[n] < 4)).sum() if len(f) else 0 for n, f in enumerate(fr)], float); nws = np.convolve(nw, np.ones(11) / 11, mode="same")
        # flowfield frames (every 100 t.u.)
        fn = run + "/flowfield.dat.gz" if os.path.exists(run + "/flowfield.dat.gz") else run + "/flowfield.dat"; op = gzip.open if fn.endswith('.gz') else open
        rows = np.loadtxt((l for l in op(fn, "rt") if len(l.split()) == 7 and l.split()[0][0].isdigit()))
        ft = np.unique(rows[:, 0]); flows = {tt: rows[rows[:, 0] == tt] for tt in ft if tt >= t[0] - 100}
        def flow_at(tn, p):
            tt = ft[np.argmin(np.abs(ft - tn))]; r = flows.get(tt)
            if r is None: return np.nan, np.nan
            x = r[:, 1] + 0.5 - Q[0]; y = r[:, 2] + 0.5 - Q[1]; d, j = tree.query(np.column_stack([x, y])); dp, jp = tree.query(p)
            sel = (d > 2) & (d < 6) & (np.abs(((s[j] - s[jp] + P / 2) % P) - P / 2) < 4)
            if not sel.any(): return np.nan, np.nan
            nn = Nn[j[sel]]; vn = r[sel, 4] * nn.real + r[sel, 5] * nn.imag; vt = -r[sel, 4] * nn.imag + r[sel, 5] * nn.real
            return np.mean(np.abs(vn)), np.mean(np.abs(vt))
        def cos_hist(n, p0):
            out = []; c = p0
            for m in range(n, n - 21, -1):
                h = fr[m]; hp = np.flatnonzero(h[:, 2] > 0) if len(h) else np.array([], int)
                if not len(hp): break
                d2 = np.hypot(*(h[hp, :2] - c).T); kk = np.argmin(d2)
                if d2[kk] > 2.5: break
                c = h[hp[kk], :2]
                if (n - m) in (5, 10, 20): _, jj = tree.query(c); out.append(((n - m), np.cos(np.angle(np.exp(1j * (h[hp[kk], 3] - np.angle(Nn[jj])))))))
            dct = dict(out); return [dct.get(5, np.nan), dct.get(10, np.nan), dct.get(20, np.nan)]
        def feats(n, i):
            f = fr[n]; p = f[i, :2]; dp, jp = tree.query(p); tau = np.angle(np.exp(1j * (f[i, 3] - np.angle(Nn[jp]))))
            w = (f[:, 2] > 0) & (D[n] < 4); coh = abs(np.mean(np.exp(1j * (f[w, 3] - np.angle(Nn[tree.query(f[w, :2])[1]])))))
            b = np.flatnonzero((f[:, 2] > 0) & (D[n] > 4)); rn, hc = 15.0, 0.0
            if len(b):
                dv = p - f[b, :2]; r = np.hypot(*dv.T); kk = np.argmin(r); rn = min(r[kk], 15.0); hc = float(np.dot([-np.cos(f[b[kk], 3]), -np.sin(f[b[kk], 3])], dv[kk]) / max(r[kk], 1e-9))
            vn, vt = flow_at(t[n], p)
            fh = [flow_at(t[max(n - dm, 0)], p) for dm in (5, 10, 20)]
            return [dp, np.cos(tau), rn, hc, coh, vn, vt] + cos_hist(n, p) + [x for v in fh for x in v]
        dep, oth, ctl = [], [], []; rng = np.random.default_rng(5)
        for n in range(21, len(fr) - 26):
            f, g = fr[n], fr[n + 1]
            if not len(f) or not len(g): continue
            gp = g[g[:, 2] > 0, :2]
            for i in np.flatnonzero((f[:, 2] > 0) & (D[n] < 4)):
                p = f[i, :2]; dd = np.hypot(*(gp - p).T) if len(gp) else np.array([99.0])
                if dd.min() > 2.5 or tree.query(gp[np.argmin(dd)])[0] < 4: continue
                cur = gp[np.argmin(dd)]; ok = True
                for m in range(n + 2, n + 11):
                    h = fr[m]; hp = h[h[:, 2] > 0, :2] if len(h) else np.zeros((0, 2))
                    if not len(hp): ok = False; break
                    d2 = np.hypot(*(hp - cur).T); kk = np.argmin(d2)
                    if d2[kk] > 2.5: ok = False; break
                    cur = hp[kk]
                if not ok or tree.query(cur)[0] <= 6: continue
                # resident? backward chain in the band for 20 frames
                c = p; res = True
                for m in range(n - 1, n - 21, -1):
                    h = fr[m]; hp = h[h[:, 2] > 0, :2] if len(h) else np.zeros((0, 2))
                    if not len(hp): res = False; break
                    d2 = np.hypot(*(hp - c).T); kk = np.argmin(d2)
                    if d2[kk] > 2.5 or tree.query(hp[kk])[0] > 4: res = False; break
                    c = hp[kk]
                swap = False
                for m in range(n + 1, n + 11):
                    h = fr[m]
                    if len(h) and ((h[:, 2] > 0) & (D[m] < 4) & (np.hypot(*(h[:, :2] - p).T) < 4) & (np.hypot(*(h[:, :2] - cur).T) > 4)).any(): swap = True; break
                pre = nws[n - 5]; post = nws[n + 5:n + 26].max()
                is_dep = (pre <= N0 + 0.25) and res and (not swap) and (post < pre - 0.25)
                (dep if is_dep else oth).append(feats(n, i))
                pw = [jj for jj in np.flatnonzero((f[:, 2] > 0) & (D[n] < 4)) if jj != i and abs(abs(f[jj, 0]) - abs(p[0])) / A < 0.1]
                if pw: ctl.append(feats(n, rng.choice(pw)))
        T = (len(fr) - 47) * 2.0
        return rel, dict(dep=np.array(dep).reshape(-1, 16), oth=np.array(oth).reshape(-1, 16), ctl=np.array(ctl).reshape(-1, 16), T=T, nw=nw[21:-26].mean(), N0=N0)
    except Exception as ex:
        return rel, None
if __name__ == "__main__":
    man = {r["dir"]: r for r in csv.DictReader(open(SW + "/manifest_shape_sweep.tsv"), delimiter="\t")}
    runs = sorted(os.path.relpath(os.path.dirname(f), SW) for f in glob.glob(SW + "/*/*/*/*/s*/DONE"))
    with Pool(8) as pool: res = pool.map(one, runs, chunksize=4)
    L = []; say = lambda x: (L.append(x), print(x, flush=True)); names = ["depth", "cos tau", "nearest bulk +1/2", "its heading", "layer coherence", "|v_n| at site", "|v_t| at site", "cos tau 10 t.u. before", "cos tau 20 t.u. before", "cos tau 40 t.u. before", "|v_n| 10 before", "|v_t| 10 before", "|v_n| 20 before", "|v_t| 20 before", "|v_n| 40 before", "|v_t| 40 before"]
    for el in (False, True):
        for R in (1.0, 1.5, 2.0, None):
            S = [o for rel, o in res if o and (float(man[rel]["aspect"]) >= 1.25) == el and (R is None or float(man[rel]["Req_la"]) == R)]
            wt = sum(o["nw"] * o["T"] for o in S) / 1000; nd = sum(len(o["dep"]) for o in S); no = sum(len(o["oth"]) for o in S)
            say(f"{'elongated' if el else 'circle-like':11s} R_eq {R if R else 'all'}: depleting escapes {nd / wt:.3f} per wall defect per 1000 t.u. ({nd}; {nd / max(nd + no, 1):.0%} of all escapes)")
        S = [o for rel, o in res if o and (float(man[rel]["aspect"]) >= 1.25) == el]
        Dp = np.concatenate([o["dep"] for o in S]); Ot = np.concatenate([o["oth"] for o in S]); Ct = np.concatenate([o["ctl"] for o in S])
        for q, nm in enumerate(names):
            say(f"     {nm:18s}: depleting {np.nanmean(Dp[:, q]):+.3f}   other escapes {np.nanmean(Ot[:, q]):+.3f}   resident controls {np.nanmean(Ct[:, q]):+.3f}")
    open(os.path.join(H, "depleting2.txt"), "w").write("\n".join(L) + "\n")
