#!/usr/bin/env python3
"""Range-resolved approach statistics (self-match excluded). Shape sweep, frames every 2 t.u. (second half).
Escape events as in channels.py. For each escape (and a random wall +1/2 at the same frame = control) the defect is followed BACKWARDS 10 frames along its own nearest +1/2 chain
(< 2.5 cells per frame); at each of these frames, all OTHER bulk +1/2 (depth > 4, not the chain member) within 15 cells are histogrammed by distance r (1-cell bins) and, for r < 10,
the cosine between the neighbour's head direction (-angle, extensile) and the direction from the neighbour to the wall defect (> 0: heading towards it).
Output per class: ratio of neighbour density escape/control vs r; mean heading cosine (escape vs control) for r < 5 and 5-10. -> approach_range.txt"""
import csv, glob, json, os, sys
import numpy as np
from multiprocessing import Pool
from scipy.spatial import cKDTree
ROOT = "/home/baloewe/Physics/Overcharged/runs"; sys.path.insert(0, ROOT + "/cavity")
from edge_criteria import read_with_angle  # noqa
from wall_fourier import wall_polyline  # noqa
SW = ROOT + "/shapes_server"; H = os.path.dirname(os.path.abspath(__file__)); RB = np.arange(0, 16, 1.0)
def back_chain(fr, n, p, nback=10):
    pos = {n: p}; cur = p
    for m in range(n - 1, n - 1 - nback, -1):
        h = fr[m]; hp = h[h[:, 2] > 0, :2] if len(h) else np.zeros((0, 2))
        if not len(hp): break
        d = np.hypot(*(hp - cur).T); k = np.argmin(d)
        if d[k] > 2.5: break
        cur = hp[k]; pos[m] = cur
    return pos
TB = np.arange(-12, 12.1, 3.0); DB = np.array([4, 6, 8, 11, 15])
def neigh(fr, D, pos):
    hist = np.zeros(len(RB) - 1); cos_near, cos_far = [], []; H2 = np.zeros((len(TB) - 1, len(DB) - 1))
    for m, w in pos.items():
        h = fr[m]
        if not len(h): continue
        sel = (h[:, 2] > 0) & (D[m] > 4) & (np.hypot(*(h[:, :2] - w).T) > 0.5)
        if not sel.any(): continue
        dv = w - h[sel, :2]; r = np.hypot(*dv.T); hist += np.histogram(r, RB)[0]
        rad = w / max(np.hypot(*w), 1e-9); tan = np.array([-rad[1], rad[0]]); off = -(dv @ tan); H2 += np.histogram2d(off[r < 12], D[m][sel][r < 12], [TB, DB])[0]
        head = np.column_stack([-np.cos(h[sel, 3]), -np.sin(h[sel, 3])]); c = np.sum(head * dv, 1) / np.maximum(r, 1e-9)
        cos_near += list(c[r < 5]); cos_far += list(c[(r >= 5) & (r < 10)])
    return hist, cos_near, cos_far, H2
def one(rel):
    try:
        run = os.path.join(SW, rel); bc = json.load(open(run + "/input.json"))["BC"][0]; Q = bc["Q"]
        W, s, P = wall_polyline(bc); tree = cKDTree(W); t, fr = read_with_angle(run, np.array(Q[:2], float), 0.5); rng = np.random.default_rng(3)
        D = [tree.query(f[:, :2])[0] if len(f) else np.zeros(0) for f in fr]
        Me = np.zeros((len(TB) - 1, len(DB) - 1)); Mc = Me.copy(); He, Hc = np.zeros(len(RB) - 1), np.zeros(len(RB) - 1); ce_n, ce_f, cc_n, cc_f = [], [], [], []; ne = nc = 0
        for n in range(11, len(fr) - 11):
            f, g = fr[n], fr[n + 1]
            if not len(f) or not len(g): continue
            pw = np.flatnonzero((f[:, 2] > 0) & (D[n] < 4)); gp = g[g[:, 2] > 0, :2]
            for i in pw:
                p = f[i, :2]; dd = np.hypot(*(gp - p).T) if len(gp) else np.array([99.0])
                if dd.min() > 2.5 or tree.query(gp[np.argmin(dd)])[0] < 4: continue
                cur = gp[np.argmin(dd)]; ok = True
                for m in range(n + 2, n + 11):
                    h = fr[m]; hp = h[h[:, 2] > 0, :2] if len(h) else np.zeros((0, 2))
                    if not len(hp): ok = False; break
                    d2 = np.hypot(*(hp - cur).T); k = np.argmin(d2)
                    if d2[k] > 2.5: ok = False; break
                    cur = hp[k]
                if not ok or tree.query(cur)[0] <= 6: continue
                h1, a, b, m2 = neigh(fr, D, back_chain(fr, n, p)); He += h1; ce_n += a; ce_f += b; ne += 1; Me = Me + m2
                A_ = np.abs(W[:, 0]).max(); same = [jj for jj in pw if jj != i and abs(abs(f[jj, 0]) - abs(p[0])) / A_ < 0.1]
                if not same: continue
                j = rng.choice(same); h1, a, b, m2 = neigh(fr, D, back_chain(fr, n, f[j, :2])); Hc += h1; cc_n += a; cc_f += b; nc += 1; Mc = Mc + m2
        return rel, dict(Me=Me, Mc=Mc, He=He, Hc=Hc, ne=ne, nc=nc, ce_n=ce_n, ce_f=ce_f, cc_n=cc_n, cc_f=cc_f)
    except Exception:
        return rel, None
if __name__ == "__main__":
    man = {r["dir"]: r for r in csv.DictReader(open(SW + "/manifest_shape_sweep.tsv"), delimiter="\t")}
    runs = sorted(os.path.relpath(os.path.dirname(f), SW) for f in glob.glob(SW + "/*/*/*/*/s*/DONE"))
    with Pool(8) as pool: res = pool.map(one, runs, chunksize=4)
    L = []; say = lambda x: (L.append(x), print(x, flush=True))
    say("neighbour density of other bulk +1/2 around escaping wall defects / around random wall defects, vs distance r (cells); 20 t.u. before the escape")
    for el in (False, True):
        S = [o for rel, o in res if o and (float(man[rel]["aspect"]) >= 1.25) == el]
        He = sum(o["He"] for o in S); Hc = sum(o["Hc"] for o in S); ne = sum(o["ne"] for o in S); nc = sum(o["nc"] for o in S)
        ratio = (He / max(ne, 1)) / np.maximum(Hc / max(nc, 1), 1e-12)
        say(f"{'elongated' if el else 'circle-like':11s} ({ne} escapes): ratio by r: " + " ".join(f"{int(a)}-{int(b)}:{x:.2f}" for a, b, x in zip(RB[:-1], RB[1:], ratio)))
        g = lambda k: np.mean(sum((o[k] for o in S), [])) if sum((o[k] for o in S), []) else np.nan
        Me = sum(o["Me"] for o in S) / max(ne, 1); Mc = sum(o["Mc"] for o in S) / max(nc, 1); rr = Me / np.maximum(Mc, 1e-12)
        say("     local map (escape/control density): rows = depth of the neighbour " + str(list(DB)) + ", columns = offset along the wall " + str(list(TB)))
        for k in range(len(DB) - 1): say(f"       depth {DB[k]}-{DB[k + 1]}: " + " ".join(f"{x:5.2f}" for x in rr[:, k]))
        say(f"     heading cosine towards the wall defect: r < 5: escapes {g('ce_n'):+.3f} vs control {g('cc_n'):+.3f};  5 <= r < 10: escapes {g('ce_f'):+.3f} vs control {g('cc_f'):+.3f}")
    open(os.path.join(H, "approach_map.txt"), "w").write("\n".join(L) + "\n")
