#!/usr/bin/env python3
"""Refinements of channels.py (user, 2026-10-11; cf. Kozhukhov, Loewe, Thijssen, Shendruk, arXiv:2503.19664, 'topological kicks'):
(1) IDENTITY of the bulk +1/2 that approaches a wall +1/2 before it escapes; (2) DIRECTION: who closes the distance.
Tracking: all defects, frames every 2 t.u. (second half), same-sign nearest-neighbour links (< 2.5 cells, greedy by distance); ids carry birth time, birth type
(PAIR: an opposite-sign defect is born within 3 cells in the same frame; WALL: born in the wall band (depth < 4); OTHER) and whether the id has ever been a wall defect.
Escape events as in channels.py (wall +1/2 leaving the band, deeper than 6 cells within 20 t.u.). For each escape, the closest bulk +1/2 (depth > 4) over the preceding 10 frames;
if closer than 3 cells it is the 'kicker', classified: DYNAMICAL (pair-born, age < 200 t.u.), ESCAPED (was a wall defect earlier), PERSISTENT (age > 500 t.u., never at the wall),
OTHER. Control: the same classification of the closest bulk +1/2 to a random wall +1/2 (if within 3 cells). Direction: over the 10 frames before the closest approach, the
change of separation split into the bulk +1/2's and the wall defect's displacement along their connecting line (fraction of the closing done by the bulk +1/2). -> kick_identity.txt"""
import csv, glob, json, os, sys
import numpy as np
from multiprocessing import Pool
from scipy.spatial import cKDTree
ROOT = "/home/baloewe/Physics/Overcharged/runs"; sys.path.insert(0, ROOT + "/cavity")
from edge_criteria import read_with_angle  # noqa
from wall_fourier import wall_polyline  # noqa
SW = ROOT + "/shapes_server"; H = os.path.dirname(os.path.abspath(__file__))
def track(fr, D, t, maxgap=3):
    """ids[n] = ids of the defects of frame n; tracks may resume after up to maxgap missed frames (distance < 2.5 * frames elapsed, capped at 5).
    info[id] = dict(birth, type, sign, wall_ever)."""
    ids, info, nxt, last = [], {}, 0, {}          # last[id] = (frame, x, y, sign)
    for n, f in enumerate(fr):
        cur = -np.ones(len(f), int)
        cand = [(k, v) for k, v in last.items() if n - v[0] <= maxgap]
        if len(f) and cand:
            for sg in (1, -1):
                a = np.flatnonzero(np.sign(f[:, 2]) == sg); b = [(k, v) for k, v in cand if v[3] == sg]
                if not len(a) or not b: continue
                bx = np.array([[v[1], v[2]] for _, v in b]); gap = np.array([n - v[0] for _, v in b])
                d = np.hypot(f[a, None, 0] - bx[None, :, 0], f[a, None, 1] - bx[None, :, 1]); lim = np.minimum(2.5 * gap, 5.0)[None, :]
                d = np.where(d <= lim, d + 0.01 * gap[None, :], np.inf); used_a, used_b = set(), set()
                for flat in np.argsort(d, axis=None):
                    i, j = divmod(flat, len(b))
                    if not np.isfinite(d[i, j]): break
                    if i in used_a or j in used_b: continue
                    used_a.add(i); used_b.add(j); cur[a[i]] = b[j][0]
        new = np.flatnonzero(cur < 0)
        for i in new:
            cur[i] = nxt; typ = "OTHER"
            if n > 0:
                if D[n][i] < 4: typ = "WALL"
                opp = [j for j in new if j != i and np.sign(f[j, 2]) != np.sign(f[i, 2]) and np.hypot(*(f[j, :2] - f[i, :2])) < 3]
                if opp: typ = "PAIR"
            info[nxt] = dict(birth=t[n], type=typ, sign=int(np.sign(f[i, 2])), wall_ever=False); nxt += 1
        for i in range(len(f)):
            last[cur[i]] = (n, f[i, 0], f[i, 1], int(np.sign(f[i, 2])))
            if f[i, 2] > 0 and D[n][i] < 4: info[cur[i]]["wall_ever"] = True
        ids.append(cur)
    return ids, info
def classify(info, idx, tnow, wall_before, life=None):
    inf = info[idx]; age = tnow - inf["birth"]
    if life is not None and life.get(idx, 0) < 10: return "TRANSIENT"
    if wall_before: return "ESCAPED"
    if inf["type"] == "PAIR" and age < 200: return "DYNAMICAL"
    if age > 500: return "PERSISTENT"
    return "OTHER"
def one(rel):
    try:
        run = os.path.join(SW, rel); bc = json.load(open(run + "/input.json"))["BC"][0]; Q = bc["Q"]
        W, s, P = wall_polyline(bc); tree = cKDTree(W); t, fr = read_with_angle(run, np.array(Q[:2], float), 0.5); t = np.array(t)
        def merge(f):
            if len(f) < 2: return f
            keep = np.ones(len(f), bool); out = []
            for i in range(len(f)):
                if not keep[i]: continue
                grp = [j for j in range(i, len(f)) if keep[j] and np.sign(f[j, 2]) == np.sign(f[i, 2]) and np.hypot(*(f[j, :2] - f[i, :2])) < 1.5]
                for j in grp: keep[j] = False
                row = f[i].copy(); row[:2] = f[grp, :2].mean(0); out.append(row)
            return np.array(out)
        fr = [merge(f) for f in fr]
        D = [tree.query(f[:, :2])[0] if len(f) else np.zeros(0) for f in fr]; ids, info = track(fr, D, t); rng = np.random.default_rng(2)
        life = {}
        for n_, a_ in enumerate(ids):
            for i_ in a_: life[i_] = life.get(i_, 0) + 1
        # per id: frames at which it was a wall defect (for 'escaped' status at a given time)
        wall_frames = {}
        for n, f in enumerate(fr):
            for i in range(len(f)):
                if f[i, 2] > 0 and D[n][i] < 4: wall_frames.setdefault(ids[n][i], []).append(n)
        out = dict(kick=[], ctrl=[], frac=[], esc=0)
        for n in range(11, len(fr) - 11):
            f, g = fr[n], fr[n + 1]
            if not len(f) or not len(g): continue
            pw = np.flatnonzero((f[:, 2] > 0) & (D[n] < 4))
            for i in pw:
                me = ids[n][i]; nxt_i = np.flatnonzero(ids[n + 1] == me)
                if not len(nxt_i) or D[n + 1][nxt_i[0]] < 4: continue
                # escape if the same id is deeper than 6 within 10 frames
                deep = False
                for m in range(n + 2, n + 12):
                    k = np.flatnonzero(ids[m] == me)
                    if not len(k): break
                    if D[m][k[0]] > 6: deep = True; break
                if not deep: continue
                out["esc"] += 1
                recent = any(len(np.flatnonzero(ids[m] == me)) and D[m][np.flatnonzero(ids[m] == me)[0]] > 6 for m in range(max(n - 20, 0), n))
                out.setdefault("recent", []).append(bool(recent))
                jc = rng.choice(pw); mc = ids[n][jc]
                out.setdefault("recent_ctrl", []).append(bool(any(len(np.flatnonzero(ids[m] == mc)) and D[m][np.flatnonzero(ids[m] == mc)[0]] > 6 for m in range(max(n - 20, 0), n))))
                def kicker(pidx, n0):
                    best = (99.0, None, None)
                    for m in range(max(n0 - 20, 0), n0 + 1):
                        h = fr[m]; kk = np.flatnonzero(ids[m] == pidx)
                        if not len(h) or not len(kk): continue
                        sel = np.flatnonzero((h[:, 2] > 0) & (D[m] > 4) & (ids[m] != pidx))
                        if not len(sel): continue
                        dd = np.hypot(*(h[sel, :2] - h[kk[0], :2]).T); j = np.argmin(dd)
                        if dd[j] < best[0]: best = (dd[j], ids[m][sel[j]], m)
                    return best
                dmin, kid, mk = kicker(me, n)
                if dmin < 3:
                    wb = any(x < mk for x in wall_frames.get(kid, []))
                    out["kick"].append(classify(info, kid, t[mk], wb, life))
                    # direction: approach phase (mk-5 -> mk) and departure phase (mk -> n+5); signed contributions along the connecting line
                    def contrib(m0, m1):
                        a0 = np.flatnonzero(ids[m0] == me); b0 = np.flatnonzero(ids[m0] == kid); a1 = np.flatnonzero(ids[m1] == me); b1 = np.flatnonzero(ids[m1] == kid)
                        if not (len(a0) and len(b0) and len(a1) and len(b1)): return None
                        w0, bp0, w1, bp1 = fr[m0][a0[0], :2], fr[m0][b0[0], :2], fr[m1][a1[0], :2], fr[m1][b1[0], :2]
                        u = (bp0 - w0) / max(np.hypot(*(bp0 - w0)), 1e-9)
                        return (-np.dot(bp1 - bp0, u), np.dot(w1 - w0, u))      # >0 = moving towards the other
                    ca = contrib(max(mk - 5, 0), mk); cd = contrib(mk, min(n + 5, len(fr) - 1))
                    out["frac"].append((ca, cd, mk - n))
                j = rng.choice(pw); dmin, kid, mk = kicker(ids[n][j], n)
                if dmin < 3:
                    wb = any(x < mk for x in wall_frames.get(kid, [])); out["ctrl"].append(classify(info, kid, t[mk], wb, life))
        out["life"] = list(life.values()); return rel, out
    except Exception as ex:
        return rel, None
if __name__ == "__main__":
    man = {r["dir"]: r for r in csv.DictReader(open(SW + "/manifest_shape_sweep.tsv"), delimiter="\t")}
    runs = sorted(os.path.relpath(os.path.dirname(f), SW) for f in glob.glob(SW + "/*/*/*/*/s*/DONE"))
    with Pool(8) as pool: res = pool.map(one, runs, chunksize=4)
    L = []; say = lambda x: (L.append(x), print(x, flush=True)); cats = ("TRANSIENT", "DYNAMICAL", "ESCAPED", "PERSISTENT", "OTHER")
    for el in (False, True):
        for R in (1.0, 1.5, 2.0):
            S = [o for rel, o in res if o and (float(man[rel]["aspect"]) >= 1.25) == el and float(man[rel]["Req_la"]) == R]
            K = sum((o["kick"] for o in S), []); C = sum((o["ctrl"] for o in S), []); F = sum((o["frac"] for o in S), []); E = sum(o["esc"] for o in S)
            ks = "  ".join(f"{c} {K.count(c) / max(len(K), 1):.0%}" for c in cats); cs = "  ".join(f"{c} {C.count(c) / max(len(C), 1):.0%}" for c in cats)
            A = np.array([x[0] for x in F if x[0] is not None]); Dp = np.array([x[1] for x in F if x[1] is not None]); lag = np.array([x[2] for x in F]) * 2.0
            Rc = sum((o.get("recent", []) for o in S), []); Rcc = sum((o.get("recent_ctrl", []) for o in S), [])
            say(f"{'elongated' if el else 'circle-like':11s} R_eq {R}: {E} escapes; escaper was in the bulk within 40 t.u. before: {np.mean(Rc):.0%} (random wall defect: {np.mean(Rcc):.0%}) | "
                f"other bulk +1/2 within 3 cells: {len(K) / max(E, 1):.0%} of escapes vs {len(C) / max(E, 1):.0%} control | kicker: {ks} | control: {cs}")
            if len(A) and len(Dp):
                say(f"      closest approach {np.median(-lag):.0f} t.u. before the escape (median); approach phase (10 t.u.): bulk +1/2 closes {A[:, 0].mean():+.2f}, wall defect closes {A[:, 1].mean():+.2f} cells;"
                    f" departure (closest approach -> escape + 10 t.u.): bulk +1/2 {-Dp[:, 0].mean():+.2f}, wall defect {-Dp[:, 1].mean():+.2f} cells away from the other")
    open(os.path.join(H, "kick_identity.txt"), "w").write("\n".join(L) + "\n")
