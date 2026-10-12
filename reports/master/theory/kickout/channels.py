#!/usr/bin/env python3
"""How are wall +1/2 lost? (user's question 2026-10-11: kick-out by bulk defects?) Shape sweep (675 runs), second half, defect frames every 2 t.u.
Wall +1/2 = depth < 4. Each wall +1/2 at frame n is linked to frame n+1 (nearest +1/2 within 2.5 cells). A wall defect with no successor is classified:
  ANNIHILATION: a -1/2 was within 2.5 cells of it at frame n and no -1/2 remains within 2.5 cells of that -1/2's position at n+1;
  otherwise LOST (tracking loss / merging; not interpreted).
A wall defect whose successor chain (over the next 10 frames = 20 t.u.) ends deeper than 6 cells: ESCAPE (counted once, at the first frame it leaves the band).
For escapes and for a control sample (random wall +1/2 at the same frames): distance to the nearest bulk +1/2 and bulk -1/2 (depth > 4) minimised over the preceding 10 frames.
SWAP: an escape after which a bulk +1/2 enters the band within 4 cells of the vacated site within 20 t.u.
Output per class (circle-like aspect <= 1.22 / elongated >= 1.25) and R_eq: rates per wall defect per 1000 t.u. of escape and annihilation; swap fraction; close-approach
fractions (< 3, < 5 cells) for escapes vs control; and the charge balance. -> channels.txt"""
import csv, glob, json, os, sys
import numpy as np
from multiprocessing import Pool
from scipy.spatial import cKDTree
ROOT = "/home/baloewe/Physics/Overcharged/runs"; sys.path.insert(0, ROOT + "/cavity")
from edge_criteria import read_with_angle  # noqa
from wall_fourier import wall_polyline  # noqa
SW = ROOT + "/shapes_server"; H = os.path.dirname(os.path.abspath(__file__))
def one(rel):
    try:
        run = os.path.join(SW, rel); bc = json.load(open(run + "/input.json"))["BC"][0]; Q = bc["Q"]
        W, s, P = wall_polyline(bc); tree = cKDTree(W); t, fr = read_with_angle(run, np.array(Q[:2], float), 0.5); rng = np.random.default_rng(1)
        D = [tree.query(f[:, :2])[0] if len(f) else np.zeros(0) for f in fr]
        n_wall_frames = 0; esc = ann = lost = swap = 0; app_e, app_c = [], []
        for n in range(10, len(fr) - 11):
            f, g = fr[n], fr[n + 1]
            if not len(f): continue
            pw = np.flatnonzero((f[:, 2] > 0) & (D[n] < 4)); n_wall_frames += len(pw)
            gp = g[g[:, 2] > 0, :2] if len(g) else np.zeros((0, 2)); gm = g[g[:, 2] < 0, :2] if len(g) else np.zeros((0, 2))
            fm = f[f[:, 2] < 0, :2]
            for i in pw:
                p = f[i, :2]
                if not len(gp) or np.min(np.hypot(*(gp - p).T)) > 2.5:          # no successor
                    if len(fm):
                        dm = np.hypot(*(fm - p).T); j = np.argmin(dm)
                        if dm[j] < 2.5 and (not len(gm) or np.min(np.hypot(*(gm - fm[j]).T)) > 2.5): ann += 1; continue
                    lost += 1; continue
                # escape: follow 10 frames; count only if it leaves the band at this step (depth at n+1 >= 4)
                k = np.argmin(np.hypot(*(gp - p).T)); dn1 = tree.query(gp[k])[0]
                if dn1 < 4: continue
                cur = gp[k]; ok = True
                for m in range(n + 2, n + 11):
                    h = fr[m]; hp = h[h[:, 2] > 0, :2] if len(h) else np.zeros((0, 2))
                    if not len(hp): ok = False; break
                    dd = np.hypot(*(hp - cur).T); kk = np.argmin(dd)
                    if dd[kk] > 2.5: ok = False; break
                    cur = hp[kk]
                if not ok or tree.query(cur)[0] <= 6: continue
                esc += 1
                def approach(pt):
                    best = [99.0, 99.0]
                    for m in range(n - 10, n + 1):
                        h = fr[m]
                        if not len(h): continue
                        b = D[m] > 4
                        for c, sg in ((0, 1), (1, -1)):
                            sel = b & (np.sign(h[:, 2]) == sg)
                            if sel.any(): best[c] = min(best[c], np.min(np.hypot(*(h[sel, :2] - pt).T)))
                    return best
                app_e.append(approach(p)); j = rng.choice(pw); app_c.append(approach(f[j, :2]))
                for m in range(n + 1, n + 11):          # swap: a bulk +1/2 enters the band near the vacated site
                    h = fr[m]
                    if not len(h): continue
                    sel = (h[:, 2] > 0) & (D[m] < 4) & (np.hypot(*(h[:, :2] - p).T) < 4) & (np.hypot(*(h[:, :2] - cur).T) > 4)
                    if sel.any(): swap += 1; break
        T = (len(fr) - 21) * 2.0
        return rel, dict(T=T, nw=n_wall_frames / max(len(fr) - 21, 1), esc=esc, ann=ann, lost=lost, swap=swap, app_e=np.array(app_e).reshape(-1, 2), app_c=np.array(app_c).reshape(-1, 2))
    except Exception as ex:
        return rel, None
if __name__ == "__main__":
    man = {r["dir"]: r for r in csv.DictReader(open(SW + "/manifest_shape_sweep.tsv"), delimiter="\t")}
    runs = sorted(os.path.relpath(os.path.dirname(f), SW) for f in glob.glob(SW + "/*/*/*/*/s*/DONE"))
    with Pool(8) as pool: res = pool.map(one, runs, chunksize=4)
    L = []; say = lambda x: (L.append(x), print(x, flush=True))
    say("rates per wall defect per 1000 t.u.; approach = min distance of the nearest BULK defect over the 20 t.u. before the event (escape) vs a random wall defect at the same frame")
    for el in (False, True):
        for R in (1.0, 1.5, 2.0):
            S = [r for rel, r in res if r and (float(man[rel]["aspect"]) >= 1.25) == el and float(man[rel]["Req_la"]) == R]
            wt = sum(r["nw"] * r["T"] for r in S) / 1000; E = sum(r["esc"] for r in S); A = sum(r["ann"] for r in S); Lo = sum(r["lost"] for r in S); Sw = sum(r["swap"] for r in S)
            ae = np.concatenate([r["app_e"] for r in S]); ac = np.concatenate([r["app_c"] for r in S])
            f = lambda a, c, d: (a[:, c] < d).mean() if len(a) else np.nan
            say(f"{'elongated' if el else 'circle-like':11s} R_eq {R}: escape {E / wt:6.2f}  annihilation {A / wt:6.2f}  untracked {Lo / wt:6.2f}  | swap {Sw / max(E, 1):.0%} of escapes | "
                f"bulk +1/2 within 3/5 cells: escapes {f(ae, 0, 3):.2f}/{f(ae, 0, 5):.2f} vs control {f(ac, 0, 3):.2f}/{f(ac, 0, 5):.2f} | bulk -1/2 within 3/5: escapes {f(ae, 1, 3):.2f}/{f(ae, 1, 5):.2f} vs control {f(ac, 1, 3):.2f}/{f(ac, 1, 5):.2f}")
    open(os.path.join(H, "channels.txt"), "w").write("\n".join(L) + "\n")
