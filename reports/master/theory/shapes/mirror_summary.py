#!/usr/bin/env python3
"""Mirror control summary (2026-10-19): +k shape sweep (shape_runs.csv, wall +1/2) vs -k mirror sweep (shape_runs_mirror.csv, wall -1/2), same cavities.
Per (family, aspect, R_eq): medians over k, rho, seeds of
  +k: vacancies 2k-2-N, bulk charge Q_b, escape prob per wall defect per 20 t.u., |circulation|
  -k: vacancies 2|k|+2-N, gain N-2|k|, bulk charge, escape prob, |circulation|, opposite-sign (+1/2) wall defects.
Delta = value minus the circle at the same R_eq (does elongation empty the layer?). Writes mirror_summary.txt."""
import csv, numpy as np
from collections import defaultdict
def load(f):
    G = defaultdict(lambda: defaultdict(list))
    for r in csv.DictReader(open(f)):
        if r.get("nwall") in (None, ""): continue
        key = (r["family"], round(float(r["aspect"]), 2), r["Req_la"][:3])
        for c in ("vac", "gain", "qbulk", "esc", "circ", "wall_opp", "two"):
            if r.get(c) not in (None, ""): G[key][c].append(float(r[c]))
    return G
P, M = load("shape_runs.csv"), load("shape_runs_mirror.csv")
med = lambda G, k, c: np.nanmedian(G[k][c]) if G[k][c] else np.nan
out = ["Mirror control: +k (wall +1/2) vs -k (wall -1/2), same 675 cavities each (3 seeds; k, rho pooled; medians). Delta = minus the circle at same R_eq.",
       f"{'family':8s} {'asp':>5s} {'Req':>4s} | {'+k vac':>6s} {'dvac':>5s} {'Qb':>5s} {'esc':>7s} {'circ':>5s} | {'-k vac':>6s} {'dvac':>5s} {'gain':>5s} {'Qb':>5s} {'esc':>7s} {'circ':>5s} {'opp':>5s} {'n':>3s}"]
for key in sorted(M, key=lambda k: (k[2], k[0] != "circle", k[0], k[1])):
    c = ("circle", 1.0, key[2])
    out.append(f"{key[0]:8s} {key[1]:5.2f} {key[2]:>4s} | {med(P,key,'vac'):6.2f} {med(P,key,'vac')-med(P,c,'vac'):5.2f} {med(P,key,'qbulk'):5.2f} {med(P,key,'esc'):7.4f} {med(P,key,'circ'):5.2f} | "
               f"{med(M,key,'vac'):6.2f} {med(M,key,'vac')-med(M,c,'vac'):5.2f} {med(M,key,'gain'):5.2f} {med(M,key,'qbulk'):5.2f} {med(M,key,'esc'):7.4f} {med(M,key,'circ'):5.2f} {med(M,key,'wall_opp'):5.2f} {len(M[key]['vac']):3d}")
# pooled by aspect class
out.append("\nPooled (all R_eq): aspect class | +k dvac, esc | -k dvac, esc")
for lo, hi, nm in ((0.99, 1.01, "circle"), (1.1, 1.3, "1.1-1.29"), (1.39, 9, ">=1.4")):
    def pool(G, c):
        v = []
        for k in G:
            if lo <= k[1] <= hi: v += [x - med(G, ("circle", 1.0, k[2]), c) if c == "vac" else x for x in G[k][c]]
        return np.nanmedian(v), len(v)
    (pv, n1), (pe, _), (mv, n2), (me, _) = pool(P, "vac"), pool(P, "esc"), pool(M, "vac"), pool(M, "esc")
    out.append(f"  {nm:9s} +k dvac {pv:5.2f} esc {pe:.4f} (n {n1}) | -k dvac {mv:5.2f} esc {me:.4f} (n {n2})")
open("mirror_summary.txt", "w").write("\n".join(out) + "\n"); print("\n".join(out))
