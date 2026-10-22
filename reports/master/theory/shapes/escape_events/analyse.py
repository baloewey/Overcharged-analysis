#!/usr/bin/env python3
"""Per-event escape analysis (events.npz from escape_events.py). Local crowding = mean arc gap to the two wall +1/2 neighbours (cells); local tilt = |defect
orientation - outward normal| (deg) at the start frame. (1) pooled escape probability by tilt x gap; (2) within-run Mantel-Haenszel odds ratios (each run split at
its own median), so run-level differences (shape, size, k, rho) cancel; (3) the same with the other variable held in within-run terciles."""
import csv, numpy as np
Z = np.load(__import__("os").environ.get("EVFILE", "events.npz")); E = Z["E"]; run = E[:, 0].astype(int); tilt = np.abs(E[:, 2]); gap = 0.5 * (E[:, 3] + E[:, 4]); import sys
OUT = sys.argv[1] if len(sys.argv) > 1 else "esc"
esc = (E[:, 7] if OUT == "esc" else E[:, 9]) > 0.5
if OUT == "late": print("OUTCOME: lagged escape (in the wall at +40 t.u., depth > 6 at +60 t.u.); tilt and gap measured >= 40 t.u. before departure")
meta = {int(r["idx"]): r for r in csv.DictReader(open("runs.csv"))}
asp = np.array([float(meta[i]["aspect"]) for i in range(len(meta))])[run]
print(f"{len(E)} wall +1/2 starts, {esc.sum()} escapes (p {esc.mean():.4f}); 20 t.u. window; {len(meta)} runs")
def table(sel, title):
    print(f"\n{title}: escape probability x1000 (n starts / 1000)")
    gb = [0, 6, 9, 12, 16, 22, 35, 1e9]; tb = [0, 30, 45, 60, 75, 90, 180]
    print("gap (cells)  | " + " | ".join(f"tilt {a:>3d}-{b:<3d}" for a, b in zip(tb[:-1], tb[1:])))
    for g0, g1 in zip(gb[:-1], gb[1:]):
        row = []
        for t0, t1 in zip(tb[:-1], tb[1:]):
            m = sel & (gap >= g0) & (gap < g1) & (tilt >= t0) & (tilt < t1)
            row.append(f"{1000 * esc[m].mean():5.1f} ({m.sum() / 1000:5.1f})" if m.sum() >= 2000 else "      -       ")
        print(f"{g0:4.0f}-{g1 if g1 < 1e8 else 999:4.0f}    | " + " | ".join(row))
table(np.ones(len(E), bool), "All shapes")
table(asp >= 1.4, "Elongated (aspect >= 1.4)")
# within-run MH odds ratios
order = np.argsort(run, kind="stable"); bounds = np.r_[0, np.flatnonzero(np.diff(run[order])) + 1, len(run)]
def mh(xvar, strat=None):
    num = den = 0.0
    for a, b in zip(bounds[:-1], bounds[1:]):
        ii = order[a:b]
        groups = [ii] if strat is None else [ii[q] for q in np.array_split(np.argsort(strat[ii]), 3)]
        for g in groups:
            if len(g) < 40: continue
            hi = xvar[g] > np.median(xvar[g]); e = esc[g]; n = len(g)
            A_, B_ = (e & hi).sum(), (~e & hi).sum(); C_, D_ = (e & ~hi).sum(), (~e & ~hi).sum()
            num += A_ * D_ / n; den += B_ * C_ / n
    return num / den
print("\nWithin-run Mantel-Haenszel odds ratio, escape | upper vs lower half of the variable (1 = no effect)")
print(f"  wider gap (sparser locally):            {mh(gap):.2f}    | holding tilt in terciles: {mh(gap, tilt):.2f}")
print(f"  larger |tilt| (more tangential):        {mh(tilt):.2f}    | holding gap in terciles:  {mh(tilt, gap):.2f}")
dep = E[:, 5]; print(f"  deeper start:                           {mh(dep):.2f}")
