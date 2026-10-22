#!/usr/bin/env python3
"""Who leaves a crowded layer? (user list item 2, 2026-10-22). Shape-sweep events (shapes/escape_events/events.npz), lagged outcome (in the wall at +40 t.u.,
deeper than 6 cells at +60 t.u.), so the start frame is >= 40 t.u. before departure. Within each run, each variable is converted to its within-run percentile;
future escapers vs all other starts: mean percentile (50 = no difference) and the escape probability in the bottom vs top quartile. Variables: smaller gap
(min of the two arc gaps), larger gap, gap asymmetry |gl - gr|/(gl + gr), depth, |curvature| x R_eq, |tilt|. Crowded (rho >= 0.12) vs sparse (rho <= 0.06),
aspect >= 1.4 and circles. -> who_leaves.txt"""
import csv, numpy as np
S = "/home/baloewe/Physics/Overcharged/reports/master/theory/shapes"; Z = np.load(S + "/escape_events/" + __import__("os").environ.get("EVFILE", "events.npz"))["E"]
run = Z[:, 0].astype(int); gl, gr = Z[:, 3], Z[:, 4]; late = Z[:, 9] > 0.5
V = {"smaller gap": np.minimum(gl, gr), "larger gap": np.maximum(gl, gr), "gap asymmetry": np.abs(gl - gr) / (gl + gr), "depth": Z[:, 5], "|curvature| R_eq": Z[:, 6], "|tilt|": np.abs(Z[:, 2])}
meta = list(csv.DictReader(open(S + "/escape_events/runs.csv"))); SR = {r["run"]: r for r in csv.DictReader(open(S + "/shape_runs.csv"))}
rho = np.array([float(SR[m["run"]]["rho"]) for m in meta]); asp = np.array([float(SR[m["run"]]["aspect"]) for m in meta])
order = np.argsort(run, kind="stable"); b = np.r_[0, np.flatnonzero(np.diff(run[order])) + 1, len(run)]
PCT = {}
for nm, x in V.items():
    pc = np.empty(len(x))
    for a, c in zip(b[:-1], b[1:]):
        ii = order[a:c]; pc[ii] = (np.argsort(np.argsort(x[ii] + 1e-9 * np.random.default_rng(0).random(len(ii)))) + 0.5) / len(ii) * 100
    PCT[nm] = pc
L = ["WHO LEAVES (lagged escapes, start >= 40 t.u. before departure; within-run percentiles; 50 = no difference)"]
for cls, sel_r in (("aspect>=1.4, crowded rho>=0.12", (asp >= 1.4) & (rho >= 0.11)), ("aspect>=1.4, sparse rho<=0.06", (asp >= 1.4) & (rho <= 0.065)),
                   ("circles, crowded", (asp < 1.05) & (rho >= 0.11)), ("circles, sparse", (asp < 1.05) & (rho <= 0.065))):
    m = np.isin(run, np.flatnonzero(sel_r)); e = m & late; L.append(f"\n{cls}: {m.sum()} starts, {e.sum()} lagged escapes (p {1000 * e.sum() / m.sum():.2f} per 1000)")
    for nm in V:
        pc = PCT[nm]; lo = m & (pc < 25); hi = m & (pc > 75)
        L.append(f"  {nm:16s}: escapers' mean percentile {pc[e].mean():5.1f} (+-{pc[e].std() / np.sqrt(max(e.sum(), 1)):.1f}) | p esc bottom quartile {1000 * late[lo].mean():6.2f}, top {1000 * late[hi].mean():6.2f} per 1000")
open("who_leaves.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
