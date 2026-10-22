#!/usr/bin/env python3
"""Recheck of the PRL exit numbers (0.068 per wall defect per 20 t.u. at 0 vacancies -> 0.028 at 5; provenance/dsqrts_exit.txt, perdef_ff36.npy, k 6, N0 10).
(1) perdef_ff36: within-run rate by frame vacancies (each run's own rates, averaged over runs with weights = frames; between-run differences cancel).
(2) shape sweep (events.npz): prompt escape (20 t.u.) by frame vacancies (2k-2 - wall +1/2 count), within run, at each rho level (fixed density)."""
import csv, numpy as np
T = "/home/baloewe/Physics/Overcharged/reports/master/theory"
E = [e for e in np.load(f"{T}/energetics/perdef_ff36.npy", allow_pickle=True) if "err" not in e and e["eps"] >= 1.3 and len(e["rows"])]
def within(groups, V):
    """groups: list of (vac array, outcome array) per run. Pooled rate and within-run rate (run's rate at v minus run's mean, plus grand mean), per v."""
    allv = np.concatenate([g[0] for g in groups]); ally = np.concatenate([g[1] for g in groups]); gm = ally.mean(); out = []
    for v in V:
        pooled = ally[allv == v]; num = den = 0.0; nr = 0
        for vv, y in groups:
            m = vv == v
            if m.sum() >= 20: num += (y[m].mean() - y.mean()) * m.sum(); den += m.sum(); nr += 1
        out.append((v, pooled.mean() if len(pooled) else np.nan, len(pooled), gm + num / den if den else np.nan, nr))
    return out
L = ["(1) perdef_ff36 (aspect 1.4-2, lambda 1.2/2/4, k 6): exit rate per wall defect per 20 t.u. by frame vacancies; pooled vs within-run (run mean removed)"]
G = [((e["N0"] - e["rows"][:, 2]).astype(int), e["rows"][:, 6].astype(float)) for e in E]
for v, p, n, w, nr in within(G, range(0, 7)):
    if n > 500: L.append(f"  vac {v}: pooled {p:.4f} | within-run {w:.4f} (n {n}, runs {nr})")
Z = np.load("events.npz"); X = Z["E"]; run = X[:, 0].astype(int); esc = X[:, 7]
meta = list(csv.DictReader(open("runs.csv"))); SR = {r["run"]: r for r in csv.DictReader(open("../shape_runs.csv"))}
k = np.array([int(SR[m["run"]]["k"]) for m in meta]); rho = np.array([float(SR[m["run"]]["rho"]) for m in meta]); asp = np.array([float(SR[m["run"]]["aspect"]) for m in meta])
lev = np.argmin(np.abs(rho[:, None] - np.array([0.04, 0.06, 0.09, 0.12, 0.14])), 1); vac = (2 * k[run] - 2) - X[:, 8].astype(int)
L.append("\n(2) shape sweep, aspect >= 1.4: prompt escape (20 t.u.) per wall +1/2 by frame vacancies, within-run, at fixed rho level (x1000; runs contributing)")
for i, r0 in enumerate((0.04, 0.06, 0.09, 0.12, 0.14)):
    rs = np.flatnonzero((lev == i) & (asp >= 1.4)); G = [(vac[run == j], esc[run == j]) for j in rs]
    L.append(f"  rho {r0:.2f} ({len(rs)} runs): " + "; ".join(f"vac {v}: {1000 * w:.1f} (pooled {1000 * p:.1f}, {nr})" for v, p, n, w, nr in within(G, range(-1, 9)) if n > 2000 and nr >= 5))
open("vacancy_check.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
