#!/usr/bin/env python3
"""3-seed vs 10-seed comparison of the main +k shape-sweep results (27 Oct). Medians over runs per class and wall-density level. -> compare_3v10.txt"""
import csv, numpy as np
LV = np.array([0.04, 0.06, 0.09, 0.12, 0.14]); L = ["3-SEED vs 10-SEED (+k shape sweep; medians over runs; '3s' / '10s')"]
def load(fn): return [r for r in csv.DictReader(open(fn)) if not r.get("err")]
S10 = {r["run"]: r for r in load("shapes/shape_runs.csv")}; S3 = {r["run"]: r for r in load("shapes/shape_runs_3seeds.csv")}
E10 = {r["run"]: r for r in load("crowd_protect/exchange/shape_esc2.csv")}; E3 = {r["run"]: r for r in load("crowd_protect/exchange/shape_esc2_3seeds.csv")}
lv = lambda r: int(np.argmin(abs(float(r["rho"]) - LV)))
def table(title, fn):
    L.append(f"\n{title}"); L.append("class          | " + " | ".join(f"rho {x:.2f}        " for x in LV))
    for nm, lo, hi in (("circles", 0.99, 1.01), ("aspect 1.1-1.3", 1.05, 1.35), ("aspect >= 1.4", 1.39, 9)):
        row = []
        for i in range(5):
            v3 = [fn(S3, E3, k) for k, r in S3.items() if lo <= float(r["aspect"]) <= hi and lv(r) == i]; v10 = [fn(S10, E10, k) for k, r in S10.items() if lo <= float(r["aspect"]) <= hi and lv(r) == i]
            v3 = [x for x in v3 if np.isfinite(x)]; v10 = [x for x in v10 if np.isfinite(x)]
            row.append(f"{np.median(v3):.4f} / {np.median(v10):.4f}")
        L.append(f"{nm:14s} | " + " | ".join(row))
f = lambda D, k, c: float(D[k][c]) if D[k][c] not in ("", "nan") else np.nan
table("Escape probability per wall +1/2 per 20 t.u. (validated tracker): 3s / 10s", lambda S, E, k: float(E[k]["p_esc"]) if k in E else np.nan)
table("Vacant fraction (2k-2-N_wall)/(2k-2): 3s / 10s", lambda S, E, k: f(S, k, "vac") / (2 * int(S[k]["k"]) - 2))
table("Tilt coherence C: 3s / 10s", lambda S, E, k: f(S, k, "coh"))
table("Bulk charge Q_b: 3s / 10s", lambda S, E, k: f(S, k, "qbulk"))
table("|circulation|: 3s / 10s", lambda S, E, k: f(S, k, "circ"))
L.append("\nCoherence gating (aspect >= 1.1): escape x1000 coherent (C > 0.85) / incoherent, 3s -> 10s")
for i, r0 in enumerate(LV):
    out = []
    for S, E in ((S3, E3), (S10, E10)):
        g = [k for k, r in S.items() if float(r["aspect"]) > 1.05 and lv(r) == i]
        c = [1000 * float(E[k]["p_esc"]) for k in g if float(S[k]["coh"]) > 0.85]; ic = [1000 * float(E[k]["p_esc"]) for k in g if float(S[k]["coh"]) < 0.85]
        out.append(f"{np.median(c):6.2f} / {np.median(ic):6.2f} (n {len(c)}/{len(ic)})")
    L.append(f"  rho {r0:.2f}: " + "  ->  ".join(out))
X3 = {r["run"]: r for r in load("crowd_protect/exchange/exchange_count2_g9_3seeds.csv")}; X10 = {r["run"]: r for r in load("crowd_protect/exchange/exchange_count2_g9.csv")}
L.append("\nClean exchange, aspect >= 1.4: escape hazard / entry hazard per 1000 t.u., 3s -> 10s")
for i, r0 in enumerate(LV):
    out = []
    for X, S in ((X3, S3), (X10, S10)):
        g = [k for k, r in X.items() if r["kind"] == "cavity" and k.split("shapes_server/")[1] in S and float(S[k.split("shapes_server/")[1]]["aspect"]) >= 1.39 and lv(S[k.split("shapes_server/")[1]]) == i]
        out.append(f"{np.median([1000 * float(X[k]['esc']) / max(float(X[k]['expW']), 1) for k in g]):5.2f} / {np.median([1000 * float(X[k]['ent']) / max(float(X[k]['expB']), 1) for k in g]):4.2f} (n {len(g)})")
    L.append(f"  rho {r0:.2f}: " + "  ->  ".join(out))
open("compare_3v10.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
