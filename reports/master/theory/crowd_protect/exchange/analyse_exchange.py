#!/usr/bin/env python3
"""Clean exchange hazards by wall density (exchange_count.csv). Escape hazard = escapes / wall-resident exposure; entry hazard = entries / bulk-resident
exposure (per defect per 1000 t.u.); wall-born rate = tracks born in the wall per wall length per 1000 t.u. (pair creation + tracker gaps); bulk charge per length.
Medians over runs. -> analyse_exchange.txt"""
import csv, numpy as np
R = [r for r in csv.DictReader(open("exchange_count.csv")) if not r.get("err")]; S = {r["run"]: r for r in csv.DictReader(open("../../shapes/shape_runs.csv"))}
f = lambda r, k: float(r[k]); LV = np.array([0.04, 0.06, 0.09, 0.12, 0.14])
def row(rs):
    he = [1000 * f(r, "esc") / max(f(r, "expW"), 1) for r in rs]; hi = [1000 * f(r, "ent") / max(f(r, "expB"), 1) for r in rs]
    bw = [1000 * f(r, "born_w") / f(r, "Lw") / f(r, "T") for r in rs]; ef = [1000 * f(r, "esc") / f(r, "Lw") / f(r, "T") for r in rs]; nf = [1000 * f(r, "ent") / f(r, "Lw") / f(r, "T") for r in rs]
    qb = [f(r, "qb") / f(r, "Lw") for r in rs]; nb = [f(r, "nB") / f(r, "Lw") for r in rs]
    return (f"{len(rs):3d} | {np.median(he):6.2f} | {np.median(hi):6.2f} | {np.median(ef):.4f} | {np.median(nf):.4f} | {np.median(bw):.4f} | {np.median(nb):.4f} | {np.median(qb):.4f}")
L = ["CLEAN EXCHANGE (tracked; escape = >= 20 t.u. in wall -> bulk, no return within 20 t.u.; entry = >= 20 t.u. in bulk -> wall, stays 20 t.u.)",
     "class | rho | n | escape hazard per wall defect /1000 t.u. | entry hazard per bulk +1/2 /1000 t.u. | escapes per length /1000 t.u. | entries per length | wall-born tracks per length | bulk +1/2 per length | bulk charge per length"]
for nm, lo, hi in (("circles", 0.99, 1.01), ("aspect>=1.4", 1.39, 9)):
    for i, r0 in enumerate(LV):
        rs = [r for r in R if r["kind"] == "cavity" and lo <= float(S[r["run"].split("shapes_server/")[1]]["aspect"]) <= hi and np.argmin(abs(f(r, "rho") - LV)) == i]
        L.append(f"{nm:11s} | {r0:.2f} | " + row(rs))
ch = [r for r in R if r["kind"] == "channel"]
for m in sorted(set(int(r["m"]) for r in ch)):
    rs = [r for r in ch if int(r["m"]) == m]; L.append(f"channel m{m:<3d}| {2 * m / 192:.3f} | " + row(rs) + "   " + ", ".join(sorted(set(r["run"].split("/")[0] for r in rs))))
open("analyse_exchange.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
