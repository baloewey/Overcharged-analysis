#!/usr/bin/env python3
"""Old (exchange_count.py: 3.5-cell mutual nearest neighbour, gaps 2) vs predictive tracker (exchange_count2.py, gates 6 and 9) for the clean exchange
hazards, by wall density; plus the run-level escape probability of shape_runs.csv (10 t.u. nearest-neighbour chain). Medians over runs.
-> compare_trackers.txt"""
import csv, numpy as np
S = {r["run"]: r for r in csv.DictReader(open("../../shapes/shape_runs.csv"))}; LV = np.array([0.04, 0.06, 0.09, 0.12, 0.14])
def load(fn): return {r["run"]: r for r in csv.DictReader(open(fn)) if not r.get("err")}
O, G6, G9 = load("exchange_count.csv"), load("exchange_count2_g6.csv"), load("exchange_count2_g9.csv")
h = lambda r, a, b: 1000 * float(r[a]) / max(float(r[b]), 1)
L = ["ESCAPE / ENTRY HAZARDS (per defect per 1000 t.u.), old tracker vs predictive (gate 6 / gate 9); medians over runs; and shape_runs p_esc (per 20 t.u., x1000)",
     "class | rho | n | escape: old, g6, g9 | entry: old, g6, g9 | shape_runs p_esc x1000 | median track length g9 (t.u.) | exact check: slope of charge beyond 4 cells (per length per 1000 t.u.)"]
for nm, lo, hi in (("circles", 0.99, 1.01), ("aspect>=1.4", 1.39, 9)):
    for i, r0 in enumerate(LV):
        rs = [k for k in G9 if k in O and k in G6 and G9[k]["kind"] == "cavity" and lo <= float(S[k.split("shapes_server/")[1]]["aspect"]) <= hi and np.argmin(abs(float(G9[k]["rho"]) - LV)) == i]
        md = lambda D, a, b: np.median([h(D[k], a, b) for k in rs])
        L.append(f"{nm:11s} | {r0:.2f} | {len(rs):3d} | {md(O, 'esc', 'expW'):6.2f}, {md(G6, 'esc', 'expW'):6.2f}, {md(G9, 'esc', 'expW'):6.2f} | {md(O, 'ent', 'expB'):5.2f}, {md(G6, 'ent', 'expB'):5.2f}, {md(G9, 'ent', 'expB'):5.2f} | "
                 f"{1000 * np.median([float(S[k.split('shapes_server/')[1]]['esc']) for k in rs]):6.2f} | {np.median([float(G9[k]['med_len']) for k in rs]):5.0f} | {np.median([float(G9[k]['qb_slope']) for k in rs]):+.5f}")
ch = sorted({int(G9[k]["m"]) for k in G9 if G9[k]["kind"] == "channel"})
for m in ch:
    rs = [k for k in G9 if G9[k]["kind"] == "channel" and int(G9[k]["m"]) == m and k in O and k in G6]
    md = lambda D, a, b: np.median([h(D[k], a, b) for k in rs])
    L.append(f"channel m{m:<3d}| {2 * m / 192:.3f} | {len(rs):3d} | {md(O, 'esc', 'expW'):6.2f}, {md(G6, 'esc', 'expW'):6.2f}, {md(G9, 'esc', 'expW'):6.2f} | {md(O, 'ent', 'expB'):5.2f}, {md(G6, 'ent', 'expB'):5.2f}, {md(G9, 'ent', 'expB'):5.2f} | - | {np.median([float(G9[k]['med_len']) for k in rs]):5.0f} | {np.median([float(G9[k]['qb_slope']) for k in rs]):+.5f}")
open("compare_trackers.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
