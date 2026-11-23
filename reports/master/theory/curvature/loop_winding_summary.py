#!/usr/bin/env python3
"""Summary of loop_winding.csv by class (frame-weighted). -> loop_winding.txt"""
import csv, os, re, numpy as np
H = os.path.dirname(os.path.abspath(__file__)); LA = 22.21
rows = list(csv.DictReader(open(os.path.join(H, "loop_winding.csv"))))
def cls(r):
    p = r["run"]
    if p.startswith("cavity/"): return f"circles (main sweep) R/l_a {('<0.8', '0.8-2.2', '>2.2')[(float(r['R']) / LA >= 0.8) + (float(r['R']) / LA > 2.2)]}"
    m = re.search(r"shapes/(\w+)/([\w.]+)/", p); kind, par = m.group(1), m.group(2)
    if kind == "circle": return "shapes: circle"
    return f"shapes: {kind} {par}"
G = {}
for r in rows: G.setdefault(cls(r), []).append(r)
L = ["DIRECT WINDING TEST (loop_winding.py; director every 100 t.u., second half; loops at depth y inside the wall)",
     "aligned = |alpha| < 75 deg all around the loop (director never near perpendicular to the loop); frame-weighted over runs",
     "class | runs | y | P(w = Q_in) | P(aligned) | P(Q_in = 1 | aligned) | P(Q_in = 1 | not aligned) | <Q_in>"]
for c in sorted(G):
    g = G[c]
    for y in ("3", "5"):
        n = np.array([float(r[f"n{y}"]) for r in g]); f = lambda key: np.array([float(r[f"{key}{y}"]) if r[f"{key}{y}"] not in ("nan", "") else np.nan for r in g])
        pal = f("Pal"); nal = n * pal; nnal = n * (1 - pal)
        q1a = np.nansum(f("Q1al") * nal) / max(np.sum(nal[np.isfinite(f("Q1al"))]), 1); q1n = np.nansum(f("Q1nal") * nnal) / max(np.sum(nnal[np.isfinite(f("Q1nal"))]), 1)
        L.append(f"{c:40s} | {len(g):4d} | {y} | {np.sum(f('Peq') * n) / n.sum():.3f} | {np.sum(pal * n) / n.sum():.2f} | {q1a:.3f} | {q1n:.2f} | {np.sum(f('Qin') * n) / n.sum():+.2f}")
open(os.path.join(H, "loop_winding.txt"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
