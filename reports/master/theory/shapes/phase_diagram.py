#!/usr/bin/env python3
"""Shape phase diagram from shape_runs.csv (seed means per (family, aspect, R_eq, k)). Panels: vacancies, bulk charge, tilt coherence, P(two vortices),
circulation vs aspect for each R_eq (ellipses and ovals); transition aspect per R_eq and family = where the tilt coherence (median over k) drops below 0.75
(linear interpolation); the transition expressed as the long semi-axis A/l_a, short semi-axis B/l_a and perimeter P/l_a (to see which length sets it).
Writes phase_shapes.png and phase_shapes.txt."""
import csv, os, json, glob
import numpy as np
from collections import defaultdict
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
H = os.path.dirname(os.path.abspath(__file__)); SW = "/home/baloewe/Physics/Overcharged/runs/shapes_server"; LA = 22.21
rows = list(csv.DictReader(open(os.path.join(H, "shape_runs.csv"))))
def fl(x):
    try: return float(x)
    except: return np.nan
# geometry per (family, aspect, Req): semi-axes and perimeter from one input.json
geo = {}
for r in rows:
    key = (r["family"], round(fl(r["aspect"]), 2), r["Req_la"][:3])
    if key in geo: continue
    bc = json.load(open(os.path.join(SW, r["run"], "input.json")))["BC"][0]
    import sys; sys.path.insert(0, "/home/baloewe/Physics/Overcharged/runs/cavity"); from wall_fourier import wall_polyline
    W, s, P = wall_polyline(bc); geo[key] = (np.abs(W[:, 0]).max() / LA, np.abs(W[:, 1]).max() / LA, P / LA)
G = defaultdict(lambda: defaultdict(list))
for r in rows:
    key = (r["family"], round(fl(r["aspect"]), 2), r["Req_la"][:3])
    for c in ("vac", "qbulk", "coh", "two", "circ", "esc"): G[key][c].append(fl(r[c]))
L = []
def say(s): print(s); L.append(s)
fig, ax = plt.subplots(1, 5, figsize=(20, 4.6)); cols = {"1.0": "C0", "1.5": "C1", "2.0": "C2"}; mk = {"ellipse": "o-", "oval": "s--", "circle": "*"}
for j, (c, lab) in enumerate((("vac", "vacancies ($2k-2$ minus wall count)"), ("qbulk", "bulk charge"), ("coh", "tilt coherence"), ("two", "P(two vortices)"), ("circ", "|circulation|"))):
    for fam in ("ellipse", "oval", "circle"):
        for req in ("1.0", "1.5", "2.0"):
            ks = sorted(k for k in G if k[0] == fam and k[2] == req)
            if not ks: continue
            xs = [k[1] for k in ks]; ys = [np.nanmedian(G[k][c]) for k in ks]
            ax[j].plot(xs, ys, mk[fam], color=cols[req], ms=4, label=f"{fam} R_eq {req}" if j == 0 else None)
    ax[j].set_xlabel("aspect ratio"); ax[j].set_title(f"({'abcde'[j]}) {lab}", fontsize=9)
fig.legend(*ax[0].get_legend_handles_labels(), fontsize=7, frameon=False, ncol=9, loc="lower center"); fig.tight_layout(rect=(0, 0.1, 1, 1)); fig.savefig(os.path.join(H, "phase_shapes.png"), dpi=150)
say("transition (tilt coherence < 0.75, median over k): aspect, and A/l_a, B/l_a, P/l_a at the transition")
for fam in ("ellipse", "oval"):
    for req in ("1.0", "1.5", "2.0"):
        ks = sorted(k for k in G if k[0] == fam and k[2] == req); xs = np.array([k[1] for k in ks]); ys = np.array([np.nanmedian(G[k]["coh"]) for k in ks])
        circ_coh = np.nanmedian(G[("circle", 1.0, req)]["coh"]) if ("circle", 1.0, req) in G else np.nan
        xs = np.r_[1.0, xs]; ys = np.r_[circ_coh, ys]; o = np.argsort(xs); xs, ys = xs[o], ys[o]
        i = np.argmax(ys < 0.75) if (ys < 0.75).any() else None
        if i:
            a_t = xs[i - 1] + (0.75 - ys[i - 1]) * (xs[i] - xs[i - 1]) / (ys[i] - ys[i - 1])
            k1 = (fam, round(xs[i - 1], 2), req) if xs[i - 1] > 1 else (fam, round(xs[i], 2), req); k2 = (fam, round(xs[i], 2), req)
            g1 = geo.get(k1, geo[k2]); g2 = geo[k2]; w = (a_t - xs[i - 1]) / (xs[i] - xs[i - 1]) if xs[i] > xs[i - 1] else 1
            A_t, B_t, P_t = [g1[m] + w * (g2[m] - g1[m]) for m in range(3)]
            say(f"  {fam:7s} R_eq {req}: transition aspect {a_t:.2f}; A/l_a {A_t:.2f}, B/l_a {B_t:.2f}, P/l_a {P_t:.2f}")
        else: say(f"  {fam:7s} R_eq {req}: no transition in range")
open(os.path.join(H, "phase_shapes.txt"), "w").write("\n".join(L) + "\n")
