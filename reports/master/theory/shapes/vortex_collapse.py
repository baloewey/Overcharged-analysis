#!/usr/bin/env python3
"""Which length sets the loss of the single vortex in elongated cavities? (user, 2026-10-20). Shape sweep (shape_runs.csv: medians over k, rho, seeds per family,
aspect, R_eq). Geometry from the actual wall polyline of one run per shape: overall length Lx and width Wy (cells), in units of l_a = 22.21. Candidate variables:
  aspect = Lx/Wy (scale-free: width-set vortices),  extra length (Lx - Wy)/l_a (an active-length criterion),  Lx/l_a.
Collapse test: for each variable, the crossing of the two-vortex fraction through 0.25 and of the tilt coherence through 0.7 at each size (R_eq 1, 1.5, 2) and family;
spread of the crossings (coefficient of variation) -- the variable whose crossings agree across sizes sets the transition. -> vortex_collapse.txt, fig_vortex_collapse.{pdf,png}"""
import csv, glob, json, sys, numpy as np
from collections import defaultdict
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, "/home/baloewe/Physics/Overcharged/runs/cavity"); from wall_fourier import wall_polyline  # noqa
LA = 22.21; SW = "/home/baloewe/Physics/Overcharged/runs/shapes_server"
G = defaultdict(lambda: defaultdict(list)); geo = {}
for r in csv.DictReader(open("shape_runs.csv")):
    key = (r["family"], round(float(r["aspect"]), 2), r["Req_la"][:3]); G[key]["coh"].append(float(r["coh"])); G[key]["two"].append(float(r["two"]) if r["two"] else np.nan)
    if key not in geo:
        bc = json.load(open(f"{SW}/{r['run']}/input.json"))["BC"][0]; W, s, P = wall_polyline(bc); geo[key] = (np.ptp(W[:, 0]), np.ptp(W[:, 1]))
rows = []
for key, v in G.items():
    Lx, Wy = geo[key]; rows.append(dict(fam=key[0], req=key[2], asp=Lx / Wy, extra=(Lx - Wy) / LA, Lx=Lx / LA, coh=np.nanmedian(v["coh"]), two=np.nanmedian(v["two"])))
VARS = [("asp", "Aspect ratio $L_x/W_y$ (width-set vortices)"), ("extra", r"Extra length $(L_x-W_y)/\ell_a$"), ("Lx", r"Length $L_x/\ell_a$")]
def crossing(xs, ys, thr, up):
    o = np.argsort(xs); xs, ys = np.array(xs)[o], np.array(ys)[o]
    for i in range(len(xs) - 1):
        a, b = ys[i] - thr, ys[i + 1] - thr
        if (a < 0 <= b) if up else (a >= 0 > b): return xs[i] + (xs[i + 1] - xs[i]) * abs(a) / (abs(a) + abs(b))
    return np.nan
L = ["Collapse test: crossings of two-vortex fraction (0.25, rising) and tilt coherence (0.7, falling) per family and size, in each candidate variable; CV = sd/mean across the 6 (family, size) curves."]
cols = {"1.0": "#9ecae1", "1.5": "#4a90c2", "2.0": "#08306b"}
fig, axs = plt.subplots(2, 3, figsize=(15, 8.2))
for j, (v, lab) in enumerate(VARS):
    for qi, (q, thr, up, ylab) in enumerate((("coh", 0.7, False, "Tilt coherence $C$"), ("two", 0.25, True, "Fraction of time with two vortices"))):
        ax = axs[qi, j]; cr = []
        for fam, mk, ls in (("ellipse", "o", "-"), ("oval", "s", "--")):
            for req in cols:
                S = sorted([r for r in rows if (r["fam"] == fam or r["fam"] == "circle") and r["req"] == req], key=lambda r: r[v])
                xs, ys = [r[v] for r in S], [r[q] for r in S]; ax.plot(xs, ys, ls=ls, marker=mk, ms=4, color=cols[req], lw=1.1); cr.append(crossing(xs, ys, thr, up))
        cr = np.array(cr); cv = np.nanstd(cr) / np.nanmean(cr) if np.isfinite(cr).sum() > 1 else np.nan
        ax.axhline(thr, color="0.7", lw=0.7, ls=":"); ax.set_xlabel(lab); ax.set_ylabel(ylab)
        ax.set_title(f"Crossings {np.nanmin(cr):.2f}–{np.nanmax(cr):.2f} (CV {cv:.2f})", fontsize=10)
        L.append(f"{v:6s} {q:4s}: crossings {np.round(cr, 2)} (order: ellipse R_eq 1/1.5/2, oval R_eq 1/1.5/2); range {np.nanmin(cr):.2f}-{np.nanmax(cr):.2f}; CV {cv:.3f}")
        if v == "asp": ax.set_xlim(0.95, 2.6)
from matplotlib.lines import Line2D
h = [Line2D([], [], color=cols[r], lw=1.5, label=fr"$R_\mathrm{{eq}}={r}\,\ell_a$") for r in cols] + [Line2D([], [], color="0.3", marker="o", label="Ellipse"), Line2D([], [], color="0.3", marker="s", ls="--", label="Cassini oval")]
fig.legend(handles=h, loc="upper center", ncol=5, frameon=False, fontsize=10); fig.subplots_adjust(top=0.9, hspace=0.38, wspace=0.3, left=0.06, right=0.98, bottom=0.08)
for i, a in enumerate(axs.flat): a.text(-0.12, 1.05, "(" + "abcdef"[i] + ")", transform=a.transAxes, fontweight="bold")
fig.savefig("fig_vortex_collapse.pdf"); fig.savefig("fig_vortex_collapse.png", dpi=150)
open("vortex_collapse.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
