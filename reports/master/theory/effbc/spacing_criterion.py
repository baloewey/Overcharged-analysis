#!/usr/bin/env python3
"""Spacing criterion for the +1 bulk (user's question, 2026-10-07), tested on the main circle sweep (act 0.07, +k, k >= 2, seed means).
Behind a wall layer of N = 2k - 2 +1/2 with spacing l = 2 pi R / N, the bulk sees an effective planar (+tilt/2 spiral) condition whose
patchiness decays as exp(-2 pi delta / l); it is smooth (< 0.1 rad) beyond ~0.34 l from the wall (effective_bc.py). Candidate controls:
  R / l_a (size), l / l_a (spacing; in this sweep l ~ 1/rho at any k), R_eff / l_a with R_eff = R - 0.34 l (radius of the region that
  sees a smooth effective condition), and R / l (number of spacings per radius ~ N / 2 pi).
Outcomes per (rho, k): bulk charge q_slow (mobility criterion), q_dens (density criterion), pinned lattice f_static (wall Fourier).
For each binary outcome (bulk +1: 0.75 <= q_slow <= 1.25; pinned: f_static >= 0.75) and each single control: best threshold accuracy
and ROC AUC. Writes spacing_criterion.txt and spacing_criterion.png."""
import csv, os
from collections import defaultdict
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
H = os.path.dirname(os.path.abspath(__file__)); CAV = "/home/baloewe/Physics/Overcharged/runs/cavity"; LA = 22.21
D = defaultdict(lambda: defaultdict(list))
for fn, cols in ((CAV + "/analysis_weekend/criteria_all.csv", ("q_slow", "q_dens")), (CAV + "/analysis_fourier/wall_fourier_all.csv", ("fstatic",))):
    for r in csv.DictReader(open(fn)):
        if r["set"] != "pos" or not r["run"].startswith("main_act0.07"): continue
        key = (round(float(r["rho"]), 3), int(r["k"])); D[key]["R"] = [float(r["R"])]
        for c in cols:
            try:
                v = float(r[c])
                if np.isfinite(v): D[key][c].append(v)
            except ValueError: pass
rows = []
for (rho, k), d in sorted(D.items()):
    if k < 2 or not all(d.get(c) for c in ("q_slow", "q_dens", "fstatic")): continue
    R = d["R"][0]; N = 2 * k - 2; ell = 2 * np.pi * R / N
    rows.append(dict(rho=rho, k=k, R=R / LA, ell=ell / LA, Reff=(R - 0.34 * ell) / LA, RoverL=R / ell,
                     q=np.mean(d["q_slow"]), qd=np.mean(d["q_dens"]), f=np.mean(d["fstatic"])))
L = []
def say(s): print(s); L.append(s)
say(f"{len(rows)} (rho, k) points, k = 2..{max(r['k'] for r in rows)}, rho {min(r['rho'] for r in rows)}..{max(r['rho'] for r in rows)}")
def auc(x, y):
    pos, neg = x[y], x[~y]
    if len(pos) == 0 or len(neg) == 0: return np.nan
    return np.mean([(p > n) + 0.5 * (p == n) for p in pos for n in neg])
def best_thr(x, y):
    best = (0, None, None)
    for t in np.unique(x):
        for sgn in (1, -1):
            acc = np.mean((sgn * (x - t) >= 0) == y)
            if acc > best[0]: best = (acc, t, sgn)
    return best
outcomes = {"bulk +1 (0.75 <= q_slow <= 1.25)": np.array([0.75 <= r["q"] <= 1.25 for r in rows]),
            "any bulk charge (q_dens >= 0.5)": np.array([r["qd"] >= 0.5 for r in rows]),
            "pinned lattice (f_static >= 0.75)": np.array([r["f"] >= 0.75 for r in rows])}
ctrl = {"R/l_a": "R", "l/l_a": "ell", "R_eff/l_a": "Reff", "R/l": "RoverL"}
for oname, y in outcomes.items():
    say(f"\n{oname}: {y.sum()} of {len(y)} points")
    for cname, key in ctrl.items():
        x = np.array([r[key] for r in rows]); acc, t, sgn = best_thr(x, y); a = auc(x, y)
        say(f"   {cname:10s}: AUC {max(a, 1 - a):.3f}; best single threshold {'>=' if sgn > 0 else '<='} {t:.2f} -> accuracy {acc:.3f}")
    # two-variable: size and spacing together (axis-aligned box rule R >= t1 and l <= t2)
    xR = np.array([r["R"] for r in rows]); xl = np.array([r["ell"] for r in rows]); best = (0,)
    for t1 in np.unique(xR):
        for t2 in np.unique(xl):
            acc = np.mean(((xR >= t1) & (xl <= t2)) == y)
            if acc > best[0]: best = (acc, t1, t2)
    say(f"   two-sided rule R/l_a >= {best[1]:.2f} and l/l_a <= {best[2]:.2f}: accuracy {best[0]:.3f}")
fig, ax = plt.subplots(1, 2, figsize=(10.5, 4))
xR = np.array([r["R"] for r in rows]); xl = np.array([r["ell"] for r in rows]); q = np.array([r["q"] for r in rows]); f = np.array([r["f"] for r in rows])
sc = ax[0].scatter(xR, xl, c=np.clip(q, 0, 2), cmap="coolwarm", vmin=0, vmax=2, s=30, edgecolor="k", lw=0.3)
plt.colorbar(sc, ax=ax[0], label="bulk charge $q_\\mathrm{slow}$"); ax[0].set_xscale("log"); ax[0].set_yscale("log")
ax[0].set_xlabel(r"$R/\ell_a$"); ax[0].set_ylabel(r"wall spacing $\ell/\ell_a$"); ax[0].set_title("(a) bulk charge", fontsize=10)
sc = ax[1].scatter(xR, xl, c=f, cmap="viridis", vmin=0, vmax=1, s=30, edgecolor="k", lw=0.3)
plt.colorbar(sc, ax=ax[1], label="$f_\\mathrm{static}$ (pinned lattice)"); ax[1].set_xscale("log"); ax[1].set_yscale("log")
ax[1].set_xlabel(r"$R/\ell_a$"); ax[1].set_ylabel(r"$\ell/\ell_a$"); ax[1].set_title("(b) pinned wall lattice", fontsize=10)
for a in ax: a.axhline(1, color="k", lw=0.6, ls="--")
fig.tight_layout(); fig.savefig(os.path.join(H, "spacing_criterion.png"), dpi=150)
open(os.path.join(H, "spacing_criterion.txt"), "w").write("\n".join(L) + "\n")
