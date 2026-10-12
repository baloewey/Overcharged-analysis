#!/usr/bin/env python3
"""
Wall-lattice dynamics in ovals and ellipses vs circles (wall_fourier_dyn.csv + wall_fourier_dyn_circles_k06.csv).
  (a) k = 6: static fraction f_static of the wall lattice vs aspect ratio; circles (aspect 1) of matching size
      in the colour of the family they match (H34 ~ R 1.43 l_a, H48 ~ 2.15, H64 ~ 2.87, ellipses ~ 1.23).
  (b) k = 6: per-defect drift of wall +1/2 along the wall, divided by the near-wall flow |v_w| (< 0 = upstream).
  (c) all runs: drift ratio vs vacancy fraction 1 - N_wall/(2k - 2) of the wall lattice (circles k = 6, all R;
      ovals and ellipses, all k); marker fill = f_static >= 0.75 (pinned) or open. Spearman rho printed.
  (d) k scan at half-length 34: f_static and drift ratio vs k, b/a 2.0 (aspect 1.29) and 1.6 (aspect 1.51).
Writes fig_wall_dynamics.png / .pdf.
"""
import csv, os, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
LA = 22.21
INK, INK2, GRID = "#15171b", "#555a62", "#e9ebee"
FAM = {"H34": ("#2a78d6", "o", "ovals, half-length 1.53 ℓa"), "H48": ("#eb6834", "s", "ovals, half-length 2.16 ℓa"),
       "H64": ("#1baf7a", "D", "ovals, half-length 2.88 ℓa"), "ell": ("#eda100", "^", "ellipses, √(A/π) = 1.27 ℓa")}
CIRC_MATCH = {"H34": 1.43, "H48": 2.15, "H64": 2.87, "ell": 1.23}


def fl(x):
    return np.nan if x in ("", "inf", "nan") else float(x)


def family(run):
    if run.startswith("W3"): return "ell"
    m = re.search(r"/(H\d+)", run)
    return m.group(1) if m and m.group(1) in FAM else None


ov = list(csv.DictReader(open(os.path.join(HERE, "wall_fourier_dyn.csv"))))
ci = list(csv.DictReader(open(os.path.join(HERE, "wall_fourier_dyn_circles_k06.csv"))))
for r in ci: r["Rla"] = round(float(r["R"]) / LA, 2)


def by_cavity(rows, key):
    G = {}
    for r in rows: G.setdefault(key(r), []).append(r)
    return G


def style(ax):
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    for s in ("left", "bottom"): ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=8); ax.grid(color=GRID, lw=0.6); ax.set_axisbelow(True)


fig, axs = plt.subplots(2, 2, figsize=(10, 7.6))
(a, b), (c, d) = axs
k6 = [r for r in ov if r["k"] == "6" and family(r["run"])]
for fam, (col, mk, lab) in FAM.items():
    rows = [r for r in k6 if family(r["run"]) == fam]
    G = by_cavity(rows, lambda r: round(float(r["aspect"]), 2))
    cr = [r for r in ci if abs(r["Rla"] - CIRC_MATCH[fam]) < 0.02]
    xs = [1.0] + sorted(G)
    for ax, col_name in ((a, "fstatic"), (b, "drift_ratio")):
        pts = [cr] + [G[x] for x in sorted(G)]
        for x, rs in zip(xs, pts):
            ax.scatter([x] * len(rs), [fl(r[col_name]) for r in rs], s=12, color=col, alpha=0.35, lw=0, zorder=2)
        ys = [np.nanmean([fl(r[col_name]) for r in rs]) for rs in pts]
        ax.plot(xs, ys, color=col, lw=2, zorder=3)
        ax.scatter(xs, ys, s=46, marker=mk, color=col, edgecolor="white", lw=1.2, zorder=4, label=lab)
a.set_ylabel("static fraction f_static (1 = pinned lattice)", color=INK, fontsize=9)
a.set_title("(a) k = 6: wall lattice pinning", loc="left", fontsize=10, color=INK)
b.axhline(0, color=INK2, lw=0.8)
b.set_ylabel("wall +½ drift / |v_w|   (< 0: against the wall flow)", color=INK, fontsize=9)
b.set_title("(b) k = 6: drift of wall +½ along the wall", loc="left", fontsize=10, color=INK)
for ax in (a, b):
    ax.set_xlabel("aspect ratio (1 = circle of matching size)", color=INK, fontsize=9); style(ax)
a.legend(frameon=False, fontsize=7.5, loc="lower left")

# (c) drift vs vacancy fraction, all runs
from scipy.stats import spearmanr
vac = lambda r: 1 - min(fl(r["Nwall"]) / (2 * int(r.get("k") or 6) - 2), 1)
allx, ally = [], []
groups = [("circ", [dict(r, k="6") for r in ci], INK2, "X", "circles k = 6, R = 0.25–2.9 ℓa")]
groups += [(fam, [r for r in ov if family(r["run"]) == fam], col, mk, lab + (" (k = 6)" if fam == "ell" else " (all k)"))
           for fam, (col, mk, lab) in FAM.items()]
for fam, rs, col, mk, lab in groups:
    x = np.array([vac(r) for r in rs]); y = np.array([fl(r["drift_ratio"]) for r in rs]); pin = np.array([fl(r["fstatic"]) >= 0.75 for r in rs])
    c.scatter(x[pin], y[pin], s=30, marker=mk, color=col, edgecolor="white", lw=0.8, zorder=4, label=lab)
    c.scatter(x[~pin], y[~pin], s=30, marker=mk, facecolor="none", edgecolor=col, lw=1.1, zorder=4)
    allx += list(x); ally += list(y)
allx, ally = np.array(allx), np.array(ally); ok = np.isfinite(ally)
rs_ = spearmanr(allx[ok], ally[ok])[0]
c.axhline(0, color=INK2, lw=0.8)
c.text(0.98, 0.97, f"Spearman ρ = {rs_:.2f} (all runs)\nfilled: pinned (f_static ≥ 0.75); open: not pinned",
       transform=c.transAxes, ha="right", va="top", fontsize=7.5, color=INK2)
c.set_xlabel("vacancy fraction of the wall lattice, 1 − N_wall/(2k − 2)", color=INK, fontsize=9)
c.set_ylabel("wall +½ drift / |v_w|", color=INK, fontsize=9)
c.set_title("(c) all runs: upstream drift follows the wall vacancies", loc="left", fontsize=10, color=INK)
c.set_ylim(-1.3, 0.35); style(c)
c.legend(frameon=False, fontsize=7, loc="lower left")

# (d) k scan
ks = [r for r in ov if r["run"].startswith("W2") or r["run"].startswith("P3/H34")]
for ba, asp, ls in (("2.0", 1.29, "-"), ("1.6", 1.51, "--")):
    rs = [r for r in ks if f"ba{ba}_" in r["run"]]
    G = by_cavity(rs, lambda r: int(r["k"]))
    kk = sorted(G)
    for col_name, col, mk, lab in (("fstatic", "#2a78d6", "o", "f_static"), ("drift_ratio", "#eb6834", "s", "drift / |v_w|")):
        y = [np.nanmean([fl(r[col_name]) for r in G[k]]) for k in kk]
        e = [np.nanstd([fl(r[col_name]) for r in G[k]]) / np.sqrt(len(G[k])) for k in kk]
        d.errorbar(kk, y, yerr=e, color=col, ls=ls, lw=1.8, marker=mk, ms=6, mec="white", capsize=2,
                   label=f"{lab}, aspect {asp}")
d.axhline(0, color=INK2, lw=0.8)
d.set_xlabel("imposed winding k", color=INK, fontsize=9); d.set_ylabel("value (seed mean ± s.e.)", color=INK, fontsize=9)
d.set_title("(d) k scan, half-length 1.53 ℓa", loc="left", fontsize=10, color=INK)
d.legend(frameon=False, fontsize=7, loc="lower right"); style(d)

fig.tight_layout()
for ext in ("png", "pdf"): fig.savefig(os.path.join(HERE, "fig_wall_dynamics." + ext), dpi=160)
print("written fig_wall_dynamics.png/.pdf")
