#!/usr/bin/env python3
"""Alternatives for PRL Fig. 1(c) (user, 2026-10-22: 'too many caveats; large k has more defects to leak, maybe a normalization'). 10-seed circle sweep
(main_act0.07, +k, k 2-10, rho <= 0.15), per-run wall-distance profiles (reports/master/theory/crowd_protect/fig1c/bands_all.npz); all counts use the same
bands, so every panel adds up: wall = d < 4 cells. Medians over seeds per (k, R).
(A) vacant fraction f_vac = 1 - N_wall/(2k-2) vs R/l_a, per k (normalized).  (B) the current two quantities on one partition: charge outside the wall band
Q_b = k - N_wall/2 (top) and f_vac (bottom); no knee.  (C) occupancy of a thick layer, N(d < 8)/(2k-2): the deficit is absorbed by thickening.
(D) where the wall charge goes, k = 6 and 10: fraction of the 2k-2 wall sites' charge found at d < 4, 4-8 and > 8 cells.  (E) charge farther than d from the
wall, k = 6, five radii (the knee picture made explicit).  (F) compact: f_vac vs R/l_a, median over k with the k range shaded."""
import numpy as np, matplotlib.pyplot as plt
from prlstyle import ROOT
Z = np.load(ROOT + "/reports/master/theory/crowd_protect/fig1c/bands_all.npz"); E = Z["edges"]; dc = 0.5 * (E[:-1] + E[1:]); K, Rr, Hp, Hm = Z["k"], Z["R"], Z["Hp"], Z["Hm"]
q = 0.5 * (Hp - Hm); keep = Z["rho"] <= 0.151
nw4 = Hp[:, dc < 4].sum(1); n8 = Hp[:, dc < 8].sum(1); fv = 1 - nw4 / (2 * K - 2); occ8 = n8 / (2 * K - 2); qb = q[:, dc >= 4].sum(1)
cm = plt.cm.viridis(np.linspace(0, 0.92, 9)); KS = range(2, 11)
def med(y, k):
    m = keep & (K == k); G = {}
    for R, v in zip(Rr[m], y[m]): G.setdefault(round(R, 2), []).append(v)
    xs = np.array(sorted(G)); return xs, np.array([np.median(G[x]) for x in xs])
fig = plt.figure(figsize=(16, 9.4)); g = fig.add_gridspec(2, 3, hspace=0.36, wspace=0.28, left=0.06, right=0.99, top=0.94, bottom=0.07)
def win(a): a.axvspan(0.8, 2.2, color="0.93", lw=0, zorder=0)
ax = fig.add_subplot(g[0, 0])
for k, c in zip(KS, cm): x, y = med(fv, k); ax.plot(x, y, "-o", ms=2.5, lw=0.9, color=c, label=f"{k}")
win(ax); ax.axhline(0, color="0.5", ls="--", lw=0.7); ax.set_ylim(-0.4, 0.35); ax.set_xlabel(r"Cavity radius $R/\ell_a$"); ax.set_ylabel(r"Vacant fraction $1-N_\mathrm{wall}/(2k-2)$")
ax.legend(title="$k$", ncol=3, fontsize=8, loc="lower right"); ax.set_title("A: normalised wall count", fontsize=11, loc="left")
sub = g[0, 1].subgridspec(2, 1, hspace=0.08); a1 = fig.add_subplot(sub[0]); a2 = fig.add_subplot(sub[1], sharex=a1)
for k, c in zip(KS, cm):
    x, y = med(qb, k); a1.plot(x, y, "-o", ms=2.5, lw=0.9, color=c); x, y = med(fv, k); a2.plot(x, y, "-o", ms=2.5, lw=0.9, color=c)
for a in (a1, a2): win(a)
a1.axhline(1, color="0.5", ls="--", lw=0.7); a1.set_ylim(0.5, 3.5); a1.set_ylabel(r"$Q_b=k-N_\mathrm{wall}/2$"); plt.setp(a1.get_xticklabels(), visible=False)
a2.axhline(0, color="0.5", ls="--", lw=0.7); a2.set_ylim(-0.4, 0.35); a2.set_ylabel(r"$f_\mathrm{vac}$"); a2.set_xlabel(r"Cavity radius $R/\ell_a$")
a1.set_title("B: current pair, one partition (no knee)", fontsize=11, loc="left")
ax = fig.add_subplot(g[0, 2])
for k, c in zip(KS, cm): x, y = med(occ8, k); ax.plot(x, y, "-o", ms=2.5, lw=0.9, color=c)
win(ax); ax.axhline(1, color="0.5", ls="--", lw=0.7); ax.set_xlabel(r"Cavity radius $R/\ell_a$"); ax.set_ylabel(r"$N(d<8)/(2k-2)$"); ax.set_title("C: thick-layer occupancy (8 cells)", fontsize=11, loc="left")
ax = fig.add_subplot(g[1, 0])
for k, ls in ((6, "-"), (10, "--")):
    tot = (2 * k - 2) / 2
    for band, sel, c in (("$d<4$", dc < 4, "#08306b"), ("$4\\leq d<8$", (dc >= 4) & (dc < 8), "#e6550d"), ("$d\\geq8$", dc >= 8, "#6a9f58")):
        y = q[:, sel].sum(1) / tot; x, yy = med(y, k); ax.plot(x, yy, ls, color=c, lw=1.5, marker="o", ms=2.5, label=f"$k={k}$, {band}")
win(ax); ax.set_xlabel(r"Cavity radius $R/\ell_a$"); ax.set_ylabel("Charge in band / $(k-1)$"); ax.legend(fontsize=7.5, ncol=2, loc="center left"); ax.set_title("D: where the wall charge sits", fontsize=11, loc="left")
ax = fig.add_subplot(g[1, 1]); RS = (0.95, 1.43, 1.91, 2.39, 2.87); bc = plt.cm.Blues(np.linspace(0.95, 0.35, 5))
for R0, c in zip(RS, bc):
    m = keep & (K == 6) & (np.abs(Rr - R0) < 0.02)
    if not m.any(): continue
    enc = 6 - np.cumsum(q[m].mean(0)); ax.plot(np.r_[0, E[1:-1]], np.r_[6, enc[:-1]], "-", color=c, lw=1.6, label=f"$R={R0}\\,\\ell_a$")
ax.axhline(1, color="0.5", ls=":", lw=0.8); ax.axvline(4, color="0.4", ls="--", lw=0.8); ax.set_xlim(0, 20); ax.set_ylim(-1, 6.3)
ax.set_xlabel("Distance from the wall $d$ (cells)"); ax.set_ylabel("Charge farther than $d$"); ax.legend(fontsize=8.5, title="$k=6$", title_fontsize=8.5); ax.set_title("E: enclosed-charge profiles", fontsize=11, loc="left")
ax = fig.add_subplot(g[1, 2]); allx = sorted({round(R, 2) for R in Rr[keep]}); lo, md, hi = [], [], []
for x in allx:
    v = [np.median(fv[keep & (K == k) & (np.round(Rr, 2) == x)]) for k in KS if (keep & (K == k) & (np.round(Rr, 2) == x)).any()]
    lo.append(min(v)); md.append(np.median(v)); hi.append(max(v))
ax.fill_between(allx, lo, hi, color="#9ecae1", alpha=0.6, lw=0, label="Range over $k$"); ax.plot(allx, md, "-", color="#08306b", lw=1.8, label="Median over $k$")
win(ax); ax.axhline(0, color="0.5", ls="--", lw=0.7); ax.set_ylim(-0.4, 0.35); ax.set_xlabel(r"Cavity radius $R/\ell_a$"); ax.set_ylabel(r"$f_\mathrm{vac}$"); ax.legend(fontsize=9, loc="lower right")
ax.set_title("F: compact single curve", fontsize=11, loc="left")
for e in ("pdf", "png"): fig.savefig(f"fig1c_alternatives.{e}", dpi=170)
