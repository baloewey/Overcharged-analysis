#!/usr/bin/env python3
"""PRL Fig. 4 (beyond circles), full page width, redrawn from cached data (2026-10-15). Run: figs/run.sh fig4_shapes.py
(a) tilt coherence vs aspect, (b) vacancies vs aspect: seed and k medians per (family, aspect, R_eq) from shapes/shape_runs.csv (as phase_diagram.py).
(c) wall law: bulk charge per unit wall length vs winding density q (energetics/perlength.npy, as fig_perlength c), channel point (runs/channel, q 0.031).
(d) wall blocks the inward jet: transmitted fraction g(d) of the inward component of a +1/2's own jet vs depth (energetics/jetmap_ff12.txt circle,
    jetmap_oval.txt ovals; parsed from the text output).
(e) untilting run-up: cos tau of escapers vs resident controls against time relative to departure (energetics/etsg_shapes.npy), elongated cavities (aspect >= 1.4) pooled.
(f) shielding: potential felt by one bulk +1/2 from the measured wall layer (one-constant, K = 32.5 kBT; bulkdyn/shield_wl/maps_*_meas.csv, layer with the
    median leftover range), circle vs aspect-2 ellipse, shared square-root colour scale.
(g) mirror control placeholder (server mirror shape sweep pending)."""
import os, re, json, csv, numpy as np
from collections import defaultdict
from prlstyle import plt, PAGE, RED, BLUE, GREY, DARK, label, save, ROOT
from matplotlib.colors import PowerNorm
from matplotlib.lines import Line2D
import matplotlib

T = ROOT + "/reports/master/theory"

# ---------------- data
rows = list(csv.DictReader(open(T + "/shapes/shape_runs.csv")))
G = defaultdict(lambda: defaultdict(list))
for r in rows:
    G[(r["family"], round(float(r["aspect"]), 2), r["Req_la"][:3])]["coh"].append(float(r["coh"]))
    G[(r["family"], round(float(r["aspect"]), 2), r["Req_la"][:3])]["vac"].append(float(r["vac"]))

A = [tuple(a) for a in np.load(T + "/energetics/perlength.npy", allow_pickle=True)]   # family, aspect, Req, P, k, nwall, Qb
CH_Q, CH_QB, CH_QBE = 0.031, 0.0132, 0.0004


def gline(tag):
    s = open(T + f"/energetics/jetmap_{tag}.txt").read().split("g(d) = ")[1].splitlines()[0]
    s = s.split(":", 1)[1] if s.startswith("<") else s
    pr = re.findall(r"([\d.]+):([\d.]+)", s.split("^2>:")[-1]); return np.array([[float(a), float(b)] for a, b in pr])


E = np.load(T + "/energetics/etsg_shapes.npy", allow_pickle=True)
kind = E[:, 1].astype(str); lag = E[:, 2].astype(int); ct = E[:, 6].astype(float); asp_e = np.round(E[:, 7].astype(float), 2)

J = json.load(open(T + "/bulkdyn/shield_wl/layers.json")); A2, B2, R1 = J["A2"], J["B2"], J["R1"]; PIK = np.pi * 32.5


def maps(name):
    a = np.loadtxt(T + f"/bulkdyn/shield_wl/maps_{name}.csv", delimiter=","); out = []
    for l in np.unique(a[:, 0]):
        b = a[a[:, 0] == l]; out.append((b[:, 1], b[:, 2], b[:, 3] * PIK, (b[:, 3] - b[:, 4]) * PIK))
    return out


def median_layer(M, lev):
    r = [np.ptp(vl[lev(x, y) < 0.6]) for x, y, v, vl in M]; return int(np.argsort(r)[len(r) // 2]), np.median(r)


def grid(x, y, v):
    xs, ys = np.unique(x), np.unique(y); Gr = np.full((len(xs), len(ys)), np.nan)
    Gr[np.searchsorted(xs, x), np.searchsorted(ys, y)] = v; return xs, ys, Gr


# ---------------- figure
fig = plt.figure(figsize=(PAGE, 5.0))
gs0 = fig.add_gridspec(2, 1, height_ratios=[1, 1.05], hspace=0.5)
gt = gs0[0].subgridspec(1, 3, wspace=0.45); gb = gs0[1].subgridspec(1, 4, width_ratios=[1, 1, 1.55, 0.75], wspace=0.5)
cols = {"1.0": "#9ecae1", "1.5": "#4a90c2", "2.0": "#08306b"}
for j, (c, yl) in enumerate((("coh", "Tilt coherence $C$"), ("vac", "Vacancies, $2k-2-N_\\mathrm{wall}$"))):
    ax = fig.add_subplot(gt[0, j])
    for fam, ls, mk in (("ellipse", "-", "o"), ("oval", "--", "s")):
        for req in ("1.0", "1.5", "2.0"):
            ks = sorted(k for k in G if k[0] == fam and k[2] == req)
            xs = [1.0] + [k[1] for k in ks]; ys = [np.nanmedian(G[("circle", 1.0, req)][c])] + [np.nanmedian(G[k][c]) for k in ks]
            ax.plot(xs, ys, ls=ls, marker=mk, ms=2.8, color=cols[req], lw=1.0)
    ax.set_xlabel("Aspect ratio"); ax.set_ylabel(yl); ax.set_xticks([1.0, 1.5, 2.0, 2.5])
    if c == "coh":
        ax.axhspan(0, 0, color="none"); ax.set_ylim(0.3, 1.0)
        h = [Line2D([], [], color=cols[r], lw=1.2, label=fr"$R_\mathrm{{eq}}={r}\,\ell_a$") for r in ("1.0", "1.5", "2.0")]
        h += [Line2D([], [], color=DARK, ls="-", marker="o", ms=2.8, label="Ellipse"), Line2D([], [], color=DARK, ls="--", marker="s", ms=2.8, label="Cassini oval")]
        ax.legend(handles=h, loc="upper right", fontsize=6, ncol=1, borderaxespad=0.1)
    label(ax, "ab"[j])

# (c) wall law
ax = fig.add_subplot(gt[0, 2])
q = np.array([a[4] / a[3] for a in A]); qb = np.array([a[6] / a[3] for a in A]); asp = np.array([a[1] for a in A]); fam = np.array([a[0] for a in A])
el = asp >= 1.5
ax.scatter(q[~el], qb[~el], s=5, color="0.75", lw=0, label="Aspect $<1.5$", rasterized=True)
ax.scatter(q[el], qb[el], s=5, color=BLUE, lw=0, alpha=0.8, label="Aspect $\\geq1.5$", rasterized=True)
ax.errorbar([CH_Q], [CH_QB], yerr=[CH_QBE], fmt="D", ms=4.5, color="w", mec="k", mew=0.8, ecolor="k", capsize=2, zorder=5, label="Channel (flat limit)")
ax.set_xlabel("Winding density $q=k/P$"); ax.set_ylabel("Bulk charge per length $Q_b/P$")
ax.set_yticks([0.005, 0.01, 0.015])
ax.legend(loc="upper right", fontsize=6, handletextpad=0.2, borderaxespad=0.1); ax.set_ylim(0.0025, 0.0195); ax.set_xlim(0.012, 0.078)
label(ax, "c")

# (d) jet transmission
ax = fig.add_subplot(gb[0, 0])
for tag, col, lab in (("ff12", RED, "Circle"), ("oval", BLUE, "Ovals")):
    g = gline(tag); ax.plot(g[:, 0], g[:, 1], "-o", ms=2.8, color=col, label=lab)
ax.axhline(1, color="0.8", lw=0.6); ax.axvspan(0, 2, color="0.92", lw=0, zorder=0)
ax.set_xlabel("Depth $d$ (cells)"); ax.set_ylabel("Inward jet transmitted $g(d)$"); ax.set_ylim(0, 1.15); ax.set_xlim(0, 12.5)
ax.legend(loc="lower right", fontsize=6); label(ax, "d")

# (e) run-up
ax = fig.add_subplot(gb[0, 1])
for kd, col, mk, lab in (("esc", RED, "o", "Escaper"), ("ctl", GREY, "s", "Control")):
    sel = (kind == kd) & (asp_e >= 1.4); ls_ = sorted(set(lag[sel]))
    m = [np.nanmean(ct[sel & (lag == l)]) for l in ls_]; e = [np.nanstd(ct[sel & (lag == l)]) / np.sqrt(np.isfinite(ct[sel & (lag == l)]).sum()) for l in ls_]
    ax.errorbar(2 * np.array(ls_), m, e, color=col, marker=mk, ms=2.8, capsize=1.5, lw=1.0, ls="-" if kd == "esc" else "--", label=lab)
ax.axvline(0, color="0.7", lw=0.6, ls=":")
ax.set_xlabel("Time to departure (t.u.)"); ax.set_ylabel(r"Tilt projection $\cos\tau$"); ax.legend(loc="upper left", fontsize=6)
label(ax, "e")

# (f) shielding maps
Mc, Me = maps("circle_meas"), maps("ellipse_meas")
ic, rc = median_layer(Mc, lambda x, y: np.hypot(x, y) / R1); ie, re_ = median_layer(Me, lambda x, y: np.hypot(x / A2, y / B2))
th = np.linspace(0, 2 * np.pi, 400); VM = 300; norm = PowerNorm(0.5, vmin=0, vmax=VM)
sub = gb[0, 2].subgridspec(1, 2, width_ratios=[1.0, 1.42], wspace=0.05)
for i, (Mm, k, lay, A_, B_, ttl) in enumerate(((Mc, ic, J["circle_meas"], R1, R1, "Circle"), (Me, ie, J["ellipse_meas"], A2, B2, "Ellipse 2"))):
    ax = fig.add_subplot(sub[0, i]); x, y, v, _ = Mm[k]; xs, ys, Gr = grid(x, y, v - np.nanmin(v))
    im = ax.pcolormesh(xs, ys, np.minimum(Gr, VM).T, cmap="magma_r", shading="nearest", norm=norm, rasterized=True)
    ax.contour(xs, ys, Gr.T, levels=np.arange(20, VM, 20), colors="w", linewidths=0.3, alpha=0.7)
    ax.plot(A_ * np.cos(th), B_ * np.sin(th), "k-", lw=0.8); z = np.array(lay[k]); ax.plot(z[:, 0], z[:, 1], "o", ms=2.6, color=BLUE, mec="w", mew=0.4)
    ax.set_aspect("equal"); ax.axis("off"); ax.set_xlim(-A2 * 1.05, A2 * 1.05) if i else ax.set_xlim(-R1 * 1.08, R1 * 1.08)
    ax.set_title(f"{ttl}\n{len(lay[k])} wall $+1/2$", fontsize=6, pad=1, linespacing=1.1)
    if i == 0: label(ax, "f", x=-0.02, y=1.12)
cax = fig.add_axes([0, 0, 1, 1]); fig.canvas.draw()
p0 = fig.axes[-3].get_position(); p1 = fig.axes[-2].get_position()
cax.set_position([p0.x0 + 0.01, p0.y0 - 0.03, p1.x1 - p0.x0 - 0.02, 0.014])
cb = fig.colorbar(im, cax=cax, orientation="horizontal", ticks=[0, 20, 80, 180, 300]); cb.ax.tick_params(labelsize=6, length=2)
cb.set_label("Layer potential on a bulk $+1/2$ ($k_BT$)", fontsize=6.5, labelpad=1)

# (g) mirror placeholder
ax = fig.add_subplot(gb[0, 3]); ax.set_xlim(1, 2.5); ax.set_ylim(0, 5); ax.set_xticks([1.0, 2.0]); ax.set_yticks([0, 2, 4])
ax.set_xlabel("Aspect ratio"); ax.set_ylabel("Vacancies"); ax.text(1.75, 2.5, "Mirror control\n$+k$ vs $-k$\n(pending)", ha="center", va="center", fontsize=6.5, color="0.5")
label(ax, "g")
save(fig, "fig4_shapes")
print(f"shielding: circle layer {ic} range0.6 {rc:.1f} kBT (n wall {len(J['circle_meas'][ic])}); ellipse layer {ie} range0.6 {re_:.1f} kBT (n wall {len(J['ellipse_meas'][ie])})")
