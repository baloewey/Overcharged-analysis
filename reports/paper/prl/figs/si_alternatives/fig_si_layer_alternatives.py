#!/usr/bin/env python3
"""Candidate SI panels: the wall layer becomes diffusive outside the small-circle limit (user TODO, 2026-10-22).
(a-d) kymographs: arclength (fraction of the perimeter) of the wall +1/2 (d < 4 cells) against time, 1000 t.u., for a small circle, a large circle,
a sparse and a crowded aspect-2 ellipse. (e) static fraction of the wall density at the lattice mode (1 = pinned lattice) against R/l_a, circles, per k
(runs/cavity/analysis_weekend/main_sweep.csv f_static, 10-seed medians). (f) the same against aspect, shape sweep, per rho (crowd_protect/stationarity).
(g) along-wall mean-squared displacement of wall +1/2 (tracked, contiguous stretches within 5 cells of the wall, lags with >= 300 samples) for the four cases. (h) layer thickening, circles: fraction of the +1/2 within
8 cells of the wall that sit at 4-8 cells, against R/l_a, per k (crowd_protect/fig1c/bands_all.npz)."""
import csv, numpy as np, matplotlib.pyplot as plt
from common import CASES, load, tracks, LA
T = "/home/baloewe/Physics/Overcharged/reports/master/theory"; CV = "/home/baloewe/Physics/Overcharged/runs/cavity/analysis_weekend"
fig = plt.figure(figsize=(19, 9)); g = fig.add_gridspec(2, 4, hspace=0.42, wspace=0.3, left=0.05, right=0.99, top=0.94, bottom=0.07)
CC = ("#08306b", "#4a90c2", "#e6550d", "#4d1a02"); msd = []
for i, ((lab, run), c) in enumerate(zip(CASES, CC)):
    D = load(run, 0.5); ax = fig.add_subplot(g[0, i]); t = D["t"]; m = t <= t[0] + 1000
    for n in np.flatnonzero(m):
        f = D["fr"][n]
        if not len(f): continue
        p = f[f[:, 2] > 0]; d, j = D["tree"].query(p[:, :2]); w = d < 4; ax.plot(np.full(w.sum(), t[n] - t[0]), D["s"][j[w]] / D["P"], ".", ms=1.2, color=c)
    ax.set_xlabel("Time (t.u.)"); ax.set_ylabel("Position along the wall, $s/P$"); ax.set_title(lab, fontsize=10.5); ax.set_ylim(0, 1); ax.set_xlim(0, 1000)
    tr = tracks(D); lags = np.arange(1, 101); acc = np.zeros(len(lags)); cnt = np.zeros(len(lags))
    for v in tr.values():
        a = np.array(v)
        for seg in np.split(a, np.flatnonzero(a[:, 3] >= 5) + 1):          # contiguous stretches within 5 cells of the wall
            seg = seg[seg[:, 3] < 5]
            if len(seg) < 20: continue
            s_ = np.unwrap(seg[:, 4] / D["P"] * 2 * np.pi) * D["P"] / (2 * np.pi)
            for q, L in enumerate(lags):
                if len(s_) > L: acc[q] += np.sum((s_[L:] - s_[:-L]) ** 2); cnt[q] += len(s_) - L
    ok = cnt >= 300; acc[~ok] = np.nan
    msd.append((lab, c, lags * (t[1] - t[0]), acc / np.maximum(cnt, 1)))
for a, s in zip(fig.axes, "abcd"): a.text(-0.12, 1.06, f"({s})", transform=a.transAxes, fontsize=14, fontweight="bold")
M = [r for r in csv.DictReader(open(CV + "/main_sweep.csv")) if r["set"] == "main_act0.07" and r["f_static"] not in ("", "nan")]
ax = fig.add_subplot(g[1, 0]); cm = plt.cm.viridis(np.linspace(0, 0.92, 9))
for k, c in zip(range(2, 11), cm):
    G = {}
    for r in M:
        if int(r["k"]) == k and float(r["rho"]) <= 0.151: G.setdefault(round(float(r["R"]) / LA, 2), []).append(float(r["f_static"]))
    xs = sorted(G); ax.plot(xs, [np.median(G[x]) for x in xs], "-o", ms=3, color=c, lw=1, label=f"{k}")
ax.axvspan(0.8, 2.2, color="0.93", zorder=0); ax.set_xlabel(r"Cavity radius $R/\ell_a$"); ax.set_ylabel("Layer static fraction"); ax.set_ylim(0, 1.05); ax.legend(title="$k$", fontsize=8, ncol=3, loc="lower left")
ax.set_title("Circles: pinned ring lost at large $R$", fontsize=10.5)
F = [r for r in csv.DictReader(open(T + "/crowd_protect/stationarity/shape_fourier.csv")) if not r.get("err")]; S = {r["run"]: r for r in csv.DictReader(open(T + "/shapes/shape_runs.csv"))}
ax = fig.add_subplot(g[1, 1]); LV = np.array([0.04, 0.06, 0.09, 0.12, 0.14]); RC = ("#fdae6b", "#f98c3a", "#e6550d", "#a63603", "#4d1a02")
for i in range(5):
    G = {}
    for r in F:
        s_ = S[r["run"]]
        if np.argmin(abs(float(s_["rho"]) - LV)) == i and s_["family"] in ("circle", "ellipse"): G.setdefault(round(float(s_["aspect"]), 2), []).append(float(r["f_lat"]))
    xs = sorted(G); ax.plot(xs, [np.median(G[x]) for x in xs], "-o", ms=3, color=RC[i], label=f"{LV[i]:.2f}")
ax.set_xlabel("Aspect ratio (ellipses)"); ax.set_ylabel("Layer static fraction"); ax.set_ylim(0, 1.05); ax.legend(title=r"$\rho$", fontsize=8, loc="lower left"); ax.set_title("Ellipses: diffusive, crowding restores", fontsize=10.5)
ax = fig.add_subplot(g[1, 2])
for lab, c, x, y in msd: ax.loglog(x, y, "-", color=c, lw=1.6, label=lab)
ax.set_xlabel("Lag (t.u.)"); ax.set_ylabel(r"Along-wall MSD (cells$^2$)"); ax.legend(fontsize=7.5, loc="upper left"); ax.set_title("Along-wall motion, wall $+\\frac{1}{2}$", fontsize=10.5)
Z = np.load(T + "/crowd_protect/fig1c/bands_all.npz"); E = Z["edges"]; dc = 0.5 * (E[:-1] + E[1:]); ax = fig.add_subplot(g[1, 3])
for k, c in zip(range(2, 11), cm):
    G = {}
    for kk, R, hp in zip(Z["k"], Z["R"], Z["Hp"]):
        if kk == k: G.setdefault(round(R, 2), []).append(hp[(dc >= 4) & (dc < 8)].sum() / max(hp[dc < 8].sum(), 1e-9))
    xs = sorted(G); ax.plot(xs, [np.median(G[x]) for x in xs], "-o", ms=3, color=c, lw=1, label=f"{k}")
ax.axvspan(0.8, 2.2, color="0.93", zorder=0); ax.set_xlabel(r"Cavity radius $R/\ell_a$"); ax.set_ylabel("Share of near-wall $+\\frac{1}{2}$ at 4–8 cells"); ax.set_title("Circles: the layer thickens at large $R$", fontsize=10.5)
for a, s in zip(fig.axes[4:], "efgh"): a.text(-0.14, 1.06, f"({s})", transform=a.transAxes, fontsize=14, fontweight="bold")
for e in ("pdf", "png"): fig.savefig(f"fig_si_layer_alternatives.{e}", dpi=170)
