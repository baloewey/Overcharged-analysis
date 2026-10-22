#!/usr/bin/env python3
"""Kinetics vs steady state with crowding (2026-10-22). (a) clean escape hazard per wall +1/2 and entry hazard per bulk +1/2 (exchange/exchange_count.csv,
tracked), elongated cavities (aspect >= 1.4) and channels (L 192, W 40). (b) bulk charge per wall length (same runs; circles added). (c) layer stationarity:
static fraction of the wall density at the lattice mode 2k-2 (stationarity/shape_fourier.csv) by aspect class. (d) Fig. 1(c) accounting, large circles:
+1/2 in the band 4-8 cells from the wall vs R/l_a (fig1c/charge_bands.txt): the wall 'vacancies' of Fig. 1(c) sit just inside the 4-cell band."""
import csv, numpy as np, matplotlib.pyplot as plt
LV = np.array([0.04, 0.06, 0.09, 0.12, 0.14]); R = [r for r in csv.DictReader(open("exchange/exchange_count.csv")) if not r.get("err")]
S = {r["run"]: r for r in csv.DictReader(open("../shapes/shape_runs.csv"))}; f = lambda r, k: float(r[k])
fig, ax = plt.subplots(1, 4, figsize=(19, 4.6)); fig.subplots_adjust(wspace=0.33, left=0.05, right=0.99, bottom=0.15, top=0.9)
def cav(lo, hi):
    out = []
    for i in range(5):
        rs = [r for r in R if r["kind"] == "cavity" and lo <= float(S[r["run"].split("shapes_server/")[1]]["aspect"]) <= hi and np.argmin(abs(f(r, "rho") - LV)) == i]
        out.append(rs)
    return out
el = cav(1.39, 9); ci = cav(0.99, 1.01); ch = {}
for r in R:
    if r["kind"] == "channel": ch.setdefault(int(r["m"]), []).append(r)
ms = sorted(ch); rch = np.array([2 * m / 192 for m in ms])
he = lambda rs: np.median([1000 * f(r, "esc") / max(f(r, "expW"), 1) for r in rs]); hi_ = lambda rs: np.median([1000 * f(r, "ent") / max(f(r, "expB"), 1) for r in rs])
qb = lambda rs: np.median([f(r, "qb") / f(r, "Lw") for r in rs])
ax[0].semilogy(LV, [he(x) for x in el], "-o", color="#a63603", label="Escape, cavities (aspect ≥ 1.4)"); mz = [m for m in ms if he(ch[m]) > 0]; ax[0].semilogy([2 * m / 192 for m in mz], [he(ch[m]) for m in mz], "-s", color="#a63603", mfc="w", label="Escape, channels")
ax[0].semilogy(LV, [hi_(x) for x in el], "-o", color="#3b6fb6", label="Entry, cavities"); ax[0].semilogy(rch, [hi_(ch[m]) for m in ms], "-s", color="#3b6fb6", mfc="w", label="Entry, channels")
ax[0].set_xlabel(r"Wall density $\rho=2k/P$"); ax[0].set_ylabel("Hazard per defect per 1000 t.u."); ax[0].legend(loc="lower left", fontsize=9); ax[0].set_ylim(2e-3, 60)
ax[1].plot(LV, [qb(x) for x in el], "-o", color="#a63603", label="Cavities, aspect ≥ 1.4"); ax[1].plot(LV, [qb(x) for x in ci], "-o", color="#3b6fb6", label="Circles")
ax[1].plot(rch, [qb(ch[m]) for m in ms], "-s", color="0.3", mfc="w", label="Channels"); ax[1].set_ylim(0, 0.011)
ax[1].set_xlabel(r"Wall density $\rho$"); ax[1].set_ylabel("Bulk charge per wall length"); ax[1].legend(loc="lower left", fontsize=9.5)
F = [r for r in csv.DictReader(open("stationarity/shape_fourier.csv")) if not r.get("err")]
for (nm, lo, hi), c in zip((("Circles", 1.0, 1.0), ("Aspect 1.1–1.3", 1.05, 1.35), ("Aspect 1.35–1.7", 1.35, 1.75), ("Aspect 1.9–2.5", 1.85, 3)), ("#08306b", "#4a90c2", "#e6550d", "#a63603")):
    y = [np.median([float(r["f_lat"]) for r in F if lo <= float(S[r["run"]]["aspect"]) <= hi and np.argmin(abs(float(S[r["run"]]["rho"]) - LV)) == i]) for i in range(5)]
    ax[2].plot(LV, y, "-o", color=c, label=nm)
ax[2].set_xlabel(r"Wall density $\rho$"); ax[2].set_ylabel("Layer stationarity (static fraction)"); ax[2].set_ylim(0, 1.05); ax[2].legend(loc="lower right", fontsize=9.5)
rows = [l.split("|") for l in open("fig1c/charge_bands.txt").read().splitlines()[1:]]
for k, c in ((6, "#4a90c2"), (8, "#e6550d"), (10, "#4d1a02")):
    rr = sorted((float(x[1]), float(x[3]), float(x[4].split(",")[0])) for x in rows if int(x[0]) == k)
    ax[3].plot([a for a, b, d in rr], [d for a, b, d in rr], "-o", color=c, label=f"$k={k}$, +1/2 at 4–8 cells"); ax[3].plot([a for a, b, d in rr], [-b for a, b, d in rr], "--", color=c, lw=1, label=f"$k={k}$, wall deficit")
ax[3].axvspan(0.8, 2.2, color="0.93", zorder=0); ax[3].set_xlabel(r"Cavity radius $R/\ell_a$"); ax[3].set_ylabel("Defects (time average)"); ax[3].legend(loc="upper left", fontsize=8.5, ncol=1)
tit = ["Exits blocked, entries not", "Steady bulk charge", "Layer stationarity", "Fig. 1(c) wall deficit"]
for a, s, t in zip(ax, "abcd", tit): a.text(-0.13, 1.04, f"({s})", transform=a.transAxes, fontsize=14, fontweight="bold"); a.set_title(t, fontsize=12)
for e in ("pdf", "png"): fig.savefig(f"fig_exchange.{e}", dpi=200)
