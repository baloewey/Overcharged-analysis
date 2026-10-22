#!/usr/bin/env python3
"""What lowers tilt coherence in crowded layers (2026-10-22). Shape sweep (675 runs; circles vs aspect >= 1.4), frames every 2 t.u., wall +1/2 tracked.
(a) coherence C = |<exp(i tilt)>| per frame (median over runs) vs rho, and C predicted from the wrong-sign fraction m and the majority tilt tau alone,
    sqrt(cos^2 tau + (1-2m)^2 sin^2 tau). (b) tilt histograms, circles, signed relative to the frame's majority sign. (c) kinks start more often and
(d) last longer when crowded (episodes of a tracked defect with the minority sign). (e) the near-wall flow does NOT get less coherent: |<v_t>|/<|v_t|> at
1-3 cells (wallflow.npz). (f) elastic cost of a one-defect kink whose defect slides to its own lock site (kink_energy.py; one constant K 44), vs N for three
tilts; markers: the measured (rho, tilt) of the circles."""
import csv, numpy as np, matplotlib.pyplot as plt
Z = np.load("kink_life.npz", allow_pickle=True); SR = {r["run"]: r for r in csv.DictReader(open("../shapes/shape_runs.csv"))}
LV = np.array([0.04, 0.06, 0.09, 0.12, 0.14]); RCOL = ("#fdae6b", "#f98c3a", "#e6550d", "#a63603", "#4d1a02"); runs = Z["runs"]; HB = Z["HB"]; hc = 0.5 * (HB[:-1] + HB[1:])
rho = np.array([float(SR[r]["rho"]) for r in runs]); asp = np.array([float(SR[r]["aspect"]) for r in runs]); lev = np.argmin(np.abs(rho[:, None] - LV), 1)
CL = (("Circles", asp < 1.05, "#3b6fb6", "o"), ("Aspect ≥ 1.4", asp >= 1.4, "#a63603", "s"))
fig, ax = plt.subplots(2, 3, figsize=(15, 9)); fig.subplots_adjust(wspace=0.32, hspace=0.42, left=0.07, right=0.98, top=0.94, bottom=0.08); ax = ax.ravel()
peak = {}
for nm, sel, c, mk in CL:
    Cm, Cp, rate, dur = [], [], [], []
    for i in range(5):
        m = np.flatnonzero((lev == i) & sel); F = [Z["fstat"][j] for j in m]; wf = np.median([1 - f[:, 1].mean() for f in F])
        h = Z["hist"][m].sum(0); tau = np.radians(hc[hc > 0][np.argmax(h[hc > 0])]); peak[(nm, i)] = np.degrees(tau)
        Cm.append(np.median([f[:, 0].mean() for f in F])); Cp.append(np.sqrt(np.cos(tau) ** 2 + (1 - 2 * wf) ** 2 * np.sin(tau) ** 2))
        rate.append(np.median(Z["rate"][m])); E = np.concatenate([Z["eps"][j] for j in m]); dur.append(E.mean())
    ax[0].plot(LV, Cm, "-" + mk, color=c, ms=6, label=f"{nm}, measured"); ax[0].plot(LV, Cp, ":" + mk, color=c, ms=5, mfc="w", label=f"{nm}, from wrong-sign fraction")
    ax[2].plot(LV, rate, "-" + mk, color=c, ms=6, label=nm); ax[3].plot(LV, dur, "-" + mk, color=c, ms=6, label=nm)
ax[0].set_ylabel("Tilt coherence $C$"); ax[0].set_ylim(0.3, 1.05); ax[0].legend(loc="center left", bbox_to_anchor=(0.0, 0.64), fontsize=8.5)
for i in range(5):
    m = np.flatnonzero((lev == i) & (asp < 1.05)); h = Z["hist"][m].sum(0); h = h / h.sum() / 5
    ax[1].semilogy(hc, np.maximum(h, 1e-6), "-", color=RCOL[i], lw=1.6, label=fr"$\rho={LV[i]:.2f}$")
ax[1].axvline(0, color="0.8", lw=0.8, zorder=0); ax[1].set_xlabel("Tilt relative to the majority sign (deg)"); ax[1].set_ylabel("Probability density (1/deg)")
ax[1].set_xlim(-180, 180); ax[1].set_xticks([-180, -90, 0, 90, 180]); ax[1].set_ylim(1e-6, 0.2); ax[1].legend(loc="upper left", fontsize=9.5, ncol=1)
ax[1].text(-95, 2e-6, "Wrong sign", ha="center", fontsize=10, color="0.35"); ax[1].text(120, 2e-6, "Majority", ha="center", fontsize=10, color="0.35")
ax[2].set_ylabel("Kink starts per defect per 1000 t.u."); ax[2].set_yscale("log"); ax[2].legend(loc="center right", fontsize=9.5)
ax[3].set_ylabel("Mean kink lifetime (t.u.)"); ax[3].legend(loc="lower right", fontsize=9.5)
W = np.load("wallflow.npz"); wr = {r: a for r, a in zip(W["runs"], W["prof"])}; dc = 0.5 * (W["DB"][:-1] + W["DB"][1:]); band = (dc > 1) & (dc < 3)
for nm, sel, c, mk in CL[:1]:
    y = [np.median([abs(wr[r][4][band].mean()) / wr[r][2][band].mean() for r in runs[(lev == i) & sel]]) for i in range(5)]
    ax[4].plot(LV, y, "-" + mk, color=c, ms=6, label=nm)
ax[4].set_ylabel(r"Near-wall flow coherence $|\langle v_t\rangle|/\langle|v_t|\rangle$"); ax[4].set_ylim(0, 1.05); ax[4].text(0.97, 0.05, "Circles, 1–3 cells from the wall", transform=ax[4].transAxes, ha="right", fontsize=10)
K = np.load("kink_energy.npz"); NS = K["NS"]; PP = 2 * np.pi * 31.830989
for tau, c in zip(K["TAUS"], ("#9ecae1", "#4a90c2", "#08306b")): ax[5].plot(NS / PP, [K[f"N{N}_t{tau}"][1] for N in NS], "-", color=c, lw=1.6, label=f"Tilt {tau}°")
for i in range(5):
    t = peak[("Circles", i)]; N = int(round(LV[i] * PP / 2) * 2); N = min(max(N, 6), 28)
    ax[5].plot(LV[i], np.interp(t, K["TAUS"], [K[f"N{N}_t{tt}"][1] for tt in K["TAUS"]]), "D", color=RCOL[i], ms=8, mec="k")
ax[5].set_xlabel(r"Wall density $\rho=N/P$"); ax[5].set_ylabel(r"Kink cost ($k_BT$)"); ax[5].legend(loc="upper right", fontsize=9.5); ax[5].set_ylim(0, 40)
ax[5].text(0.03, 0.04, "Diamonds: measured tilt of the circles", transform=ax[5].transAxes, fontsize=9.5)
for a in (ax[0], ax[2], ax[3], ax[4]): a.set_xlabel(r"Wall density $\rho=2k/P$"); a.set_xticks(LV)
tit = ["Coherence and its cause", "Tilt distribution, circles", "Kink frequency", "Kink lifetime", "Flow at the wall stays coherent", "Elastic cost of a kink"]
for a, s, t in zip(ax, "abcdef", tit): a.text(-0.14, 1.04, f"({s})", transform=a.transAxes, fontsize=14, fontweight="bold"); a.set_title(t, fontsize=12)
for e in ("pdf", "png"): fig.savefig(f"fig_coherence.{e}", dpi=200)
