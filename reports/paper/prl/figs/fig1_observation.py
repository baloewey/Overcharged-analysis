#!/usr/bin/env python3
"""PRL Fig. 1 (observation). Full page width. Run: figs/run.sh fig1_observation.py
(a) Time-averaged defect density, k = 6: passive (R = 1.43 l_a) and active at R = 0.86, 1.43, 2.15 l_a (reports/master/data/dist_*.npz, written by
    reports/master/scripts/fig_distributions.py from runs/cavity/server_main; +1/2 blue, -1/2 orange, shared log scale, defects per cell^2 per frame).
(b) Anchoring and the lock, k = 4: wall director theta = k phi, the 2k - 2 = 6 homeotropic points, and a ring of 6 wall +1/2 sharing a common tilt tau
    by a registry shift tau/(2k - 2) (lock tau_j = (2k - 2) phi_j + const).
(c) 10-seed circle sweep (act 0.07, lambda 2, +k, rho <= 0.15, k = 2-10; l_a = 22.21): bulk charge (knee of Q(r), primary criterion;
    runs/cavity/analysis_weekend/criteria_all.csv q_knee, mean +- s.e. over the seeds that have a knee; points with a knee in < 5 of 10 seeds left out)
    and wall +1/2 count minus (2k - 2) (main_sweep.csv N_wall), against R/l_a; window 0.8-2.2 shaded.
(d) Snapshot of an active circle, k = 6, R = 1.43 l_a (main sweep rho 0.06, s05, t = 9000): flow vorticity (flow smoothed over 2 cells), flow arrows,
    +1/2 comets (red, dot = head), -1/2 tripods (blue)."""
import os, sys, csv, json, numpy as np
from matplotlib.colors import LogNorm
from matplotlib.path import Path
from scipy.ndimage import gaussian_filter
from prlstyle import plt, PAGE, RED, BLUE, GREY, DARK, ROOT, label, save
sys.path.insert(0, ROOT + "/runs/cavity"); sys.path.insert(0, ROOT + "/reports/master/theory/energetics")
from edge_criteria import read_with_angle  # noqa
from selfflow import load_flow  # noqa
LA = 22.21; DATA = ROOT + "/reports/master/data"; CAV = ROOT + "/runs/cavity"


def sweep():
    f = lambda r, c: float(r[c]) if r[c] not in ("", "nan") else np.nan
    C = [r for r in csv.DictReader(open(CAV + "/analysis_weekend/criteria_all.csv")) if r["set"] == "pos"]
    M = [r for r in csv.DictReader(open(CAV + "/analysis_weekend/main_sweep.csv")) if r["set"] == "main_act0.07"]
    nw = {r["run"]: f(r, "N_wall") for r in M}
    G = {}
    for r in C:
        k = int(r["k"])
        if not 2 <= k <= 10 or f(r, "rho") > 0.151: continue
        G.setdefault((k, round(f(r, "R"), 2)), []).append((f(r, "q_knee") if f(r, "has_knee") else np.nan, nw.get(r["run"], np.nan) - (2 * k - 2)))
    out = []
    for (k, R), v in sorted(G.items()):
        v = np.array(v); q = v[:, 0][np.isfinite(v[:, 0])]; w = v[:, 1][np.isfinite(v[:, 1])]
        se = lambda a: a.std(ddof=1) / np.sqrt(len(a)) if len(a) > 1 else 0.0
        out.append((k, R / LA, q.mean() if len(q) >= 0.5 * len(v) else np.nan, se(q) if len(q) >= 0.5 * len(v) else np.nan, w.mean() if len(w) else np.nan, se(w) if len(w) else np.nan, len(v)))
    return np.array(out)


def maps(axs):
    items = [("Passive\n$R=1.43\\,\\ell_a$", "main_passive__rho0.06__k06__s01"), ("Active\n$R=0.86\\,\\ell_a$", "main_act0.07__rho0.10__k06__s01"),
             ("Active\n$R=1.43\\,\\ell_a$", "main_act0.07__rho0.06__k06__s01"), ("Active\n$R=2.15\\,\\ell_a$", "main_act0.07__rho0.04__k06__s01")]
    D = [np.load(f"{DATA}/dist_cavity__server_main__{p}.npz") for _, p in items]
    vmax = max(max(d["Hp"].max(), d["Hm"].max()) for d in D); hx = max(np.ptp(d["W"][:, 0]) for d in D) / 2 + 1.5; ims = []
    for ax, (lab, _), d in zip(axs, items, D):
        ext = [d["xe"][0], d["xe"][-1], d["ye"][0], d["ye"][-1]]
        for H, cm in ((d["Hm"], "Oranges"), (d["Hp"], "Blues")):
            ims.append(ax.imshow(np.where(H > 0, H, np.nan).T, origin="lower", extent=ext, cmap=cm, interpolation="nearest", norm=LogNorm(vmin=vmax * 1e-3, vmax=vmax), rasterized=True))
        ax.plot(np.r_[d["W"][:, 0], d["W"][0, 0]], np.r_[d["W"][:, 1], d["W"][0, 1]], color="k", lw=0.6)
        ax.set_xlim(-hx, hx); ax.set_ylim(-hx, hx); ax.set_aspect("equal"); ax.axis("off")
        ax.text(0, hx * 1.02, lab, ha="center", va="bottom", fontsize=6.5)
        ax.text(0, -hx * 1.02, f"{float(d['Np']):.1f} $+\\frac{{1}}{{2}}$, {float(d['Nm']):.1f} $-\\frac{{1}}{{2}}$", ha="center", va="top", fontsize=6)
    return ims, vmax


def lock(ax, k=4, tau=np.radians(60)):
    ph = np.linspace(0, 2 * np.pi, 300); ax.plot(np.cos(ph), np.sin(ph), color=DARK, lw=0.9)
    for p in np.linspace(0, 2 * np.pi, 40, endpoint=False):
        th = k * p; c = 1.09 * np.array([np.cos(p), np.sin(p)]); d = 0.065 * np.array([np.cos(th), np.sin(th)])
        ax.plot([c[0] - d[0], c[0] + d[0]], [c[1] - d[1], c[1] + d[1]], color="0.55", lw=0.8)
    N = 2 * k - 2; hp = 2 * np.pi * np.arange(N) / N
    ax.scatter(1.09 * np.cos(hp), 1.09 * np.sin(hp), s=16, facecolor="none", edgecolor=DARK, lw=0.8, zorder=5)
    pj = hp + tau / N; r = 0.8                                           # registry shift by tau/(2k-2) -> every tail tilted by tau
    for p in pj:
        x, y = r * np.cos(p), r * np.sin(p); a = p + tau                 # tail direction = normal + tau
        ax.plot([x, x + 0.2 * np.cos(a)], [y, y + 0.2 * np.sin(a)], color=RED, lw=1.3); ax.plot(x, y, "o", color=RED, ms=2.8, mec="w", mew=0.3)
    ax.set_xlim(-1.3, 1.3); ax.set_ylim(-1.3, 1.3); ax.set_aspect("equal"); ax.axis("off")
    ax.text(0, 0.12, f"$k={k}$: {N} sites", ha="center", fontsize=6.5); ax.text(0, -0.12, "$\\tau_j=(2k-2)\\phi_j+c$", ha="center", fontsize=6.5)


def snapshot(ax, run=CAV + "/server_main/main_act0.07/rho0.06/k06/s05", tt=9000.0):
    bc = json.load(open(run + "/input.json"))["BC"][0]; Q = np.array(bc["Q"][:2], float); R = bc["R"]
    rows, ft, idx = load_flow(run); nx, ny = rows[:, 1].astype(int).max() + 1, rows[:, 2].astype(int).max() + 1
    xs, ys = np.arange(nx) + 0.5 - Q[0], np.arange(ny) + 0.5 - Q[1]; XX, YY = np.meshgrid(xs, ys, indexing="ij"); inside = np.hypot(XX, YY) < R
    ii = idx[ft[np.argmin(np.abs(ft - tt))]]; U = np.zeros((nx, ny)); V = np.zeros((nx, ny))
    U[rows[ii, 1].astype(int), rows[ii, 2].astype(int)] = rows[ii, 4]; V[rows[ii, 1].astype(int), rows[ii, 2].astype(int)] = rows[ii, 5]
    w8 = gaussian_filter(inside.astype(float), 2.0) + 1e-9; Us = gaussian_filter(np.where(inside, U, 0), 2.0) / w8; Vs = gaussian_filter(np.where(inside, V, 0), 2.0) / w8
    om = np.gradient(Vs, axis=0) - np.gradient(Us, axis=1); om[~inside] = np.nan; vm = np.nanpercentile(np.abs(om), 97)
    pm = ax.pcolormesh(xs, ys, om.T, cmap="PuOr_r", vmin=-vm, vmax=vm, shading="nearest", rasterized=True)
    m = inside[2::4, 2::4] & (np.hypot(XX, YY)[2::4, 2::4] < R - 2)
    ax.quiver(XX[2::4, 2::4][m], YY[2::4, 2::4][m], Us[2::4, 2::4][m], Vs[2::4, 2::4][m], color="0.25", scale=14, width=0.004, headwidth=3.5)
    ph = np.linspace(0, 2 * np.pi, 300); ax.plot(R * np.cos(ph), R * np.sin(ph), "k-", lw=1.0)
    t, fr = read_with_angle(run, Q, 0.0); f = fr[int(np.argmin(np.abs(np.array(t) - tt)))]
    for x, y, q, a in f:
        if q > 0:
            h = a + np.pi; ax.plot([x - 3 * np.cos(h), x], [y - 3 * np.sin(h), y], color=RED, lw=1.6); ax.plot(x, y, "o", ms=3.2, color=RED, mec="w", mew=0.4)
        else:
            for j in range(3): b = a + 2 * np.pi * j / 3; ax.plot([x, x + 2.0 * np.cos(b)], [y, y + 2.0 * np.sin(b)], color=BLUE, lw=1.3)
    ax.set_xlim(-R - 1, R + 1); ax.set_ylim(-R - 1, R + 1); ax.set_aspect("equal"); ax.axis("off")
    return pm, sum(1 for x, y, q, a in f if q > 0 and R - np.hypot(x, y) < 4)


if __name__ == "__main__":
    fig = plt.figure(figsize=(PAGE, 4.65))
    top = fig.add_gridspec(1, 6, left=0.0, right=1.0, top=0.97, bottom=0.60, width_ratios=[1, 1, 1, 1, 0.12, 1.15], wspace=0.08)
    axm = [fig.add_subplot(top[0, i]) for i in range(4)]; ims, vmax = maps(axm); label(axm[0], "a", x=0.05, y=1.13)
    cax = fig.add_axes([0.07, 0.585, 0.45, 0.013]); cb = fig.colorbar(ims[-1], cax=cax, orientation="horizontal")
    cb.set_label("Time-averaged $+\\frac{1}{2}$ density (defects per cell$^2$ per frame); $-\\frac{1}{2}$ in orange, same scale", fontsize=6.5, labelpad=1); cb.ax.tick_params(labelsize=6)
    axl = fig.add_subplot(top[0, 5]); lock(axl); label(axl, "b", x=0.05, y=1.13)
    bot = fig.add_gridspec(2, 2, left=0.08, right=1.0, top=0.47, bottom=0.07, width_ratios=[1.55, 1], hspace=0.12, wspace=0.12)
    S = sweep(); cols = plt.cm.viridis(np.linspace(0, 0.92, 9))
    aq = fig.add_subplot(bot[0, 0]); aw = fig.add_subplot(bot[1, 0], sharex=aq)
    for k, c in zip(range(2, 11), cols):
        s = S[S[:, 0] == k]; s = s[np.argsort(s[:, 1])]
        aq.errorbar(s[:, 1], s[:, 2], s[:, 3], fmt="o-", ms=2.3, lw=0.8, capsize=1, color=c, label=f"{k}")
        aw.errorbar(s[:, 1], s[:, 4], s[:, 5], fmt="o-", ms=2.3, lw=0.8, capsize=1, color=c)
    for a_ in (aq, aw): a_.axvspan(0.8, 2.2, color="0.92", lw=0, zorder=0)
    aq.axhline(1, color="0.4", lw=0.6, ls="--"); aw.axhline(0, color="0.4", lw=0.6, ls="--")
    aq.set_ylim(0.6, 2.2); aq.set_ylabel("Bulk charge $Q_b$", labelpad=2); plt.setp(aq.get_xticklabels(), visible=False)
    aw.set_ylim(-3, 1.5); aw.set_ylabel("Wall $+\\frac{1}{2}$ $-(2k-2)$", labelpad=2); aw.set_xlabel("Cavity radius $R/\\ell_a$"); aw.set_xlim(0.2, 2.95)
    aq.legend(title="$k$", ncol=9, fontsize=5.8, title_fontsize=6, loc="upper center", bbox_to_anchor=(0.5, 1.02), handlelength=1.0, columnspacing=0.8, handletextpad=0.3)
    label(aq, "c", x=-0.1, y=1.0)
    ad = fig.add_subplot(bot[:, 1]); pm, nw = snapshot(ad); label(ad, "d", x=0.05, y=1.0)
    cax2 = fig.add_axes([0.66, 0.05, 0.3, 0.012]); cb2 = fig.colorbar(pm, cax=cax2, orientation="horizontal"); cb2.set_label("Flow vorticity (1/t.u.)", fontsize=6.5, labelpad=1); cb2.ax.tick_params(labelsize=6)
    print("snapshot wall +1/2:", nw); save(fig, "fig1_observation")
    win = S[(S[:, 1] >= 0.8) & (S[:, 1] <= 2.2)]
    print(f"window 0.8-2.2: Q_b mean {np.nanmean(win[:, 2]):.3f} (range {np.nanmin(win[:, 2]):.2f}-{np.nanmax(win[:, 2]):.2f}); wall excess mean {np.nanmean(win[:, 4]):+.3f} (range {np.nanmin(win[:, 4]):+.2f}..{np.nanmax(win[:, 4]):+.2f}); points {len(win)}")
