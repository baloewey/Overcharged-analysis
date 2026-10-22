#!/usr/bin/env python3
"""Defect heights from birth to death (user, 2026-10-22: 'defect height as a function of time for each defect'), periodic channels L 192, W 24, 6000 t.u.
(runs/channel_current W24_mixed_s3: bottom +6 (+1/2 layer), top -6 (-1/2 layer); W24_sym_s3: +6/+6). Tracks from tracker2.build (predictive Hungarian,
gate 9 cells, 2-frame gaps). (a,b) height y of every track against time in a 400 t.u. window (red +1/2, blue -1/2; dashed: mid-line; grey bands: 4-cell wall
bands). (c) charge carried up across the mid-line: tracking-free Q_up(t) - Q_up(0) (black), the new tracked count (crossings + pairs straddling the line;
orange) and the old nearest-neighbour count (2.5-cell links, no prediction; grey) for the mixed channel. (d) the decomposition over the whole run: +1/2 and
-1/2 crossings up and down (counts) and the net charge each carries, plus pairs born or annihilated across the line."""
import sys, numpy as np, matplotlib.pyplot as plt
from tracker2 import build
R = "/home/baloewe/Physics/Overcharged/runs/channel_current/"
M = build(R + "W24_mixed_s3", gate=9.0); S = build(R + "W24_sym_s3", gate=9.0); O = build(R + "W24_mixed_s3", gate=2.5, predict=False)
fig = plt.figure(figsize=(17, 10)); g = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1.1], hspace=0.45, wspace=0.22, left=0.06, right=0.99, top=0.95, bottom=0.07)
T0, T1 = 3000, 3400
for i, (B, ttl) in enumerate(((M, "Mixed walls: bottom +6 (+1/2 layer), top −6 (−1/2 layer)"), (S, "Symmetric control: +6 / +6 (both +1/2 layers)"))):
    ax = fig.add_subplot(g[i, :]); t = B["t"]
    for tr in B["tracks"]:
        tt = t[np.array(tr["n"])]; m = (tt >= T0) & (tt <= T1)
        if m.sum() < 2: continue
        ax.plot(tt[m], np.array(tr["xy"])[m, 1], "-", color="#c0392b" if tr["sg"] > 0 else "#2e6fb5", lw=0.7, alpha=0.8)
    ax.axhline(B["ym"], color="0.3", ls="--", lw=0.9); ax.axhspan(B["y0"], B["y0"] + 4, color="0.9", zorder=0); ax.axhspan(B["y1"] - 4, B["y1"], color="0.9", zorder=0)
    ax.set_xlim(T0, T1); ax.set_ylim(B["y0"], B["y1"]); ax.set_xlabel("Time (t.u.)"); ax.set_ylabel("Height $y$ (cells)"); ax.set_title(ttl + "  —  one line per defect, birth to death", fontsize=11, loc="left")
    ax.text(-0.045, 1.04, "(" + "ab"[i] + ")", transform=ax.transAxes, fontsize=14, fontweight="bold")
ax = fig.add_subplot(g[2, 0]); t = M["t"]; q = M["Qup"] - np.nanmean(M["Qup"][:25])
ax.plot(t, O["cum"], "-", color="0.6", lw=1.3, label="Old count (nearest neighbour, 2.5-cell links)"); ax.plot(t, M["cum"], "-", color="#e6550d", lw=1.3, label="New tracked count (crossings + straddling pairs)")
ax.plot(t, q, "-", color="k", lw=1.0, label="Tracking-free: $Q_\\mathrm{up}(t)-Q_\\mathrm{up}(0)$")
ax.set_xlabel("Time (t.u.)"); ax.set_ylabel("Charge carried up"); ax.legend(fontsize=9, loc="upper left"); ax.set_title("Mixed channel: net charge across the mid-line", fontsize=11, loc="left")
ax.text(-0.1, 1.04, "(c)", transform=ax.transAxes, fontsize=14, fontweight="bold")
ax = fig.add_subplot(g[2, 1]); labs = ["$+\\frac{1}{2}$ up", "$+\\frac{1}{2}$ down", "$-\\frac{1}{2}$ up", "$-\\frac{1}{2}$ down"]; x = np.arange(4)
for j, (B, c, nm) in enumerate(((M, "#e6550d", "Mixed"), (S, "#3b6fb6", "Symmetric"))):
    P = B["parts"]; ax.bar(x + (j - 0.5) * 0.38, [P["p_up"], P["p_down"], P["m_up"], P["m_down"]], 0.36, color=c, label=nm)
    net = 0.5 * (P["p_up"] - P["p_down"]) - 0.5 * (P["m_up"] - P["m_down"])
    ax.text(0.98, 0.93 - 0.12 * j, f"{nm}: net by $+\\frac{{1}}{{2}}$ {0.5 * (P['p_up'] - P['p_down']):+.1f}, by $-\\frac{{1}}{{2}}$ {-0.5 * (P['m_up'] - P['m_down']):+.1f}, "
            f"straddling pairs {P['pair_birth'] + P['pair_death']:+.1f}; total {net + P['pair_birth'] + P['pair_death']:+.1f}", transform=ax.transAxes, ha="right", fontsize=8.5, color=c)
ax.set_xticks(x); ax.set_xticklabels(labs); ax.set_ylabel("Mid-line crossings in 6000 t.u."); ax.legend(fontsize=9, loc="center right", bbox_to_anchor=(1.0, 0.55)); ax.set_ylim(0, 1050)
ax.set_title("Traffic across the mid-line (new tracker)", fontsize=11, loc="left"); ax.text(-0.1, 1.04, "(d)", transform=ax.transAxes, fontsize=14, fontweight="bold")
for e in ("pdf", "png"): fig.savefig(f"fig_heights.{e}", dpi=170)
