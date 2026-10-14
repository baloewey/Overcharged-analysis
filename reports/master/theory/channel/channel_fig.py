#!/usr/bin/env python3
"""Figure for the theory report: winding channels as the flat limit of the cavities (2026-10-13).
(a) passive and (b) active channel snapshots (L 192/96, W 40, director + defects); (c) net charge across the width for W 20/40/60 (active, q 0.031);
(d) bulk share of the imposed charge vs width, with the aspect-2 ellipse (k 6) and the q-matched elongated-cavity value. -> fig_channel.png"""
import sys, json, glob, os, gzip, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, "/home/baloewe/Physics/Overcharged/runs/cavity")
from edge_criteria import read_with_angle  # noqa
R = "/home/baloewe/Physics/Overcharged/runs"; H = os.path.dirname(os.path.abspath(__file__))
def opn(fn): return gzip.open(fn + ".gz", "rt") if os.path.exists(fn + ".gz") else open(fn)
def snapshot(ax, run, title, xmax=None):
    J = json.load(open(run + "/input.json")); L, Hd = J["domain"][:2]; y0, y1 = J["BC"][2]["Q"][1], J["BC"][3]["Q"][1]
    D = np.loadtxt((l for l in opn(run + "/directorfield.dat") if len(l.split()) >= 7 and l.split()[0][0].isdigit())); tl = D[:, 0].max(); D = D[(D[:, 0] == tl) & (D[:, 2] >= y0) & (D[:, 2] < y1)]
    if xmax: D = D[D[:, 1] < xmax]
    D = D[(D[:, 1] % 2 == 0) & (D[:, 2] % 2 == 0)]
    ax.quiver(D[:, 1] + 0.5, D[:, 2] + 0.5, D[:, 4], D[:, 5], headwidth=0, headlength=0, headaxislength=0, pivot="mid", scale=55, width=0.002, color="0.45")
    t, fr = read_with_angle(run, np.array([0.0, 0.0]), 0.0); f = fr[-1]
    for x, y, q, a in f:
        if xmax and x > xmax: continue
        if q > 0: ax.plot(x, y, "o", color="#c2452d", ms=5, zorder=5); ax.annotate("", (x - 2.5 * np.cos(a), y - 2.5 * np.sin(a)), (x, y), arrowprops=dict(arrowstyle="-|>", color="#c2452d", lw=1.2), zorder=5)
        else: ax.plot(x, y, "^", color="#1f5fa8", ms=5, zorder=5)
    ax.axhline(y0, color="k", lw=2.5); ax.axhline(y1, color="k", lw=2.5); ax.set_xlim(0, xmax or L); ax.set_ylim(y0 - 1, y1 + 1); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([]); ax.set_title(title, fontsize=9)
fig = plt.figure(figsize=(13, 8.4)); gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1.35], hspace=0.42, wspace=0.32)
snapshot(fig.add_subplot(gs[0, :]), f"{R}/channel_test/passive_margin", "(a) Passive ($L=96$, $W=40$): all $+1/2$ at the walls")
snapshot(fig.add_subplot(gs[1, :]), f"{R}/channel/W40_m6_L192_s1", "(b) Active ($L=192$, $W=40$, $q=0.031$): $+1/2$ released into the bulk", xmax=192)
ax = fig.add_subplot(gs[2, 0]); col = {20: "#1f5fa8", 40: "#7a3b8f", 60: "#c2452d"}; share = {}
for W in (20, 40, 60):
    runs = sorted(glob.glob(f"{R}/channel/W{W}_m6_L192_s*")); nb = 10; e = np.linspace(0, 1, nb + 1); hp = np.zeros(nb); hm = np.zeros(nb); nf = 0; qb = []
    for run in runs:
        J = json.load(open(run + "/input.json")); y0, y1 = J["BC"][2]["Q"][1], J["BC"][3]["Q"][1]
        t, fr = read_with_angle(run, np.array([0.0, 0.0]), 0.5)
        for f in fr[::5]:
            if not len(f): continue
            nf += 1; u = (f[:, 1] - y0) / (y1 - y0); hp += np.histogram(u[f[:, 2] > 0], e)[0]; hm += np.histogram(u[f[:, 2] < 0], e)[0]
            d = np.minimum(f[:, 1] - y0, y1 - f[:, 1]); qb.append(0.5 * np.sign(f[d >= 4, 2]).sum())
    net = 0.5 * (hp - hm) / nf; share[W] = np.mean(qb) / 12
    ax.plot(0.5 * (e[1:] + e[:-1]), net / (W / 10 * 192) * 1e3, "-o", ms=3, color=col[W], label=f"W = {W}")
ax.set_xlabel("position across the channel, y/W"); ax.set_ylabel("net charge (10⁻³ per cell²)"); ax.legend(frameon=False, fontsize=8); ax.set_title("(c) where the charge sits", fontsize=9)
ax2 = fig.add_subplot(gs[2, 1]); Ws = sorted(share); ax2.plot(Ws, [100 * share[w] for w in Ws], "o-", color="k", label="channel (q = 0.031)")
ax2.axhline(40, color="#c2452d", ls="--", lw=1, label="aspect-2 ellipse, k = 6 (40 %)"); ax2.axhline(41, color="0.6", ls=":", lw=1, label="elongated cavities at q ≈ 0.03 (41 %)")
ax2.set_ylim(0, 60); ax2.set_xlabel("channel width W (cells)"); ax2.set_ylabel("bulk share of the charge (%)"); ax2.legend(frameon=False, fontsize=8); ax2.set_title("(d) the same share at every width", fontsize=9)
for a_ in (ax, ax2): [a_.spines[s_].set_visible(False) for s_ in ("top", "right")]
fig.savefig(os.path.join(H, "fig_channel.png"), dpi=180, bbox_inches="tight"); print({w: round(100 * v, 1) for w, v in share.items()})
