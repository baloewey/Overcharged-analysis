#!/usr/bin/env python3
"""Layer vs bulk +1/2 in winding channels (user, 2026-10-13): by distance from the wall, number per frame, median displacement over 10 t.u., and distance to the nearest
homeotropic anchoring point (in site spacings L/2m; 0 = on site, 0.25 = random). runs/channel W 40 and W 60 (m 6), second half, every 3rd frame. -> fig_layer_profile.png, layer_profile.txt"""
import sys, json, glob, os, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, "/home/baloewe/Physics/Overcharged/runs/cavity")
from edge_criteria import read_with_angle  # noqa
R = "/home/baloewe/Physics/Overcharged/runs/channel"; H = os.path.dirname(os.path.abspath(__file__)); E = np.arange(0, 20.5, 1.0); out = {}
for W in (40, 60):
    cnt = np.zeros(len(E) - 1); mob = [[] for _ in E[:-1]]; site = [[] for _ in E[:-1]]; nf = 0
    for run in sorted(glob.glob(f"{R}/W{W}_m6_L192_s*")):
        J = json.load(open(run + "/input.json")); y0, y1 = J["BC"][2]["Q"][1], J["BC"][3]["Q"][1]; L = J["domain"][0]; sp = L / (2 * J["BC"][2]["customAnch"])
        t, fr = read_with_angle(run, np.array([0.0, 0.0]), 0.5); lag = int(round(10 / np.median(np.diff(t))))
        for n in range(0, len(fr) - lag, 3):
            f = fr[n]
            if not len(f): continue
            nf += 1; d = np.minimum(f[:, 1] - y0, y1 - f[:, 1]); pl = f[:, 2] > 0; cnt += np.histogram(d[pl], E)[0]
            g = fr[n + lag]; gp = g[g[:, 2] > 0, :2] if len(g) else np.zeros((0, 2))
            for i in np.flatnonzero(pl):
                b = np.searchsorted(E, d[i], side="right") - 1
                if b < 0 or b >= len(E) - 1: continue
                ph = (f[i, 0] % sp) / sp; site[b].append(min(ph, 1 - ph))
                if len(gp):
                    dx = (gp[:, 0] - f[i, 0] + L / 2) % L - L / 2; dd = np.hypot(dx, gp[:, 1] - f[i, 1]); kk = np.argmin(dd)
                    if dd[kk] < 5: mob[b].append(dd[kk])
    out[W] = (cnt / nf, [np.median(x) if x else np.nan for x in mob], [np.mean(x) if x else np.nan for x in site])
xc = 0.5 * (E[1:] + E[:-1]); fig, ax = plt.subplots(1, 3, figsize=(13, 3.8)); col = {40: "#1f5fa8", 60: "#c2452d"}; L_ = []
for W, (c, m, s) in out.items():
    lim = xc < W / 2
    ax[0].plot(xc[lim], c[lim], "-o", ms=3, color=col[W], label=f"W = {W}"); ax[1].plot(xc[lim], np.array(m)[lim], "-o", ms=3, color=col[W]); ax[2].plot(xc[lim], np.array(s)[lim], "-o", ms=3, color=col[W])
    L_ += [f"W {W}: distance | +1/2 per frame | median displacement 10 t.u. | distance to nearest site"] + [f"  {E[b]:.0f}-{E[b+1]:.0f} | {c[b]:.2f} | {m[b]:.2f} | {s[b]:.3f}" for b in range(len(c)) if lim[b]]
for a_, yl, ti in zip(ax, ("+½ per frame (1-cell bins)", "median displacement in 10 t.u. (cells)", "distance to nearest site (site spacings)"),
                      ("(a) depth distribution", "(b) mobility", "(c) locking to the anchoring sites")):
    a_.set_xlabel("distance from the wall (cells)"); a_.set_ylabel(yl); a_.set_title(ti, fontsize=10); a_.axvspan(0, 4, color="0.92", zorder=0); [a_.spines[s_].set_visible(False) for s_ in ("top", "right")]
ax[2].axhline(0.25, color="0.6", ls=":", lw=1); ax[0].legend(frameon=False, fontsize=8); ax[0].text(2, 3.2, "wall band\n(< 4 cells)", ha="center", fontsize=7, color="0.4")
fig.tight_layout(); fig.savefig(os.path.join(H, "fig_layer_profile.png"), dpi=200); open(os.path.join(H, "layer_profile.txt"), "w").write("\n".join(L_) + "\n"); print("\n".join(L_[:8]))
