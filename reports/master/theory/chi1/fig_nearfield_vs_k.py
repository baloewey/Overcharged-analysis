#!/usr/bin/env python3
"""Depth of the wall defects' near-field zone vs k, against the fixed band window 2-3 cells (user, 2026-11-24).
Near-field edge = d + a/2pi (d = wall-defect depth, a = spacing 2 pi r/(2k-2)). Two estimates: depth law of the elastic ring with a central +1
(rs = R((N+3)/(3N+3))^(1/2N) - 0.61, as tilt_energy_landscape/ring) and the MEASURED MPCD depth and spacing (tiltcrowd/defect_flow_main10.txt, active circles,
10 seeds; R 31.83 only). R 31.83 (chi 0.5, R/l_a 1.43) and R 21.6 (chi 1 at the same R/l_a; depth law only).
-> fig_nearfield_vs_k.pdf/.png, fig_nearfield_vs_k.txt"""
import os, re, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
H = os.path.dirname(os.path.abspath(__file__)); TC = os.path.join(H, "..", "tiltcrowd")


def law(R, k):
    N = 2 * k - 2; rs = R * ((N - 1 + 4) / (3 * N - 1 + 4)) ** (1 / (2 * N)) - 0.61; return R - rs, 2 * np.pi * rs / N


meas = {}
for line in open(os.path.join(TC, "defect_flow_main10.txt")):
    m = re.match(r"rho([\d.]+)/k(\d+)/s(\d+)\s*\|", line)
    if not m: continue
    k = int(m.group(2)); R = k / (np.pi * float(m.group(1)))
    if abs(R - 31.83) > 0.6: continue
    f = [x.strip() for x in line.split("|")]; meas.setdefault(k, []).append((float(f[3]), float(f[2])))
ks = np.arange(3, 21)
fig, ax = plt.subplots(figsize=(9.5, 7.2)); L = ["k | R 31.83 law: d, a, d + a/2pi | measured: d, a, d + a/2pi | R 21.6 law: d + a/2pi"]
ax.axhspan(2, 3, color="#1f5fa8", alpha=0.15, lw=0, label="Band window (2–3 cells)")
for R, c, lab in ((31.83, "#c0392b", r"$R$ = 31.8 ($\chi$ = 0.5)"), (21.6, "#7d3c98", r"$R$ = 21.6 ($\chi$ = 1)")):
    d = np.array([law(R, k)[0] for k in ks]); a = np.array([law(R, k)[1] for k in ks]); edge = d + a / (2 * np.pi)
    ax.plot(ks, edge, "-", color=c, lw=2.5, label=f"Near-field edge d + a/2π, {lab}")
    ax.plot(ks, d, ":", color=c, lw=2, label=f"Defect depth d, {lab}")
    cross = ks[np.argmax(edge < 3)] if (edge < 3).any() else None; cross2 = ks[np.argmax(edge < 2)] if (edge < 2).any() else None
    if cross: ax.axvline(cross, color=c, ls="--", lw=1.2); ax.text(cross + 0.15, 10.6 if R > 30 else 9.0, f"edge < 3 cells\nfrom k = {cross}", color=c, fontsize=11, va="top")
    L.append(f"R {R}: edge drops below 3 cells at k = {cross}, below 2 cells at k = {cross2}; defect depth d below 2 cells from k = {ks[np.argmax(d < 2)] if (d < 2).any() else None}")
km = sorted(meas); dm = np.array([np.mean([x[0] for x in meas[k]]) for k in km]); am = np.array([np.mean([x[1] for x in meas[k]]) for k in km])
ax.plot(km, dm + am / (2 * np.pi), "o", color="#c0392b", ms=8, mec="white", label="Measured MPCD (R 31.8): d + a/2π")
ax.plot(km, dm, "s", color="#c0392b", mfc="white", mew=1.8, ms=7, label="Measured MPCD (R 31.8): d")
for k, dd, aa in zip(km, dm, am): L.append(f"  measured k {k:2d}: d {dd:.2f}, a {aa:.1f}, d + a/2pi {dd + aa / (2 * np.pi):.2f}")
ax.set_xlabel("Winding k"); ax.set_ylabel("Depth below the wall (cells)"); ax.set_xlim(2.5, 20.5); ax.set_ylim(0, 11)
ax.set_title("Near-field zone of the wall defects vs the fixed band window", fontsize=13); ax.legend(fontsize=9.5, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=2)
fig.tight_layout(); fig.savefig(os.path.join(H, "fig_nearfield_vs_k.pdf")); fig.savefig(os.path.join(H, "fig_nearfield_vs_k.png"), dpi=200)
open(os.path.join(H, "fig_nearfield_vs_k.txt"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
