#!/usr/bin/env python3
"""Figures: (1) fluct_mc.png: passive circles R = 31.83, measured tangential/normal rms vs exact-potential Monte Carlo (mc_exact.txt) with
K = K3 = 23 kBT, wall offset 0.61, positions rounded to the cell grid as in the detector; harmonic K = 32.5 for reference.
(2) active_wander.png: sigma_gap (= equivalent sigma_t of the stage-2 ensembles) and wall occupancy, active circles k = 6 vs R/l_a, and
active/passive ellipses."""
import os, re
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
H = os.path.dirname(os.path.abspath(__file__))
txt = open(os.path.join(H, "mc_exact.txt")).read().splitlines()
meas, mc, mch = {}, {}, {}
k = None
for l in txt:
    m = re.match(r"k=\s*(\d+) N=\s*(\d+) measured: depth ([\d.]+), tangential internal ([\d.]+), normal ([\d.]+)", l)
    if m: k = int(m.group(2)); meas[k] = (float(m.group(4)), float(m.group(5))); continue
    m = re.match(r"\s+K ([\d.]+), wall offset ([\d.]+), (grid-rounded|exact)\s*: depth ([\d.]+), tangential ([\d.]+), normal ([\d.]+)", l)
    if m and float(m.group(2)) > 0:
        key = (float(m.group(1)), m.group(3)); (mc.setdefault(key, {}))[k] = (float(m.group(5)), float(m.group(6)))
N = sorted(meas)
fig, ax = plt.subplots(1, 2, figsize=(9.5, 3.6))
for j, nm in enumerate(("tangential (internal)", "normal")):
    a = ax[j]
    a.plot(N, [meas[n][j] for n in N], "o", color="k", label="measured (passive)")
    a.plot(N, [mc[(23.0, "grid-rounded")][n][j] for n in N], "-", color="C0", label=r"exact potential, $K=K_3=23$, grid-rounded")
    a.plot(N, [mc[(23.0, "exact")][n][j] for n in N], "--", color="C0", lw=0.8, label=r"exact potential, $K=K_3$, exact positions")
    a.plot(N, [mc[(32.5, "grid-rounded")][n][j] for n in N], "-", color="C1", lw=0.8, label=r"exact potential, $K=32.5$, grid-rounded")
    a.set_xlabel("wall defects $N=2k$"); a.set_ylabel("rms (cells)"); a.set_title(f"({'ab'[j]}) {nm}, circles $R=31.83$", fontsize=10); a.set_ylim(0, None)
ax[0].legend(fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(H, "fluct_mc.png"), dpi=150)
# active wander
rows = []
for l in open(os.path.join(H, "active_wander.txt")):
    m = re.match(r"(\S.*?)\s+k=\s*6 wall count\s+([\d.]+)±([\d.]+) \(mode\s+(\d+), (\d+)% of frames\), spacing\s+([\d.]+); sigma_gap ([\d.]+)", l)
    if m: rows.append((m.group(1).strip(), float(m.group(2)), float(m.group(7))))
la = 22.2
fig, ax = plt.subplots(1, 2, figsize=(9.5, 3.6))
circ = [(float(re.search(r"R ([\d.]+)", r[0]).group(1)), r[1], r[2]) for r in rows if r[0].startswith("active circle")]
circ = np.array(circ)
ax[0].plot(circ[:, 0] / la, circ[:, 2], "o", label=r"active circles $k=6$ (3 seeds)")
pc = [r for r in rows if r[0].startswith("passive circle")]
ax[0].axhspan(min(r[2] for r in pc), max(r[2] for r in pc), color="gray", alpha=0.3, label="passive circles $k=6$")
ax[0].set_yscale("log"); ax[0].set_xlabel(r"$R/\ell_a$"); ax[0].set_ylabel(r"$\sigma_\mathrm{gap}=\mathrm{sd(gap)}/\sqrt{2}$ (cells)")
ax[0].set_title("(a) circles: wall-lattice disorder", fontsize=10); ax[0].legend(fontsize=7)
eps = {"1.2": 1.2, "1.5": 1.5, "2.0": 2.0}
for lab, mk, c in (("active ellipse", "o", "C3"), ("passive ellipse", "s", "gray")):
    pts = [(eps[re.search(r"eps([\d.]+)", r[0]).group(1)], r[2], r[1]) for r in rows if r[0].startswith(lab)]
    pts = np.array(pts); ax[1].plot(pts[:, 0], pts[:, 1], mk, color=c, label=lab + ("s (2 seeds)" if "active" in lab else "s"))
ax[1].set_yscale("log"); ax[1].set_xlabel(r"aspect ratio $\varepsilon$"); ax[1].set_ylabel(r"$\sigma_\mathrm{gap}$ (cells)")
ax[1].set_title("(b) ellipses, area 2500, $k=6$", fontsize=10); ax[1].legend(fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(H, "active_wander.png"), dpi=150)
print("ok")
