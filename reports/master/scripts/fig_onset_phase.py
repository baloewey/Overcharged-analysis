#!/usr/bin/env python3
"""
Circles (+k, k >= 2, rho <= 0.15, 3-seed means): onset of the layer/bulk structure and a phase diagram (master document).

Inputs (all per run, main sweep): runs/cavity/analysis_weekend/criteria_all.csv (q_slow = bulk charge from the
mobility criterion, q_dens = density-minimum criterion, q_knee), analysis_flow/flow_tilt.csv (vth_abs, vrms,
tilt_flow), analysis_fourier/wall_fourier_all.csv (fstatic). Means over seeds per (rho, k).
Onset: for each rho, the R/l_a at which an observable first crosses a threshold (linear interpolation between
neighbouring sizes): circulation fraction |<v_theta>|/v_rms = 0.5, tilt signed with the flow = 30 deg, q_slow = 0.5,
f_static = 0.75, q_dens = 0.5.
Phases (per (rho, k)): no bulk q_slow < 0.5; forming 0.5 <= q_slow < 0.75; +1 bulk + pinned layer 0.75 <= q_slow <= 1.25
and f_static >= 0.75; less-charged pinned layer q_slow > 1.25 and f_static >= 0.75; no pinned layer f_static < 0.75
(q_slow >= 0.75).
Writes ../figs/onset_circles.png, ../figs/phase_circles.png, ../data/onset_circles.tsv.
"""
import csv, os
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, ListedColormap
from matplotlib.patches import Patch

HERE = os.path.dirname(os.path.abspath(__file__))
CAV = os.path.abspath(os.path.join(HERE, "..", "..", "..", "runs", "cavity"))
FIGS = os.path.join(HERE, "..", "figs"); DATA = os.path.join(HERE, "..", "data")
LA = 22.21
INK, INK2, GRID = "#15171b", "#555a62", "#e9ebee"


def load(path, cols):
    D = defaultdict(lambda: defaultdict(list))
    for r in csv.DictReader(open(path)):
        if r["set"] != "pos": continue
        key = (round(float(r["rho"]), 3), int(r["k"]))
        D[key]["R"] = [float(r["R"])]
        for c in cols:
            try: D[key][c].append(float(r[c]))
            except (ValueError, KeyError): pass
    return D


C = load(os.path.join(CAV, "analysis_weekend", "criteria_all.csv"), ["q_slow", "q_dens", "q_knee"])
F = load(os.path.join(CAV, "analysis_flow", "flow_tilt.csv"), ["vth_abs", "vrms", "tilt_flow"])
W = load(os.path.join(CAV, "analysis_fourier", "wall_fourier_all.csv"), ["fstatic"])
mean = lambda d, c: float(np.nanmean(d[c])) if d.get(c) else np.nan
rows = []
for key in sorted(C):
    rho, k = key
    if rho > 0.151 or k < 2: continue          # k = 1: no wall lattice (2k - 2 = 0 sites), bulk charge trivially 1
    c, f, w = C[key], F.get(key, {}), W.get(key, {})
    circ = mean(f, "vth_abs") / mean(f, "vrms") if f and f.get("vrms") else np.nan
    rows.append(dict(rho=rho, k=k, x=c["R"][0] / LA, q_slow=mean(c, "q_slow"), q_dens=mean(c, "q_dens"), q_knee=mean(c, "q_knee"),
                     circ=circ, tilt=np.degrees(mean(f, "tilt_flow")) if f else np.nan, fstat=mean(w, "fstatic") if w else np.nan))
RHOS = sorted({r["rho"] for r in rows})
OBS = [("circ", "circulation fraction |⟨v_θ⟩|/v_rms", 0.5, (0, 1.05)),
       ("tilt", "tilt of wall +½ with the flow (deg)", 30, (-5, 75)),
       ("q_slow", "bulk charge, mobility criterion", 0.5, (-0.1, 2.6)),
       ("fstat", "wall lattice pinning f_static", 0.75, (0, 1.05)),
       ("q_dens", "bulk charge, density-minimum criterion", 0.5, (-0.1, 2.6))]


def cross(pts, thr):
    for (x0, y0), (x1, y1) in zip(pts, pts[1:]):
        if np.isfinite(y0) and np.isfinite(y1) and y0 < thr <= y1:
            return x0 + (thr - y0) * (x1 - x0) / (y1 - y0)
    return np.nan


def style(ax):
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"): ax.spines[sp].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=8); ax.grid(color=GRID, lw=0.6); ax.set_axisbelow(True)


# ---------------------------------------------------------------- onset figure
onset = {o[0]: {} for o in OBS}
fig, axs = plt.subplots(2, 3, figsize=(13.5, 8))
cmap = plt.cm.viridis
for ax, (key, lab, thr, yl) in zip(axs.flat, OBS):
    for i, rho in enumerate(RHOS):
        pts = sorted((r["x"], r[key]) for r in rows if r["rho"] == rho)
        ax.plot([p[0] for p in pts], [p[1] for p in pts], "-o", ms=2.5, lw=1.1, color=cmap(i / (len(RHOS) - 1)), label=f"ρ = {rho:.2f}")
        onset[key][rho] = cross(pts, thr)
    xs = np.array([v for v in onset[key].values() if np.isfinite(v)])
    med = np.median(xs); q1, q3 = np.percentile(xs, [25, 75])
    ax.axhline(thr, color=INK2, lw=0.8, ls=":")
    ax.axvspan(q1, q3, color="#2a78d6", alpha=0.12, lw=0); ax.axvline(med, color="#2a78d6", lw=1.2)
    ax.text(med + 0.04, yl[1] - 0.06 * (yl[1] - yl[0]), f"onset {med:.2f} ℓa\n(IQR {q1:.2f}–{q3:.2f})", fontsize=7.5, color="#1c5cab", va="top")
    ax.set_xlim(0, 3); ax.set_ylim(*yl); style(ax)
    ax.set_xlabel("R / ℓa", fontsize=9, color=INK); ax.set_title(lab, fontsize=9.5, color=INK, loc="left")
axs.flat[0].legend(frameon=False, fontsize=6.5, ncol=2, loc="lower right")
ax = axs.flat[5]
names = {"circ": "circulation", "tilt": "tilt", "q_slow": "bulk charge (mobility)", "fstat": "pinning", "q_dens": "bulk charge (density gap)"}
order = sorted(onset, key=lambda k: np.nanmedian(list(onset[k].values())))
for j, key in enumerate(order):
    xs = np.array([v for v in onset[key].values() if np.isfinite(v)])
    ax.scatter(xs, np.full(len(xs), j) + np.linspace(-0.18, 0.18, len(xs)), s=14, color="#7c818a", zorder=2)
    med = np.median(xs); q1, q3 = np.percentile(xs, [25, 75])
    ax.plot([q1, q3], [j, j], color="#2a78d6", lw=5, alpha=0.5, solid_capstyle="round"); ax.plot(med, j, "o", color="#1c5cab", ms=8, zorder=3)
ax.set_yticks(range(len(order))); ax.set_yticklabels([names[k] for k in order], fontsize=8.5, color=INK)
ax.set_xlim(0.2, 1.3); ax.set_xlabel("onset R / ℓa (one dot per ρ; bar: IQR; big dot: median)", fontsize=9, color=INK)
ax.set_title("sequence of onsets", fontsize=9.5, color=INK, loc="left"); style(ax)
fig.suptitle("Circles, +k, k ≥ 2, ρ ≤ 0.15 (3-seed means): onset of circulation, tilt, bulk charge and pinning with cavity size. "
             "Dotted: threshold; blue line and band: median and IQR over ρ of the crossing size", fontsize=9, color=INK2)
fig.tight_layout(); fig.savefig(os.path.join(FIGS, "onset_circles.png"), dpi=170); plt.close(fig)
os.makedirs(DATA, exist_ok=True)
with open(os.path.join(DATA, "onset_circles.tsv"), "w") as f:
    f.write("rho\t" + "\t".join(o[0] for o in OBS) + "\n")
    for rho in RHOS: f.write(f"{rho:.2f}\t" + "\t".join(f"{onset[o[0]][rho]:.3f}" for o in OBS) + "\n")

# ---------------------------------------------------------------- phase diagram
PH = ["no bulk", "forming", "+1 bulk, pinned layer", "less-charged pinned layer", "no pinned layer"]
PCOL = ["#cfd3d8", "#9ec5f4", "#2a78d6", "#eb6834", "#eda100"]


def phase(r):
    q, fs = r["q_slow"], r["fstat"]
    if not np.isfinite(q): return None
    if q < 0.5: return 0
    if q < 0.75: return 1
    if np.isfinite(fs) and fs < 0.75: return 4
    return 2 if q <= 1.25 else 3


def cells(ax, key, cmap=None, norm=None, cat=False):
    for rho in RHOS:
        sub = sorted((r for r in rows if r["rho"] == rho), key=lambda r: r["x"])
        xs = [r["x"] for r in sub]
        edges = [0] + [0.5 * (a + b) for a, b in zip(xs, xs[1:])] + [xs[-1] + 0.5 * (xs[-1] - xs[-2]) if len(xs) > 1 else xs[-1] + 0.1]
        i = RHOS.index(rho)
        lo = rho - 0.5 * (rho - RHOS[i - 1]) if i else rho - 0.005
        hi = rho + 0.5 * (RHOS[i + 1] - rho) if i < len(RHOS) - 1 else rho + 0.005
        for r, e0, e1 in zip(sub, [max(edges[0], xs[0] - 0.5 * (edges[1] - xs[0]))] + edges[1:-1], edges[1:]):
            if cat:
                p = phase(r)
                if p is None: continue
                col = PCOL[p]
            else:
                v = r[key]
                if not np.isfinite(v): continue
                col = cmap(norm(v))
            ax.add_patch(plt.Rectangle((lo, e0), hi - lo, e1 - e0, facecolor=col, edgecolor="white", lw=0.4))
    ax.set_xlim(RHOS[0] - 0.006, RHOS[-1] + 0.006); ax.set_ylim(0, 3)
    ax.set_xlabel("ρ = 2k / P", fontsize=9, color=INK); ax.set_ylabel("R / ℓa", fontsize=9, color=INK); style(ax)


fig, axs = plt.subplots(1, 3, figsize=(15, 5.6), gridspec_kw=dict(width_ratios=[1.25, 1, 1]))
cells(axs[0], None, cat=True)
med = {k: np.nanmedian(list(onset[k].values())) for k in onset}
axs[0].axhline(med["circ"], color=INK, lw=1, ls="--"); axs[0].text(RHOS[0] - 0.004, med["circ"] + 0.03, f"circulation onset {med['circ']:.2f}", fontsize=7, va="bottom", color=INK)
axs[0].legend(handles=[Patch(facecolor=c, edgecolor="white", label=l) for c, l in zip(PCOL, PH)], frameon=False, fontsize=7.5,
              loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3)
axs[0].set_title("(a) phases: bulk charge (mobility) + pinning of the wall lattice", fontsize=9.5, color=INK, loc="left")
cm = plt.cm.RdBu_r; nrm = TwoSlopeNorm(vmin=0, vcenter=1, vmax=2.5)
for ax, key, ti in ((axs[1], "q_slow", "(b) bulk charge, mobility criterion"), (axs[2], "q_dens", "(c) bulk charge, density-minimum criterion")):
    cells(ax, key, cm, nrm); ax.set_title(ti, fontsize=9.5, color=INK, loc="left")
sm = plt.cm.ScalarMappable(cmap=cm, norm=nrm); fig.colorbar(sm, ax=axs[1:], fraction=0.03, pad=0.02, label="bulk charge (3-seed mean)")
fig.suptitle("Circles, +k, k ≥ 2, ρ ≤ 0.15. Phases: no bulk q < 0.5; forming 0.5–0.75; +1 bulk with pinned layer 0.75–1.25 and f_static ≥ 0.75; "
             "less-charged pinned layer q > 1.25, f_static ≥ 0.75; no pinned layer f_static < 0.75", fontsize=8.5, color=INK2)
fig.savefig(os.path.join(FIGS, "phase_circles.png"), dpi=170, bbox_inches="tight"); plt.close(fig)
counts = defaultdict(int)
for r in rows:
    p = phase(r)
    if p is not None: counts[PH[p]] += 1
print("onset medians:", {k: round(v, 2) for k, v in med.items()}); print("phase counts:", dict(counts))
