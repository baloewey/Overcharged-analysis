#!/usr/bin/env python3
"""New graphics for the workshop talk (27 Nov 2026). Run with reports/figstyle/run.sh talk_figs.py.
  assets/anchoring.png  winding anchoring theta = k phi on a circle (k = 4): director sticks at the wall, 2k-2 homeotropic sites
  assets/rings.png      the lock tau_j = (2k-2) phi_j + c: a ring of 2k-2 tilts as a whole by sliding; a ring of 2k cannot (tilts wind twice)
  assets/coreflow.png   flow-at-core result (MPCD runs/coreflow, aspect-2 ellipse, k 5 vs 14, 5 seeds; depth 1-2 cells): inward flow at the core,
                        response of the core to it (OLS slope v_n on u_n), escape-frame fraction; mean +- s.e. over seeds (coreflow.one_mpcd per seed)"""
import os, sys, glob, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from multiprocessing import Pool
H = os.path.dirname(os.path.abspath(__file__)); A = os.path.join(H, "assets")
sys.path.insert(0, os.path.join(H, "..", "master", "theory", "coreflow"))
RED, GREY, NAVY = "#c0392b", "#555555", "#000014"
plt.rcParams.update({"font.size": 14, "axes.spines.top": False, "axes.spines.right": False})


def wall_sticks(ax, k, R=1.0, n=48, L=0.13):
    ph = np.linspace(0, 2 * np.pi, n, endpoint=False); th = k * ph
    for p, t in zip(ph, th):
        c = (R + 0.11) * np.array([np.cos(p), np.sin(p)]); d = 0.5 * L * np.array([np.cos(t), np.sin(t)])
        ax.plot(*np.c_[c - d, c + d], color=GREY, lw=2, solid_capstyle="round")
    ax.add_patch(plt.Circle((0, 0), R, fill=False, lw=2.5, color="#333333"))
    sites = np.arange(2 * k - 2) * np.pi / (k - 1)
    ax.plot((R + 0.27) * np.cos(sites), (R + 0.27) * np.sin(sites), "o", ms=9, mfc="white", mec=NAVY, mew=2)
    return sites


def defect(ax, p, tau, R=0.86, L=0.16):
    x, y = R * np.cos(p), R * np.sin(p); a = p + tau                       # tail direction: outward normal rotated by the tilt tau
    ax.plot([x, x + L * np.cos(a)], [y, y + L * np.sin(a)], color=RED, lw=4, solid_capstyle="round"); ax.plot(x, y, "o", color=RED, ms=9, mec="white", mew=1.5)


def anchoring():
    fig, ax = plt.subplots(figsize=(4.6, 4.6)); k = 4; wall_sticks(ax, k)
    ax.text(0, 0.08, f"Winding $k={k}$", ha="center", fontsize=13); ax.text(0, -0.2, r"Total charge $k$", ha="center", fontsize=11, color=GREY)
    ax.text(0, -1.62, "Open circles: $2k-2$ homeotropic sites", ha="center", fontsize=10, color=NAVY)
    ax.set_xlim(-1.5, 1.5); ax.set_ylim(-1.75, 1.45); ax.set_aspect("equal"); ax.axis("off")
    fig.savefig(os.path.join(A, "anchoring.png"), dpi=300, bbox_inches="tight", facecolor="white"); plt.close(fig)


def rings():
    fig, axs = plt.subplots(1, 2, figsize=(10, 5.2)); k = 4
    ax = axs[0]; sites = wall_sticks(ax, k); slide = 0.09                    # 2k-2 ring slid by the same arc: common tilt (2k-2)*slide
    for s in sites: defect(ax, s + slide, (2 * k - 2) * slide)
    ax.set_title(f"$2k-2={2 * k - 2}$: slides as a whole\n→ common tilt, no elastic cost", fontsize=16, color="#1a7a3a")
    ax = axs[1]; wall_sticks(ax, k); ph = np.arange(2 * k) * 2 * np.pi / (2 * k) + 0.09
    for p in ph: defect(ax, p, (2 * k - 2) * p)
    ax.set_title(f"$2k={2 * k}$: tilts wind twice\n→ no common tilt possible", fontsize=16, color=RED)
    for ax in axs: ax.set_xlim(-1.5, 1.5); ax.set_ylim(-1.5, 1.5); ax.set_aspect("equal"); ax.axis("off")
    fig.text(0.5, 0.02, r"Lock: $\tau_j=(2k-2)\,\phi_j+c$   (tilt of each wall $+\frac{1}{2}$ set by its position; $k=4$)", ha="center", fontsize=15)
    fig.savefig(os.path.join(A, "rings.png"), dpi=300, bbox_inches="tight", facecolor="white"); plt.close(fig)


def per_seed(run):
    import coreflow as C
    r = C.one_mpcd(run); m = (r[:, 6] >= 1) & (r[:, 6] < 2); x = r[m]
    return dict(k=int(os.path.basename(run).split("_k")[1][:2]), u=x[:, 1].mean(), slope=np.polyfit(x[:, 1], x[:, 2], 1)[0], esc=np.mean(x[:, 0] == 2), n=len(x))


def coreflow():
    runs = sorted(glob.glob("/home/baloewe/Physics/Overcharged/runs/coreflow/ell2_k*_s*"))
    with Pool(3) as p: S = p.map(per_seed, runs)
    cols = {5: "#f6a04d", 14: "#4a1c0c"}; lab = {5: "Sparse wall\n($\\rho=0.04$)", 14: "Crowded wall\n($\\rho=0.12$)"}
    panels = [("u", "Inward flow at the core\n(cells / t.u.)"), ("slope", "Response of the core\nto that flow"), ("esc", "Fraction of frames\njust before an escape")]
    fig, axs = plt.subplots(1, 3, figsize=(12, 4.4)); lines = ["COREFLOW TALK FIGURE (depth 1-2 cells; per seed)"]
    for ax, (key, title) in zip(axs, panels):
        for i, k in enumerate((5, 14)):
            v = np.array([s[key] for s in S if s["k"] == k]); m, e = v.mean(), v.std(ddof=1) / np.sqrt(len(v))
            ax.bar(i, m, 0.62, color=cols[k], yerr=e, capsize=6, error_kw=dict(lw=1.8)); lines.append(f"{key} k{k}: {m:.4f} +- {e:.4f} (seeds {np.round(v, 4).tolist()})")
            ax.text(i, m + e, f" {m:.3f}" if key != "esc" else f" {m:.4f}", ha="center", va="bottom", fontsize=14)
        ax.set_xticks([0, 1]); ax.set_xticklabels([lab[5], lab[14]], fontsize=13); ax.set_title(title, fontsize=16); ax.set_ylim(0, ax.get_ylim()[1] * 1.18)
    fig.text(0.5, -0.04, "MPCD, aspect-2 ellipse, wall $+\\frac{1}{2}$ at 1–2 cells from the wall; mean ± s.e. over 5 seeds (preliminary)", ha="center", fontsize=13, color=GREY)
    fig.tight_layout(); fig.savefig(os.path.join(A, "coreflow.png"), dpi=300, bbox_inches="tight", facecolor="white"); plt.close(fig)
    open(os.path.join(A, "coreflow_talk.txt"), "w").write("\n".join(lines) + "\n"); print("\n".join(lines))


if __name__ == "__main__":
    anchoring(); rings(); coreflow()
