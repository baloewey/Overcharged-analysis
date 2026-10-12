#!/usr/bin/env python3
"""
Flow, enstrophy and wall-defect tilt in circles (main sweep, zeta = 0.07, rho <= 0.15), from flow_tilt.py.

  page 1  mean flow: |<v_theta>|, v_rms and the circulation fraction |<v_theta>| / v_rms vs R/l_a (one curve per rho;
          red +k, blue -k) and as (rho, R/l_a) maps
  page 2  enstrophy <omega^2/2> vs R/l_a and maps
  page 3  tilt of the wall-sign defects relative to the outward normal (polar angle), signed with the circulation
          sense: curves, maps and director-field diagrams of what the angle means
Writes analysis_flow/flow_maps.pdf and one PNG/PDF per page.
"""
import collections, csv, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.patches import Arc, FancyArrowPatch

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from plots_criteria import tiles, mean, LA, INK, INK2, GRID  # noqa: E402
OUT = os.path.join(HERE, "analysis_flow")
DEG = 180 / np.pi


def load():
    G = collections.defaultdict(list)
    for r in csv.DictReader(open(os.path.join(OUT, "flow_tilt.csv"))):
        r["phi"] = float(r["vth_abs"]) / float(r["vrms"])
        r["tilt_flow_deg"] = float(r["tilt_flow"]) * DEG; r["tilt_abs_deg"] = float(r["tilt_abs"]) * DEG
        for c in r:
            if c not in ("set", "run"): r[c] = str(r[c])
        G[(r["set"], round(float(r["rho"]), 2), abs(int(r["k"])))].append(r)
    return G


def curves(ax, G, key, sets=("pos", "neg"), legend=False):
    rhos = sorted({k[1] for k in G})
    for sset, cm in (("pos", plt.get_cmap("Reds")), ("neg", plt.get_cmap("Blues"))):
        if sset not in sets: continue
        for i, rho in enumerate(rhos):
            ks = sorted(k[2] for k in G if k[0] == sset and k[1] == rho)
            if len(ks) < 2: continue
            R = [float(G[(sset, rho, k)][0]["R"]) / LA for k in ks]
            q = [mean(G[(sset, rho, k)], key)[0] for k in ks]
            ax.plot(R, q, "-o", ms=2.5, lw=1.1, color=cm(0.3 + 0.7 * i / max(len(rhos) - 1, 1)),
                    label=f"{'+k' if sset == 'pos' else '−k'}, ρ = {rho:.2f}" if legend else None)
    ax.grid(color=GRID, lw=0.6); ax.spines[["top", "right"]].set_visible(False); ax.set_xlabel("activity number R / ℓa")


def maprow(fig, axs, G, key, norm, cmap, label):
    for ax, (sset, lab) in zip(axs, (("pos", "winding +k"), ("neg", "winding −|k|"))):
        tiles(ax, G, sset, key, norm, cmap)
        ax.set_title(f"{label}: {lab}", fontsize=9); ax.set_xlabel("ρ = 2k / P", fontsize=8); ax.set_ylabel("R / ℓa", fontsize=8)
        fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, pad=0.02, fraction=0.05)


def save(fig, pdf, name):
    fig.savefig(os.path.join(OUT, name + ".png"), dpi=150); fig.savefig(os.path.join(OUT, name + ".pdf")); pdf.savefig(fig); plt.close(fig)


def page_flow(pdf, G):
    fig = plt.figure(figsize=(17, 12)); gs = fig.add_gridspec(3, 3, height_ratios=[1.1, 1, 1])
    specs = [("vth_abs", "|⟨v_θ⟩|  (mean azimuthal flow)", Normalize(0, 1.0), "YlOrBr"),
             ("vrms", "v_rms = ⟨|v|²⟩^½", Normalize(0, 1.1), "YlOrBr"),
             ("phi", "circulation fraction |⟨v_θ⟩| / v_rms", Normalize(0, 1), "viridis")]
    for col, (key, name, norm, cmap) in enumerate(specs):
        ax = fig.add_subplot(gs[0, col]); curves(ax, G, key, legend=(col == 0)); ax.set_title(name, fontsize=10)
        ax.set_ylim(0, 1.15 if col < 2 else 1.02)
        if col == 0: ax.legend(frameon=False, fontsize=6, ncol=2, loc="upper left"); ax.set_ylabel("MPCD units (cell / time)")
        maprow(fig, [fig.add_subplot(gs[1, col]), fig.add_subplot(gs[2, col])], G, key, norm, plt.get_cmap(cmap), name.split("  ")[0])
    fig.suptitle("Mean flow in circles (ζ = 0.07, ρ ≤ 0.15, means over 3 seeds). ⟨·⟩ = spatial mean over cells with r < R − 1 at each "
                 "flow frame (every 100 time units); |·| taken per frame, then time-averaged (first 20% dropped), so circulation reversals do not cancel",
                 fontsize=9.5)
    fig.tight_layout(rect=(0, 0, 1, 0.965)); save(fig, pdf, "flow_p1_mean_flow")


def page_ens(pdf, G):
    fig = plt.figure(figsize=(17, 5.2)); gs = fig.add_gridspec(1, 3)
    ax = fig.add_subplot(gs[0, 0]); curves(ax, G, "ens", legend=True); ax.set_ylim(0, None)
    ax.set_title("enstrophy ⟨ω²/2⟩ vs activity number", fontsize=10); ax.set_ylabel("⟨ω²/2⟩ (time⁻²)")
    ax.legend(frameon=False, fontsize=6, ncol=2, loc="upper right")
    maprow(fig, [fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[0, 2])], G, "ens", Normalize(0.002, 0.014), plt.get_cmap("PuRd"), "⟨ω²/2⟩")
    fig.suptitle("Enstrophy in circles: ω = ∂x v_y − ∂y v_x (centred differences, cells away from the wall), spatial mean of ω²/2 per frame, "
                 "time-averaged; ζ = 0.07, ρ ≤ 0.15, 3 seeds", fontsize=9.5)
    fig.tight_layout(rect=(0, 0, 1, 0.93)); save(fig, pdf, "flow_p2_enstrophy")


def director_sketch(ax, charge, delta, title):
    """Wall arc (outward normal = up), a wall-sign defect below it, director field drawn as headless segments;
    the defect orientation psi = 90 deg + delta (radial arm of the defect)."""
    psi = np.pi / 2 + delta
    xs, ys = np.meshgrid(np.linspace(-1.6, 1.6, 15), np.linspace(-1.5, 1.0, 12))
    ph = np.arctan2(ys, xs)
    th = 0.5 * (ph - psi) + psi if charge > 0 else -0.5 * (ph - psi) + psi
    rr = np.hypot(xs, ys); m = (rr > 0.18) & (ys < 1.0)
    ax.quiver(xs[m], ys[m], np.cos(th[m]), np.sin(th[m]), headwidth=0, headlength=0, headaxislength=0, pivot="middle",
              scale=17, width=0.008, color=INK2)
    ax.add_patch(Arc((0, -9.0), 20.4, 20.4, theta1=78, theta2=102, color=INK, lw=2.4))      # wall, curvature centre below
    ax.annotate("", (0, 1.65), (0, 1.15), arrowprops=dict(arrowstyle="->", color=INK, lw=1, ls="--"))
    ax.text(0.08, 1.5, "normal", fontsize=7, color=INK)
    L = 1.05
    ax.annotate("", (L * np.cos(psi), L * np.sin(psi)), (0, 0), arrowprops=dict(arrowstyle="-|>", color="#c0392b" if charge > 0 else "#1f5fa8", lw=1.8))
    if charge > 0:
        ax.plot(0, 0, "o", color="#c0392b", ms=7)
    else:
        ax.plot(0, 0, "^", color="#1f5fa8", ms=7)
        for j in (1, 2):
            a = psi + j * 2 * np.pi / 3
            ax.plot([0, 0.6 * np.cos(a)], [0, 0.6 * np.sin(a)], color="#1f5fa8", lw=1, ls=":")
    if delta != 0:
        ax.add_patch(Arc((0, 0), 1.0, 1.0, theta1=90, theta2=90 + delta * DEG, color=INK, lw=1))
        ax.text(0.12, 0.62, f"Δ = {delta * DEG:+.0f}°", fontsize=7.5, ha="left", color=INK)
    ax.add_patch(FancyArrowPatch((1.3, -1.75), (-1.3, -1.75), arrowstyle="-|>", mutation_scale=10, color="#2e8b57", lw=1.4))
    ax.text(0, -2.05, "flow along the wall, v_θ > 0 (ccw)", fontsize=6.5, ha="center", color="#2e8b57")
    ax.set_xlim(-1.8, 1.8); ax.set_ylim(-2.2, 1.8); ax.set_aspect("equal"); ax.axis("off"); ax.set_title(title, fontsize=8)


def page_tilt(pdf, G):
    fig = plt.figure(figsize=(17, 10)); gs = fig.add_gridspec(2, 4, width_ratios=[1, 1, 0.55, 0.55])
    a = fig.add_subplot(gs[0, 0]); curves(a, G, "tilt_flow_deg", sets=("pos",), legend=True); a.set_ylim(0, 80)
    a.set_title("+k: tilt of the wall +½ (signed with the circulation)", fontsize=9); a.set_ylabel("Δ (degrees)")
    a.legend(frameon=False, fontsize=6, ncol=2, loc="lower right")
    b = fig.add_subplot(gs[0, 1]); curves(b, G, "tilt_flow_deg", sets=("neg",), legend=True); b.set_ylim(0, 30)
    b.set_title("−k: tilt of the wall −½ (three-fold, signed with the circulation)", fontsize=9); b.set_ylabel("Δ (degrees)")
    b.legend(frameon=False, fontsize=6, ncol=2, loc="lower right")
    for ax in (a, b): ax.axhline(0, color=INK2, lw=0.6)
    ax1 = fig.add_subplot(gs[1, 0]); tiles(ax1, G, "pos", "tilt_flow_deg", Normalize(0, 75), plt.get_cmap("Reds"))
    ax2 = fig.add_subplot(gs[1, 1]); tiles(ax2, G, "neg", "tilt_flow_deg", Normalize(0, 25), plt.get_cmap("Blues"))
    for ax, t, nm in ((ax1, "+k: Δ of the wall +½ (deg)", Normalize(0, 75)), (ax2, "−k: Δ of the wall −½ (deg)", Normalize(0, 25))):
        ax.set_title(t, fontsize=9); ax.set_xlabel("ρ = 2k / P", fontsize=8); ax.set_ylabel("R / ℓa", fontsize=8)
        fig.colorbar(plt.cm.ScalarMappable(norm=nm, cmap=plt.get_cmap("Reds" if ax is ax1 else "Blues")), ax=ax, pad=0.02, fraction=0.05)
    med_p = np.nanmedian([float(r["tilt_flow"]) for g in G.values() for r in g if r["set"] == "pos"])
    med_n = np.nanmedian([float(r["tilt_flow"]) for g in G.values() for r in g if r["set"] == "neg"])
    director_sketch(fig.add_subplot(gs[0, 2]), +1, 0.0, "+½, Δ = 0:\nradial arm points at the wall")
    director_sketch(fig.add_subplot(gs[0, 3]), +1, med_p, f"+½, Δ = {med_p * DEG:.0f}° (median):\narm leans downstream")
    director_sketch(fig.add_subplot(gs[1, 2]), -1, 0.0, "−½, Δ = 0:\none arm points at the wall")
    director_sketch(fig.add_subplot(gs[1, 3]), -1, med_n, f"−½, Δ = {med_n * DEG:.0f}° (median):\narms rotate downstream")
    fig.suptitle("Tilt Δ = ψ − φ of the wall-sign defects (within 6 cells of the wall): ψ = defect orientation (code's angle, the radial arm), "
                 "φ = polar angle of the defect (outward normal). Per frame: circular mean over the wall defects; signed with the circulation "
                 "sense sign⟨v_θ⟩ at the same time stamp; averaged over flow frames and 3 seeds. Δ > 0: rotated in the flow direction.", fontsize=9, wrap=True)
    fig.tight_layout(rect=(0, 0, 1, 0.93)); save(fig, pdf, "flow_p3_tilt")


def main():
    G = load()
    with PdfPages(os.path.join(OUT, "flow_maps.pdf")) as pdf:
        page_flow(pdf, G); page_ens(pdf, G); page_tilt(pdf, G)
    print("written", os.path.join(OUT, "flow_maps.pdf"))


if __name__ == "__main__":
    main()
