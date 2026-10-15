#!/usr/bin/env python3
"""
Figures comparing the edge/bulk criteria for circles (main sweep, rho <= 0.15), one family per PDF page.

Input: analysis_weekend/criteria_all.csv (criteria_all.py). Means over seeds per (winding sign, rho, |k|).
  page 1  bulk charge maps (rho, R/l_a) for every criterion; rows = +k, -k
  page 2  r* (Q = +1 crossing), f* = 1 - r*^2, knee depth 1 - r_knee/R; schematic cavities with the shaded region
  page 3  radial defect density across the knee -> no-knee transitions (+k at large R, +k at small R, -k with |k|)
  page 4  bulk charge vs R/l_a, one panel per criterion; red tones = +k, blue tones = -k, one curve per rho
Writes analysis_weekend/criteria_maps.pdf and one PNG per page.
"""
import csv, collections, glob, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap, Normalize
from matplotlib.patches import Rectangle, Circle, Wedge

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "cassini"))
LA = 22.21
INK, INK2, GRID = "#15171b", "#555a62", "#e9ebee"
CRIT = [("q_knee", "knee (plateau of Q)"), ("q_fixed", "fixed layer: r < 11R/12"), ("q_site", "pinned-site rule"),
        ("q_dens", "density minimum"), ("q_slow", "mobility (slow)"), ("q_tilt", "slow + tilted, |tilt| < 1.8 (+k only)")]
DIV = LinearSegmentedColormap.from_list("q", ["#104281", "#6da7ec", "#f2f2f0", "#eb6834", "#8f3413"])


def load_rows():
    rows = list(csv.DictReader(open(os.path.join(HERE, "analysis_weekend", "criteria_all.csv"))))
    G = collections.defaultdict(list)
    for r in rows:
        G[(r["set"], round(float(r["rho"]), 2), abs(int(r["k"])))].append(r)
    return G


def mean(rr, key):
    v = [float(r[key]) for r in rr if r[key] not in ("", "nan") and np.isfinite(float(r[key]))]
    return (np.mean(v), len(v)) if v else (np.nan, 0)


def tiles(ax, G, sset, key, norm, cmap, frac_min=0.67):
    rhos = sorted({k[1] for k in G if k[0] == sset})
    edges_x = np.r_[rhos[0] - 0.005, 0.5 * (np.array(rhos[1:]) + np.array(rhos[:-1])), rhos[-1] + 0.005]
    for i, rho in enumerate(rhos):
        ks = sorted(k[2] for k in G if k[0] == sset and k[1] == rho)
        Rs = np.array([float(G[(sset, rho, k)][0]["R"]) / LA for k in ks])
        ey = np.r_[Rs[0] - 0.5 * (Rs[1] - Rs[0]) if len(Rs) > 1 else Rs[0] - 0.1, 0.5 * (Rs[1:] + Rs[:-1]),
                   Rs[-1] + 0.5 * (Rs[-1] - Rs[-2]) if len(Rs) > 1 else Rs[-1] + 0.1]
        for j, k in enumerate(ks):
            rr = G[(sset, rho, k)]
            m, n = mean(rr, key)
            x0, x1, y0, y1 = edges_x[i], edges_x[i + 1], max(ey[j], 0), ey[j + 1]
            if n == 0:
                ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, facecolor="#d9dbe0", edgecolor="none", hatch="///", lw=0))
            else:
                ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, facecolor=cmap(norm(m)), edgecolor="white", lw=0.3,
                                       hatch="..." if n < frac_min * len(rr) else None))
    ax.set_xlim(edges_x[0], edges_x[-1]); ax.set_ylim(0, 3.0)


def page_maps(pdf, G):
    fig, axs = plt.subplots(2, 6, figsize=(17, 6.6), sharex=True, sharey=True)
    norm = TwoSlopeNorm(vmin=-0.5, vcenter=1.0, vmax=2.5)
    for row, (sset, lab) in enumerate((("pos", "winding +k"), ("neg", "winding −|k|"))):
        for col, (key, name) in enumerate(CRIT):
            ax = axs[row, col]
            if sset == "neg" and key == "q_tilt":
                ax.text(0.5, 0.5, "tilt not defined\nfor −½ wall defects", ha="center", va="center", transform=ax.transAxes, color=INK2, fontsize=8)
                ax.set_xlim(0.025, 0.155); ax.set_ylim(0, 3); continue
            tiles(ax, G, sset, key, norm, DIV)
            if row == 0: ax.set_title(name, fontsize=9)
            if col == 0: ax.set_ylabel(f"{lab}\nR / ℓa", fontsize=9)
            if row == 1: ax.set_xlabel("ρ = 2k / P", fontsize=8)
    sm = plt.cm.ScalarMappable(norm=norm, cmap=DIV)
    cb = fig.colorbar(sm, ax=axs, shrink=0.8, pad=0.01); cb.set_label("bulk charge q_b (mean over seeds)")
    fig.suptitle("Bulk charge in circles by criterion (ζ = 0.07, ρ ≤ 0.15, 3 seeds). Grey hatched: criterion undefined in all seeds "
                 "(no knee / no plateau); dotted: defined in fewer than 2 of 3 seeds", fontsize=10)
    fig.savefig(os.path.join(HERE, "analysis_weekend", "crit_p1_bulk_maps.png"), dpi=150); fig.savefig(os.path.join(HERE, "analysis_weekend", "crit_p1_bulk_maps.pdf")); pdf.savefig(fig); plt.close(fig)


def schematic(ax, rstar, title, color="#eb6834"):
    ax.add_patch(Circle((0, 0), 1, facecolor="#f2f2f0", edgecolor=INK, lw=1.2))
    ax.add_patch(Wedge((0, 0), 1, 0, 360, width=1 - rstar, facecolor=color, alpha=0.55, edgecolor="none"))
    ax.add_patch(Circle((0, 0), rstar, fill=False, edgecolor=INK2, lw=0.8, ls="--"))
    ax.set_xlim(-1.1, 1.1); ax.set_ylim(-1.1, 1.1); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(title, fontsize=8)


def page_rstar(pdf, G):
    fig = plt.figure(figsize=(17, 7.4))
    gs = fig.add_gridspec(2, 5, width_ratios=[1, 1, 1, 0.55, 0.55])
    specs = [("rstar", "r* / R  (where Q crosses +1)", Normalize(0.6, 1.0), plt.get_cmap("viridis")),
             ("fstar", "f* = 1 − (r*/R)²  (area fraction outside r*)", Normalize(0.0, 0.6), plt.get_cmap("magma_r")),
             ("kneed", "knee depth  1 − r_knee / R", Normalize(0.0, 0.35), plt.get_cmap("cividis_r"))]
    for row, (sset, lab) in enumerate((("pos", "winding +k"), ("neg", "winding −|k|"))):
        for col, (key, name, norm, cmap) in enumerate(specs):
            ax = fig.add_subplot(gs[row, col])
            if key == "kneed":
                for g in G.values():
                    for r in g:
                        r["kneed"] = 1 - float(r["knee_r"]) if r["knee_r"] not in ("", "nan") and np.isfinite(float(r["knee_r"])) else "nan"
            tiles(ax, G, sset, key, norm, cmap)
            sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap); fig.colorbar(sm, ax=ax, pad=0.02, fraction=0.05)
            if row == 0: ax.set_title(name, fontsize=9)
            if col == 0: ax.set_ylabel(f"{lab}\nR / ℓa", fontsize=9)
            if row == 1: ax.set_xlabel("ρ = 2k / P", fontsize=8)
    # schematics
    a1 = fig.add_subplot(gs[0, 3]); schematic(a1, 0.87, "circle with a layer:\nr* ≈ 0.87, f* ≈ 0.25\n(shaded: outside r*)")
    a2 = fig.add_subplot(gs[0, 4]); schematic(a2, 0.70, "no layer:\nr* ≈ 0.70, f* ≈ 0.5", color="#8f3413")
    a3 = fig.add_subplot(gs[1, 3]); schematic(a3, 0.88, "knee at r/R ≈ 0.88:\nknee depth ≈ 0.12", color="#2a78d6")
    a4 = fig.add_subplot(gs[1, 4]); a4.axis("off")
    a4.text(0, 0.9, "r*: radius where the enclosed\ncharge first reaches +1 (scanning\ninwards from the wall).\n\nf*: area fraction between r* and\nthe wall (shaded).\n\nknee depth: distance of the\nplateau of Q(r) from the wall,\nin units of R.\n\nMirror runs: r* only where the\nbulk reaches +1 (|k| ≳ 13).",
            fontsize=8, va="top", color=INK)
    fig.suptitle("Where the bulk starts: r*, f* and knee depth (circles, ζ = 0.07, ρ ≤ 0.15). Grey hatched: undefined", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(os.path.join(HERE, "analysis_weekend", "crit_p2_rstar_fstar_knee.png"), dpi=150); fig.savefig(os.path.join(HERE, "analysis_weekend", "crit_p2_rstar_fstar_knee.pdf")); pdf.savefig(fig); plt.close(fig)


def page_density(pdf, G):
    from lobe_profiles import profile
    from charge_knee import knee
    seqs = [("+k, ρ = 0.05: growing R (layer lost at large R)", "main_act0.07", "rho0.05", [4, 6, 8, 9, 10], "k"),
            ("+k, ρ = 0.12: growing R (no layer → layer at small R)", "main_act0.07", "rho0.12", [3, 4, 5, 6, 8, 10], "k"),
            ("−k, ρ = 0.08: growing |k| (bulk builds up to +1)", "main_act0.07_neg", "rho0.08", [4, 6, 8, 10, 13, 16], "km")]
    fig, axs = plt.subplots(2, 3, figsize=(16, 8.0))
    for col, (title, st, rho, ks, pre) in enumerate(seqs):
        cm = plt.get_cmap("Reds" if st == "main_act0.07" else "Blues")
        for i, k in enumerate(ks):
            runs = sorted(glob.glob(os.path.join(HERE, "server_main", st, rho, f"{pre}{k:02d}", "s[01][0-9]")))
            if not runs: continue
            P = [profile(r, 120, 0.2) for r in runs]
            x = P[0]["x"]; npl = np.mean([p["npl"] for p in P], 0); nmi = np.mean([p["nmi"] for p in P], 0)
            Q = np.mean([np.r_[p["Q"]] for p in P], 0)
            R = float(open(os.path.join(runs[0], "input.json")).read().split('"R":')[1].split(",")[0])
            col_ = cm(0.35 + 0.65 * i / max(len(ks) - 1, 1))
            kn = [knee(p) for p in P] if st == "main_act0.07" else [None]
            has = sum(1 for q in kn if q) if st == "main_act0.07" else None
            ls = "-" if (has is None or has >= 2) else "--"
            wall = npl if st == "main_act0.07" else nmi
            axs[0, col].plot(x, np.convolve(wall, np.ones(3) / 3, "same"), color=col_, lw=1.6, ls=ls,
                             label=f"|k| = {k}, R = {R / LA:.2f} ℓa" + ("" if has is None else f" (knee {has}/3)"))
            axs[1, col].plot(P[0]["edges"][1:], Q, color=col_, lw=1.6, ls=ls)
        axs[0, col].set_yscale("log"); axs[0, col].set_ylim(1e-4, 0.3); axs[0, col].set_title(title, fontsize=9)
        axs[0, col].set_ylabel(("+½" if st == "main_act0.07" else "−½") + " density (per cell² per frame)", fontsize=8)
        axs[1, col].axhline(1, color=INK2, lw=0.7, ls=":"); axs[1, col].set_xlabel("r / R"); axs[1, col].set_ylabel("Q(r)", fontsize=8)
        axs[0, col].legend(frameon=False, fontsize=7, loc="upper left")
        for a in axs[:, col]:
            a.grid(color=GRID, lw=0.6); a.spines[["top", "right"]].set_visible(False); a.set_xlim(0, 1)
    fig.suptitle("Radial density of the wall-sign defects (top) and enclosed charge (bottom) across the knee → no-knee transitions "
                 "(means over 3 seeds; dashed: knee in fewer than 2 seeds)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(os.path.join(HERE, "analysis_weekend", "crit_p3_density_transitions.png"), dpi=150); fig.savefig(os.path.join(HERE, "analysis_weekend", "crit_p3_density_transitions.pdf")); pdf.savefig(fig); plt.close(fig)


def page_vsR(pdf, G):
    fig, axs = plt.subplots(2, 3, figsize=(16, 8.4), sharex=True)
    axs = axs.ravel()
    rhos = sorted({k[1] for k in G})
    reds = plt.get_cmap("Reds"); blues = plt.get_cmap("Blues")
    for ax, (key, name) in zip(axs, CRIT):
        for sset, cm in (("pos", reds), ("neg", blues)):
            for i, rho in enumerate(rhos):
                ks = sorted(k[2] for k in G if k[0] == sset and k[1] == rho)
                if len(ks) < 2: continue
                R = [float(G[(sset, rho, k)][0]["R"]) / LA for k in ks]
                q = [mean(G[(sset, rho, k)], key)[0] for k in ks]
                if np.all(np.isnan(q)): continue
                ax.plot(R, q, "-o", ms=3, lw=1.2, color=cm(0.3 + 0.7 * i / max(len(rhos) - 1, 1)),
                        label=(f"{'+k' if sset == 'pos' else '−k'}, ρ = {rho:.2f}") if key == "q_knee" else None)
        ax.axhline(1, color=INK2, lw=0.8, ls=":"); ax.set_title(name, fontsize=9); ax.set_ylim(-1.0, 3.0)
        ax.grid(color=GRID, lw=0.6); ax.spines[["top", "right"]].set_visible(False)
    for ax in axs[3:]:
        ax.set_xlabel("activity number R / ℓa")
    axs[0].set_ylabel("bulk charge q_b"); axs[3].set_ylabel("bulk charge q_b")
    axs[0].legend(frameon=False, fontsize=6, ncol=2, loc="lower right")
    fig.suptitle("Bulk charge vs activity number for each linear density (red tones: +k, light → dark = low → high ρ; blue tones: −|k|). "
                 "Circles, ζ = 0.07, means over 3 seeds", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(os.path.join(HERE, "analysis_weekend", "crit_p4_bulk_vs_R.png"), dpi=150); fig.savefig(os.path.join(HERE, "analysis_weekend", "crit_p4_bulk_vs_R.pdf")); pdf.savefig(fig); plt.close(fig)


def main():
    G = load_rows()
    with PdfPages(os.path.join(HERE, "analysis_weekend", "criteria_maps.pdf")) as pdf:
        page_maps(pdf, G); page_rstar(pdf, G); page_density(pdf, G); page_vsR(pdf, G)
    print("written analysis_weekend/criteria_maps.pdf")


if __name__ == "__main__":
    main()
