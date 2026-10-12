#!/usr/bin/env python3
"""
Wall Fourier figures for circles (main sweep, zeta = 0.07, rho <= 0.15), from wall_fourier_all.py.
  page 1  static: m*/m_lat, occupancy N_wall/m_lat, instantaneous order S(m*)/N_wall, static fraction f_static
  page 2  dynamic: lattice memory F(m*, 1000), decay time of the fluctuating part tau_fl, sliding speed / |<v_theta>|
Curves vs R/l_a (red +k, blue -k, one per rho) on top, (rho, R/l_a) maps for +k and -k below.
Writes analysis_fourier/fourier_maps.pdf and one PNG/PDF per page.
"""
import collections, csv, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import Normalize, LogNorm, TwoSlopeNorm

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from plots_flow import curves, maprow  # noqa: E402
OUT = os.path.join(HERE, "analysis_fourier")


def load():
    G = collections.defaultdict(list)
    for r in csv.DictReader(open(os.path.join(OUT, "wall_fourier_all.csv"))):
        N, ml = float(r["Nwall"]), float(r["mlat"])
        r["occ"] = str(N / ml if ml > 0 else np.nan); r["sn"] = str(float(r["Sstar"]) / max(N, 1e-9))
        r["ltaufl"] = str(np.log10(max(float(r["taufl"]), 1.0)))
        G[(r["set"], round(float(r["rho"]), 2), abs(int(r["k"])))].append(r)
    return G


def page(pdf, G, specs, title, name):
    n = len(specs)
    fig = plt.figure(figsize=(4.6 * n, 12)); gs = fig.add_gridspec(3, n, height_ratios=[1.1, 1, 1])
    for col, (key, lab, norm, cmap, ylim) in enumerate(specs):
        ax = fig.add_subplot(gs[0, col]); curves(ax, G, key, legend=(col == 0)); ax.set_title(lab, fontsize=10); ax.set_ylim(*ylim)
        if col == 0: ax.legend(frameon=False, fontsize=5.5, ncol=2, loc="lower right")
        maprow(fig, [fig.add_subplot(gs[1, col]), fig.add_subplot(gs[2, col])], G, key, norm, plt.get_cmap(cmap), lab.split(" (")[0])
    fig.suptitle(title, fontsize=9.5); fig.tight_layout(rect=(0, 0, 1, 0.955))
    fig.savefig(os.path.join(OUT, name + ".png"), dpi=150); fig.savefig(os.path.join(OUT, name + ".pdf")); pdf.savefig(fig); plt.close(fig)


def main():
    G = load()
    with PdfPages(os.path.join(OUT, "fourier_maps.pdf")) as pdf:
        page(pdf, G, [("mratio", "m* / m_lat (lattice period)", TwoSlopeNorm(1, 0.5, 2), "PuOr_r", (0.4, 2.2)),
                      ("occ", "occupancy N_wall / m_lat", Normalize(0.5, 1.1), "viridis", (0.4, 1.15)),
                      ("sn", "instantaneous order S(m*) / N_wall", Normalize(0, 1), "viridis", (0, 1.0)),
                      ("fstatic", "static fraction f_static (pinned)", Normalize(0, 1), "viridis", (0, 1.02))],
             "Wall Fourier, static (circles, ζ = 0.07, ρ ≤ 0.15, 3 seeds). Wall-sign defects within 4 cells of the wall, ρ_m = Σ exp(i m φ_j) per frame. "
             "m* = peak of S(m) = ⟨|ρ_m|²⟩/⟨N⟩; m_lat = wall count for a +1 bulk (2k − 2 for +k, 2|k| + 2 for −k).\n"
             "S(m*)/N = 1: perfect instantaneous lattice; f_static = |⟨ρ_m*⟩|²/⟨|ρ_m*|²⟩ = 1: the lattice sits at fixed sites", "fourier_p1_static")
        page(pdf, G, [("F1000", "lattice memory |F(m*, t = 1000)|", Normalize(0, 1), "viridis", (0, 1.05)),
                      ("ltaufl", "log10 τ_fl (fluctuating part, time units)", Normalize(0.5, 3.3), "magma", (0.3, 3.4)),
                      ("v_ratio", "sliding speed / |⟨v_θ⟩|", TwoSlopeNorm(0, -0.1, 0.3), "RdBu_r", (-0.1, 0.3))],
             "Wall Fourier, dynamic at m*. F(t) = |⟨ρ(n+t) ρ*(n)⟩|/⟨|ρ|²⟩ (co-moving: a rigidly sliding lattice keeps F = 1); τ_fl = 1/e time of "
             "ρ − ⟨ρ⟩ (rearrangements around the pinned lattice);\nsliding speed = R · arg Σ ρ(n+1)ρ*(n) / (m* dt) per 100-time-unit window, signed "
             "with the circulation (> 0: lattice slides with the flow), divided by the mean flow |⟨v_θ⟩|", "fourier_p2_dynamic")
    print("written", os.path.join(OUT, "fourier_maps.pdf"))


if __name__ == "__main__":
    main()
