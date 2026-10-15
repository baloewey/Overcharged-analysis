"""Shared style for the PRL composite figures (2026-10-15). Figures are drawn at their printed size (PRL: column 3.375 in, page 7.0 in),
so fonts are given in points and NOT rescaled: run the figure scripts with FIGSTYLE_SCALE=1 (figs/run.sh does this) to keep the house
capitalisation without the x1.3 font scale. Every quantity on an axis is named and has units; panel labels (a), (b), ... in bold at the top left.
Palette: + wall/+1/2 red, -1/2 blue, passive/control grey; k = 4/6/8 blue/red/green as in the report; q (winding density) viridis.
Output: vector PDF + PNG (300 dpi) next to the script; large collections rasterised by the figstyle hook."""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

COL, PAGE = 3.375, 7.0                     # inches
RED, BLUE, GREEN, GREY, DARK = "#c2452d", "#1f5fa8", "#2e7d3c", "0.55", "0.25"
KCOL = {4: BLUE, 6: RED, 8: GREEN}
H = os.path.dirname(os.path.abspath(__file__))
ROOT = "/home/baloewe/Physics/Overcharged"

plt.rcParams.update({
    "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 6.5,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "lines.linewidth": 1.1, "lines.markersize": 3.5, "legend.frameon": False, "legend.handlelength": 1.6,
    "axes.spines.top": False, "axes.spines.right": False, "savefig.dpi": 300, "pdf.fonttype": 42, "font.family": "DejaVu Sans",
    "mathtext.fontset": "dejavusans",
})


def label(ax, s, x=-0.02, y=1.02, **kw):
    """panel label (a), (b), ... at the top left of an axes (axes coordinates)"""
    ax.text(x, y, f"({s})", transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom", ha="right", **kw)


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(H, f"{name}.{ext}"), bbox_inches="tight", pad_inches=0.02)
    print("wrote", os.path.join(H, name + ".pdf"))
