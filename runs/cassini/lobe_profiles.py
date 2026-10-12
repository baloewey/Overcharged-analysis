#!/usr/bin/env python3
"""
Defect profiles in Cassini lobe coordinates (circle: a = 0, then c = r).

Lobe coordinate of a point: c = sqrt(d1 d2) (d1, d2 distances to the foci at +-a), so c = 0 at the foci,
c = a on the lemniscate (below it the level set is two separate lobe ovals), c = b on the wall.
Per run (first --tfrac of defect frames dropped):
  n_+(c), n_-(c)   number density in shells [c, c + dc] = defects per frame / shell area (area from a fine grid)
  Q(c)             accumulated charge inside level c (net charge of defects with sqrt(d1 d2) < c)
  Q_L(c), Q_R(c)   the same per lobe (x < 0, x > 0), meaningful for c < a
Shells: --nshell bins uniform in c/b. Writes PREFIX.png and PREFIX.csv; prints Q at a few c/b values.

Usage: lobe_profiles.py RUN [RUN ...] [--labels ...] [--out PREFIX]
"""
import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "cavity"))
from defect_spatial import load  # noqa: E402

PLUS, MINUS, INK, INK2, GRID = "#2a78d6", "#eb6834", "#15171b", "#555a62", "#e9ebee"


def geometry(run):
    bc = json.load(open(os.path.join(run, "input.json")))["BC"][0]
    a = bc.get("cassiniA", 0.0) if bc.get("shape", 0) == 1 else 0.0
    return a, bc["R"]


def cval(a, x, y):
    return np.sqrt(np.hypot(x - a, y) * np.hypot(x + a, y))


def shell_areas(a, b, edges, h=0.1):
    X = np.arange(-np.hypot(a, b) - 1, np.hypot(a, b) + 1, h)
    xx, yy = np.meshgrid(X, X, indexing="ij")
    c = cval(a, xx, yy).ravel()
    return np.histogram(c, edges)[0] * h * h


def ellipse_eps(run):
    """Axis ratio eps = A/B of an elliptical wall (shape 0 with aInv[0] != aInv[1]); None otherwise."""
    bc = json.load(open(os.path.join(run, "input.json")))["BC"][0]
    if bc.get("shape", 0) == 0 and abs(bc["aInv"][0] - bc["aInv"][1]) > 1e-9:
        return bc["aInv"][0] / bc["aInv"][1], bc["R"] * bc["aInv"][1]
    return None


def profile(run, nshell, tfrac):
    a, b = geometry(run)
    d = load(run, tfrac)
    fr = d["frames"]
    allf = np.vstack([f for f in fr if len(f)])
    ell = ellipse_eps(run)
    if ell is not None:
        # ellipse: levels are scaled copies of the wall, c = sqrt((x/eps)^2 + y^2) (c = b on the wall, b = semi-minor axis);
        # area inside level c = pi eps c^2
        eps, b = ell; a = 0.0
        edges = np.linspace(0, b, nshell + 1)
        A = np.pi * eps * np.diff(edges ** 2)
        c = np.hypot(allf[:, 0] / eps, allf[:, 1])
    else:
        edges = np.linspace(0, b, nshell + 1)
        A = shell_areas(a, b, edges)
        c = cval(a, allf[:, 0], allf[:, 1])
    nf = len(fr)
    npl = np.histogram(c[allf[:, 2] > 0], edges)[0] / nf / np.maximum(A, 1e-9)
    nmi = np.histogram(c[allf[:, 2] < 0], edges)[0] / nf / np.maximum(A, 1e-9)
    q = np.histogram(c, edges, weights=allf[:, 2])[0] / nf
    qL = np.histogram(c[allf[:, 0] < 0], edges, weights=allf[allf[:, 0] < 0, 2])[0] / nf
    qR = np.histogram(c[allf[:, 0] > 0], edges, weights=allf[allf[:, 0] > 0, 2])[0] / nf
    return dict(a=a, b=b, k=d["k"], x=0.5 * (edges[1:] + edges[:-1]) / b, edges=edges / b, npl=npl, nmi=nmi,
                Q=np.cumsum(q), QL=np.cumsum(qL), QR=np.cumsum(qR))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--labels", nargs="+")
    ap.add_argument("--nshell", type=int, default=60)
    ap.add_argument("--tfrac", type=float, default=0.2)
    ap.add_argument("--out", default="lobe_profiles")
    a = ap.parse_args()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    labels = a.labels or [os.path.basename(os.path.normpath(r)) for r in a.runs]
    P = [profile(r, a.nshell, a.tfrac) for r in a.runs]
    cuts = [0.5, 0.7, 0.8, 0.85, 0.9, 0.95]
    print("accumulated charge Q(c) inside Cassini level c (c/b = 1 is the wall; Q(1) = k)")
    print(f"{'run':<28} {'k':>3} {'a/b':>5}  " + "  ".join(f"c/b<{x:<4}" for x in cuts) + "   lobes at c = a: Q_L  Q_R")
    with open(a.out + ".csv", "w") as fo:
        fo.write("run,label,k,a,b,c_over_b,n_plus,n_minus,Q,Q_left,Q_right\n")
        for r, p, lab in zip(a.runs, P, labels):
            for i in range(len(p["x"])):
                fo.write(f"{r},{lab},{p['k']},{p['a']},{p['b']},{p['x'][i]:.4f},{p['npl'][i]:.6g},{p['nmi'][i]:.6g},"
                         f"{p['Q'][i]:.4f},{p['QL'][i]:.4f},{p['QR'][i]:.4f}\n")
            Qc = [np.interp(x, p["edges"][1:], p["Q"]) for x in cuts]
            if p["a"] > 0:
                ia = np.searchsorted(p["edges"][1:], p["a"] / p["b"])
                lob = f"   {p['QL'][min(ia, len(p['QL']) - 1)]:5.2f} {p['QR'][min(ia, len(p['QR']) - 1)]:5.2f}"
            else:
                lob = ""
            print(f"{lab:<28} {p['k']:>3} {p['a'] / p['b']:5.2f}  " + "  ".join(f"{q:8.2f}" for q in Qc) + lob)
    n = len(P)
    fig, axs = plt.subplots(2, n, figsize=(3.0 * n, 5.6), squeeze=False, sharex=True)
    for col, (p, lab) in enumerate(zip(P, labels)):
        ax = axs[0, col]
        ax.plot(p["x"], p["npl"], color=PLUS, lw=1.6, label="+½")
        ax.plot(p["x"], p["nmi"], color=MINUS, lw=1.6, label="−½")
        ax.set_yscale("log"); ax.set_ylim(1e-5, 0.3); ax.set_title(lab, fontsize=9, color=INK)
        ax2 = axs[1, col]
        ax2.plot(p["x"], p["Q"], color=INK, lw=1.8, label="whole cavity")
        if p["a"] > 0:
            m = p["x"] <= p["a"] / p["b"]
            ax2.plot(p["x"][m], p["QL"][m], color="#1baf7a", lw=1.4, label="left lobe")
            ax2.plot(p["x"][m], p["QR"][m], color="#7c5cc4", lw=1.4, ls="--", label="right lobe")
        ax2.axhline(1, color=INK2, lw=0.7, ls=":"); ax2.set_xlabel("c / b  (0 foci, 1 wall)")
        for x in (ax, ax2):
            x.grid(color=GRID, lw=0.6); x.spines[["top", "right"]].set_visible(False)
            if p["a"] > 0:
                x.axvline(p["a"] / p["b"], color=MINUS, lw=0.8, ls="--")
        ax2.set_ylim(-1.5, p["k"] + 0.5)
    axs[0, 0].set_ylabel("defects per cell² per frame"); axs[1, 0].set_ylabel("charge inside level c")
    axs[0, 0].legend(frameon=False, fontsize=8); axs[1, -1].legend(frameon=False, fontsize=7)
    fig.tight_layout(); fig.savefig(a.out + ".png", dpi=150)
    print("written", a.out + ".png", a.out + ".csv")


if __name__ == "__main__":
    main()
