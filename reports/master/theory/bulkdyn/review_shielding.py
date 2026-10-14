#!/usr/bin/env python3
"""Review of the shielding maps (TO REVIEW, RESEARCH.md 2026-10-14), circle checks (analytic one-constant disk, same-sign wall images).
(1) Missing wells/peaks near some layer defects: plot up to 0.99 of the radius, mark each defect with its depth, and separate the probe's OWN wall term
    (-q^2 ln(R^2 - r^2), repulsion from the wall, diverges at the wall) from the LAYER term (probe x layer charges + their images).
(2) Asymmetry of the interior field: for measured full (10) layers, range of the layer term inside 0.6 R for (a) the measured layer, (b) the same defects re-spaced
    evenly in angle at their measured depths, (c) measured angles at the mean depth, (d) ideal (even angles, mean depth).
-> fig_review_shielding.png (+pdf), review_shielding.txt"""
import os, sys, json, numpy as np
from scipy.spatial import cKDTree
R_ = "/home/baloewe/Physics/Overcharged/runs"; H = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, R_ + "/cavity"); sys.path.insert(0, os.path.join(H, "..", "energetics"))
from edge_criteria import read_with_angle  # noqa
import energetics as En  # noqa
PIK = np.pi * 32.5; RUNS = [f"{R_}/fineflow/eps1.0_s1", f"{R_}/fineflow2/eps1.0_s3", f"{R_}/fineflow2/eps1.0_s4"]


def layer_term(zs, qs, R, X, qt=0.5):
    v = np.zeros(X.shape)
    for z, q in zip(zs, qs): v += -2 * qt * q * (np.log(np.abs(X - z)) + np.log(np.abs(R - X * np.conj(z) / R)))
    return v * PIK


def self_term(R, X, qt=0.5):
    return -qt * qt * np.log(np.maximum(R ** 2 - np.abs(X) ** 2, 1e-9)) * PIK


if __name__ == "__main__":
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm
    F = []
    for run in RUNS:
        bc = json.load(open(run + "/input.json"))["BC"][0]; R = bc["R"]; W, s, P = En.wall_polyline(bc); tree = cKDTree(W)
        t, fr = read_with_angle(run, np.array(bc["Q"][:2], float), 0.5)
        for f in fr[::3]:
            if not len(f): continue
            d = tree.query(f[:, :2])[0]; m = (d <= 4) & (f[:, 2] > 0)
            if m.sum() in (9, 10): F.append((int(m.sum()), f[m, 0] + 1j * f[m, 1], R))
    R = F[0][2]; xs = np.arange(-R, R + 0.01, 0.5); XX, YY = np.meshgrid(xs, xs, indexing="ij"); X = XX + 1j * YY
    ins99 = np.abs(X) < 0.99 * R; inn = np.abs(X) < 0.6 * R
    # (2) decomposition over full-layer frames
    rows = []
    for n, z, _ in F:
        if n != 10: continue
        r, ph = np.abs(z), np.angle(z); o = np.argsort(ph); r, ph = r[o], ph[o]
        ph_even = ph.mean() + np.unwrap(ph - ph.mean()).mean() * 0 + (np.arange(10) - 4.5) * 2 * np.pi / 10
        ph_even = ph[0] + np.arange(10) * 2 * np.pi / 10; ph_even += np.mean(np.angle(np.exp(1j * (ph - ph_even))))
        q = np.full(10, 0.5); rng = lambda zz: np.ptp(layer_term(zz, q, R, X[inn]))
        rows.append((rng(z), rng(r * np.exp(1j * ph_even)), rng(r.mean() * np.exp(1j * ph)), rng(r.mean() * np.exp(1j * ph_even)), r.std(), np.std(np.angle(np.exp(1j * (ph - ph_even))))))
    A = np.array(rows)
    L = ["REVIEW of the shielding maps, circles (k 6, area 2500, measured full layers of 10, %d frames). Range of the LAYER term inside 0.6 R (kBT), median [IQR]:" % len(A)]
    for i, lab in enumerate(("measured layer", "measured depths, even angles", "measured angles, mean depth", "ideal (even angles, mean depth)")):
        L.append(f"  {lab:34s}: {np.median(A[:,i]):6.1f} [{np.percentile(A[:,i],25):.1f}, {np.percentile(A[:,i],75):.1f}]")
    L.append(f"  layer irregularity: depth s.d. {np.median(A[:,4]):.2f} cells, angular s.d. {np.degrees(np.median(A[:,5])):.1f} deg (site spacing 36 deg)")
    # (1) maps: one full and one vacancy frame, total vs layer term, up to 0.99 R
    fig, axs = plt.subplots(2, 3, figsize=(15, 9.4))
    for row, nn in enumerate((10, 9)):
        n, z, _ = next(x for x in F if x[0] == nn)
        lt = layer_term(z, np.full(len(z), 0.5), R, X); st = self_term(R, X); tot = lt + st
        for col, (V, ttl) in enumerate(((tot, "total (layer + probe's own wall term)"), (lt, "layer term only"), (st, "probe's own wall term only"))):
            ax = axs[row, col]; Vm = np.where(ins99, V - np.nanmin(np.where(inn, V, np.nan)), np.nan); hi = np.nanpercentile(np.abs(Vm[ins99]), 95)
            im = ax.pcolormesh(xs, xs, Vm.T, cmap="RdBu_r", shading="nearest", norm=TwoSlopeNorm(0, -hi, hi))
            th = np.linspace(0, 2 * np.pi, 300); ax.plot(R * np.cos(th), R * np.sin(th), "k-", lw=1); ax.plot(0.9 * R * np.cos(th), 0.9 * R * np.sin(th), "k:", lw=0.6)
            ax.plot(z.real, z.imag, "o", ms=5, color="#1f5fa8", mec="w", mew=0.6)
            if col == 0:
                for zz in z: ax.annotate(f"{R - abs(zz):.1f}", (zz.real, zz.imag), xytext=(4, 4), textcoords="offset points", fontsize=7)
            ax.set_aspect("equal"); ax.axis("off"); ax.set_title(f"({'abcdef'[3*row+col]}) {nn} wall $+1/2$: {ttl}", fontsize=9.5)
            fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02, shrink=0.8).set_label("$V - V_\\mathrm{min, inner\\ 60\\%}$ ($k_BT$)", fontsize=8)
    fig.text(0.5, 0.01, "Dotted circle: 0.9 R, the edge of the earlier maps. Numbers in (a, d): depth of each wall defect (cells).", ha="center", fontsize=9)
    fig.tight_layout(rect=(0, 0.03, 1, 1)); fig.savefig(os.path.join(H, "fig_review_shielding.png"), dpi=180)
    print("\n".join(L)); open(os.path.join(H, "review_shielding.txt"), "w").write("\n".join(L) + "\n")
