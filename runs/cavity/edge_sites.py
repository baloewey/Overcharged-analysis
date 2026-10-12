#!/usr/bin/env python3
"""
Edge/bulk classification from pinned sites, and seed-pooled ring-resolved pair correlations.

Sites: local maxima of the time-averaged +1/2 density map (0.5-cell bins, Gaussian smoothing
--sigma cells) above --thresh (defects per cell^2 per frame). No radial cut and no circle
geometry is used, so the same rule applies to any cavity shape. A +1/2 defect is "edge" if it
lies within --dsite cells of a site; everything else (and every -1/2) is bulk.

Per run: number of sites, mean edge occupancy N_E, bulk charge q_b = net charge of bulk defects
(<q_b>, P(q_b = +1)), histogram of the distance to the nearest site, and agreement with the
density cut (s > --dcut is bulk): fraction of +1/2 defects classified the same way.

Pooled pair correlations: g_ab(r | s) (see defect_spatial.ring_gr) summed over all runs of the
same (k, R) (seeds), and pooled into zones (edge ring s < 4, outer bulk, centre) over all k.

Usage:
  edge_sites.py RUN [RUN ...] [--dsite 0 (auto)] [--out PREFIX]
"""
import argparse
import collections
import os
import sys

import numpy as np
from scipy.ndimage import gaussian_filter, maximum_filter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from defect_spatial import load, ring_gr  # noqa: E402

PLUS, MINUS, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, GRID = "#15171b", "#555a62", "#e9ebee"


def find_sites(frames, R, sigma, thresh, binw=0.5):
    allp = np.vstack([f[f[:, 2] > 0, :2] for f in frames if len(f)])
    edges = np.arange(-R - binw, R + binw * 1.01, binw)
    H, _, _ = np.histogram2d(allp[:, 0], allp[:, 1], [edges, edges])
    H = gaussian_filter(H / (len(frames) * binw ** 2), sigma / binw)
    peaks = (H == maximum_filter(H, size=int(round(2 * sigma / binw)) * 2 + 1)) & (H > thresh)
    ix, iy = np.nonzero(peaks)
    c = 0.5 * (edges[1:] + edges[:-1])
    return np.column_stack([c[ix], c[iy]]), H[ix, iy]


def classify(d, sigma, thresh, dsite, dcut, dfrac=0.4):
    R, fr = d["R"], d["frames"]
    sites, hval = find_sites(fr, R, sigma, thresh)
    if dsite <= 0:   # data-driven: minimum of the smoothed distance histogram between 1.5 and 8 cells
        dall = []
        for f in fr:
            p = f[:, 2] > 0
            if len(sites) and p.any():
                dall.append(np.hypot(*(f[p, None, :2] - sites[None, :, :]).transpose(2, 0, 1)).min(1))
        dall = np.concatenate(dall) if dall else np.zeros(0)
        e = np.arange(0, 12.01, 0.25)
        h = np.histogram(dall, e)[0].astype(float)
        h = np.convolve(h, np.ones(5) / 5, "same")
        c = 0.5 * (e[1:] + e[:-1])
        m = (c >= 1.5) & (c <= 8)
        dsite = c[m][np.argmin(h[m])] if m.any() and h[m].any() else 3.0
    # cap at a fraction of the median nearest-neighbour site spacing: where the bulk +1/2 ring overlaps the
    # layer there is no clear dip and the minimum runs to the window edge, capturing bulk defects
    capped = False
    if dfrac > 0 and len(sites) > 1:
        dd = np.hypot(*(sites[:, None, :] - sites[None, :, :]).transpose(2, 0, 1))
        np.fill_diagonal(dd, np.inf)
        cap = dfrac * np.median(dd.min(1))
        if dsite > cap:
            dsite, capped = cap, True
    qb, ne, dist, agree, nplus = [], [], [], 0, 0
    for f in fr:
        p = f[:, 2] > 0
        if len(sites) and p.any():
            dd = np.hypot(*(f[p, None, :2] - sites[None, :, :]).transpose(2, 0, 1)).min(1)
        else:
            dd = np.full(p.sum(), np.inf)
        dist.append(dd)
        edge = np.zeros(len(f), bool)
        edge[np.where(p)[0]] = dd < dsite
        qb.append(f[~edge, 2].sum())
        ne.append(edge.sum())
        s = R - np.hypot(f[p, 0], f[p, 1])
        agree += int(((dd < dsite) == (s <= dcut)).sum())
        nplus += int(p.sum())
    qb, ne = np.array(qb), np.array(ne)
    nb = 10
    blocks = np.array_split(qb, nb)
    err = np.std([b.mean() for b in blocks], ddof=1) / np.sqrt(nb)
    return dict(sites=sites, nsite=len(sites), site_s=R - np.hypot(sites[:, 0], sites[:, 1]) if len(sites) else np.zeros(0),
                qb=qb.mean(), qb_err=err, P1=np.mean(np.isclose(qb, 1)), NE=ne.mean(),
                dist=np.concatenate(dist), agree=agree / max(nplus, 1), dsite=dsite, capped=capped)


def style(ax):
    ax.tick_params(colors=INK2, labelsize=8)
    for sp in ax.spines.values():
        sp.set_color(INK2)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(color=GRID, lw=0.6)
    ax.set_axisbelow(True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--tfrac", type=float, default=0.2)
    ap.add_argument("--sigma", type=float, default=1.0)
    ap.add_argument("--thresh", type=float, default=0.02)
    ap.add_argument("--dsite", type=float, default=0.0, help="capture radius; <= 0: from the minimum of the distance histogram")
    ap.add_argument("--dcut", type=float, default=4.2, help="density cut for the comparison")
    ap.add_argument("--dfrac", type=float, default=0.4, help="cap the capture radius at dfrac x median site spacing (0 = no cap)")
    ap.add_argument("--no-gr", action="store_true")
    ap.add_argument("--out", default="edge_sites")
    a = ap.parse_args()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows, groups, dists = [], collections.defaultdict(list), {}
    for run in a.runs:
        d = load(run, a.tfrac)
        c = classify(d, a.sigma, a.thresh, a.dsite, a.dcut, a.dfrac)
        k = d["k"]
        rows.append((run, k, d["R"], d["act"], c))
        groups[(k, round(d["R"], 2), d["act"])].append(d)
        dists.setdefault(k, []).append(c["dist"])
        print(f"{run:<34} k={k:>3}  sites={c['nsite']:>3} (2k-2={2*k-2:>3})  site s={np.mean(c['site_s']) if c['nsite'] else np.nan:4.1f}"
              f"  N_E={c['NE']:5.2f}  q_b={c['qb']:5.2f}±{c['qb_err']:.2f}  P1={c['P1']:.2f}  agree(d={a.dcut:g})={c['agree']:.3f}  dsite={c['dsite']:.2f}{' (capped)' if c['capped'] else ''}", flush=True)

    with open(a.out + ".csv", "w") as fo:
        fo.write("run,k,R,act,nsite,site_s_mean,dsite,N_edge,qb,qb_err,P1,agree_density_cut\n")
        for run, k, R, act, c in rows:
            fo.write(f"{run},{k},{R},{act},{c['nsite']},{np.mean(c['site_s']) if c['nsite'] else ''},{c['dsite']:.2f},"
                     f"{c['NE']:.4f},{c['qb']:.4f},{c['qb_err']:.4f},{c['P1']:.4f},{c['agree']:.4f}\n")

    # distance-to-nearest-site histograms (per k, all seeds)
    ks = sorted(dists)
    fig, ax = plt.subplots(figsize=(6.4, 3.8)); style(ax)
    cmap = plt.get_cmap("Blues")
    for i, k in enumerate(ks):
        dd = np.concatenate(dists[k]); dd = dd[np.isfinite(dd)]
        h, e = np.histogram(dd, np.arange(0, 20.25, 0.25), density=True)
        ax.plot(0.5 * (e[1:] + e[:-1]), h, color=cmap(0.35 + 0.65 * i / max(len(ks) - 1, 1)), lw=1.5, label=f"k = {k}")
    ax.axvline(a.dsite, color=INK2, lw=0.8, ls=":")
    ax.set_yscale("log"); ax.set_xlabel("distance of a +½ defect to the nearest site (cells)", fontsize=9, color=INK2)
    ax.set_ylabel("probability density", fontsize=9, color=INK2)
    ax.legend(frameon=False, fontsize=7, ncol=2)
    fig.tight_layout(); fig.savefig(a.out + "_sitedist.png", dpi=150); plt.close(fig)

    if a.no_gr:
        return
    # pooled g(r | s) per (k, R, act) over seeds
    pooled = {}
    for key, ds in sorted(groups.items()):
        O, E = None, None
        for d in ds:
            from defect_spatial import analyse  # profile needed for the expectation
            r = analyse(d, 0.25, 3.0, a.dcut)
            ss, rc, obs, exp, nref = ring_gr(d["frames"], d["R"], r["sc"], r["prof"], raw=True)
            O = obs if O is None else {ab: O[ab] + obs[ab] for ab in O}
            E = exp if E is None else {ab: E[ab] + exp[ab] for ab in E}
        pooled[key] = (ss, rc, O, E, len(ds))
    names = {"++": "+½ / +½", "+-": "+½ / −½", "--": "−½ / −½"}
    n = len(pooled)
    fig, axs = plt.subplots(3, n, figsize=(2.2 * n, 7.0), squeeze=False)
    from matplotlib.colors import LinearSegmentedColormap
    div = LinearSegmentedColormap.from_list("d", ["#2a78d6", "#f2f2f0", "#eb6834"])
    for c, (key, (ss, rc, O, E, ns)) in enumerate(pooled.items()):
        for row, ab in enumerate(("++", "+-", "--")):
            ax = axs[row, c]
            G = np.where(E[ab] > 5, O[ab] / np.maximum(E[ab], 1e-12), np.nan)
            im = ax.pcolormesh(rc, ss, np.log2(G), cmap=div, vmin=-2, vmax=2, shading="nearest")
            ax.set_ylim(key[1], 0); ax.set_xlim(0, 2 * key[1])
            ax.tick_params(colors=INK2, labelsize=7)
            ax.set_title((f"k = {key[0]} ({ns} seeds)\n" if row == 0 else "") + names[ab], fontsize=8, color=INK)
            if row == 2:
                ax.set_xlabel("r", fontsize=8, color=INK2)
            if c == 0:
                ax.set_ylabel("reference ring s", fontsize=8, color=INK2)
    cb = fig.colorbar(im, ax=axs, shrink=0.4, pad=0.01)
    cb.set_label("log₂ g (0 = uncorrelated)", fontsize=8, color=INK2)
    fig.savefig(a.out + "_gr_pooled.png", dpi=150); plt.close(fig)

    # zones pooled over all k (same R): edge s < dcut, outer bulk dcut..R/2, centre s > R/2
    fig, axs = plt.subplots(1, 3, figsize=(10.5, 3.4), sharey=True)
    cols = {"++": PLUS, "+-": AQUA, "--": MINUS}
    R0 = list(pooled)[0][1]
    zones = [("edge ring, s < %g" % a.dcut, 0, a.dcut), ("outer bulk, %g < s < R/2" % a.dcut, a.dcut, R0 / 2), ("centre, s > R/2", R0 / 2, R0)]
    for ax, (zn, s0, s1) in zip(axs, zones):
        style(ax)
        for ab in ("++", "+-", "--"):
            Osum = sum(O[ab][(ss >= s0) & (ss < s1)].sum(0) for (ss, rc, O, E, ns) in pooled.values())
            Esum = sum(E[ab][(ss >= s0) & (ss < s1)].sum(0) for (ss, rc, O, E, ns) in pooled.values())
            g = np.where(Esum > 20, Osum / np.maximum(Esum, 1e-12), np.nan)
            err = np.where(Osum > 0, g / np.sqrt(np.maximum(Osum, 1)), np.nan)
            ax.plot(rc, g, color=cols[ab], lw=1.8, label=names[ab])
            ax.fill_between(rc, g - err, g + err, color=cols[ab], alpha=0.18, lw=0)
        ax.axhline(1, color=INK2, lw=0.7, ls="--")
        ax.set_title(zn, fontsize=9, color=INK); ax.set_xlabel("pair distance r", fontsize=8, color=INK2)
        ax.set_xlim(0, 40)
    axs[0].set_ylabel("g(r | zone of reference)", fontsize=8, color=INK2)
    axs[0].legend(frameon=False, fontsize=8)
    fig.suptitle(f"Pooled over k and seeds (band: Poisson error)", fontsize=9, color=INK)
    fig.tight_layout(); fig.savefig(a.out + "_gr_zones.png", dpi=150); plt.close(fig)
    print("written", a.out + ".csv", a.out + "_sitedist.png", a.out + "_gr_pooled.png", a.out + "_gr_zones.png")


if __name__ == "__main__":
    main()
