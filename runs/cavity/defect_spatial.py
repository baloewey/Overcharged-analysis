#!/usr/bin/env python3
"""
Spatial distribution of defects in circular cavities (edge/bulk distinction).

Per run directory (input.json + defects.dat[.gz]), after dropping the first --tfrac of frames:
  2D maps       time-averaged number density of +1/2 and -1/2 defects (per cell^2 per frame)
  radial        n_±(s): number density vs distance from the wall s = R - r (per unit area)
  speed         mean speed of defects vs s: each defect is matched to the nearest same-sign defect
                in the next frame (if closer than --rlink); speed = displacement / dt
  ring g        g_ab(s, dphi): same-ring angular pair correlation in rings of width 1 (see ring_gphi)
  ring g(r)     g_ab(r | s): reference in ring s, partner anywhere, radial profile divided out (see ring_gr)
  RDF           g_ab(r) for ++, --, +- pairs: observed pair-distance histogram divided by the
                reference histogram for the same one-body density without correlations
                (pairs between decorrelated frames, i and i + nf/2; keeps the radial profile); computed for all defects and for
                bulk defects only (s > --dbulk, disc of radius R - dbulk)

Usage:
  defect_spatial.py RUN [RUN ...] [--labels A B ...] [--tfrac 0.2] [--dbulk 4.2] [--out PREFIX]
Writes PREFIX_maps.png, PREFIX_radial.png, PREFIX_rdf.png and PREFIX_profiles.csv.
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from defect_audit import read_frames  # noqa: E402

PLUS, MINUS = "#2a78d6", "#eb6834"      # categorical slots 1, 2
INK, INK2, GRID = "#15171b", "#555a62", "#e9ebee"


def load(run, tfrac):
    j = json.load(open(os.path.join(run, "input.json")))
    bc = j["BC"][-1]
    R, C = bc["R"], np.array(bc["Q"][:2], float)
    dt = j.get("defectsOut", 50) * j.get("dt", 0.1)
    t, frames = read_frames(os.path.join(run, "defects.dat"))
    frames = frames[int(tfrac * len(frames)):]
    frames = [np.column_stack([f[:, :2] - C, f[:, 2]]) if len(f) else np.zeros((0, 3)) for f in frames]
    act = j["species"][0].get("act", 0.0)
    return dict(run=run, R=R, k=bc.get("customAnch"), act=act, dt=dt, frames=frames)


def disc_pair_pdf(r, R):
    """pdf of the distance between two uniform points in a disc of radius R."""
    x = np.clip(r / (2 * R), 0, 1)
    return 4 * r / (np.pi * R ** 2) * (np.arccos(x) - x * np.sqrt(1 - x ** 2))


def _pairs(a, b, same):
    if same:
        n = len(a)
        if n < 2:
            return np.zeros(0)
        return np.hypot(*(a[:, None, :] - a[None, :, :]).transpose(2, 0, 1))[np.triu_indices(n, 1)]
    if len(a) == 0 or len(b) == 0:
        return np.zeros(0)
    return np.hypot(*(a[:, None, :] - b[None, :, :]).transpose(2, 0, 1)).ravel()


def rdf(frames, R, rbins, smin=None):
    """g_ab(r) for ++, --, +- in the inhomogeneous disc.

    Observed pair distances within a frame, divided by the pair distances between frame i and
    frame i + nf/2 (decorrelated), which keeps the one-body density profile n_a(x) but removes
    correlations: g = 1 means no correlation beyond the radial profile. Same-sign cross-frame
    pairs are scaled by N(N-1)/2 / (N_i N_j) to count pairs like the within-frame ones.
    Restricted to s = R - r > smin if given."""
    Rd = R if smin is None else R - smin
    fr = [f[np.hypot(f[:, 0], f[:, 1]) < Rd] if smin is not None else f for f in frames]
    nf, lag = len(fr), len(fr) // 2
    obs = {key: np.zeros(len(rbins) - 1) for key in ("++", "--", "+-")}
    ref = {key: np.zeros(len(rbins) - 1) for key in obs}
    for i, f in enumerate(fr):
        g = fr[(i + lag) % nf]
        P, M = f[f[:, 2] > 0, :2], f[f[:, 2] < 0, :2]
        P2, M2 = g[g[:, 2] > 0, :2], g[g[:, 2] < 0, :2]
        for key, a, b, a2, b2 in (("++", P, P, P, P2), ("--", M, M, M, M2), ("+-", P, M, P, M2)):
            same = key[0] == key[1]
            obs[key] += np.histogram(_pairs(a, b, same), rbins)[0]
            d2 = _pairs(a2, b2, False)
            if len(d2):
                w = (len(a) * (len(a) - 1) / 2) / (len(a2) * len(b2)) if same else len(a) * len(b) / (len(a2) * len(b2))
                ref[key] += w * np.histogram(d2, rbins)[0]
    rc = 0.5 * (rbins[1:] + rbins[:-1])
    return rc, {key: np.where(ref[key] > 5, obs[key] / np.maximum(ref[key], 1e-12), np.nan) for key in obs}


def layer_sites(frames, R, dlayer, nbin=360, frac=0.15):
    """Number of fixed angular sites of +1/2 defects with s < dlayer (peaks of the angle histogram)."""
    th = np.concatenate([np.arctan2(f[f[:, 2] > 0, 1], f[f[:, 2] > 0, 0])[R - np.hypot(f[f[:, 2] > 0, 0], f[f[:, 2] > 0, 1]) < dlayer] for f in frames])
    h = np.histogram(th, nbin, (-np.pi, np.pi))[0].astype(float)
    ker = np.exp(-0.5 * (np.arange(-6, 7) / 2.0) ** 2)
    hs = np.convolve(np.r_[h[-6:], h, h[:6]], ker / ker.sum(), "valid")
    above = hs > frac * hs.max()
    # count runs of 'above' on the periodic ring
    starts = above & ~np.roll(above, 1)
    return int(starts.sum()) if not above.all() else 0, len(th) / len(frames)


def ring_gphi(frames, R, width=1.0, nphi=90):
    """Same-ring angular pair correlation g_ab(s, dphi) for ++, --, +-.

    The cavity is split into rings of the given width in s = R - r. For each ring, pairs of
    defects in the same ring and frame are histogrammed in the angle difference dphi in [0, pi];
    dividing by the uniform expectation (pairs * dphi_bin / pi) gives g = 1 for no angular
    correlation. Needs only within-frame pairs and no model of the radial profile."""
    nr = int(np.ceil(R / width))
    pb = np.linspace(0, np.pi, nphi + 1)
    obs = {key: np.zeros((nr, nphi)) for key in ("++", "--", "+-")}
    npair = {key: np.zeros(nr) for key in obs}
    for f in frames:
        if len(f) < 2:
            continue
        ring = np.minimum(((R - np.hypot(f[:, 0], f[:, 1])) / width).astype(int).clip(0), nr - 1)
        phi = np.arctan2(f[:, 1], f[:, 0])
        for i in np.unique(ring):
            m = ring == i
            q, ph = f[m, 2], phi[m]
            P, M = ph[q > 0], ph[q < 0]
            for key, a, b in (("++", P, P), ("--", M, M), ("+-", P, M)):
                if key[0] == key[1]:
                    if len(a) < 2:
                        continue
                    d = np.abs(a[:, None] - a[None, :])[np.triu_indices(len(a), 1)]
                else:
                    if len(a) == 0 or len(b) == 0:
                        continue
                    d = np.abs(a[:, None] - b[None, :]).ravel()
                d = np.minimum(d, 2 * np.pi - d)
                obs[key][i] += np.histogram(d, pb)[0]
                npair[key][i] += len(d)
    g = {key: np.where(npair[key][:, None] >= 30, obs[key] / np.maximum(npair[key][:, None], 1) * nphi, np.nan) for key in obs}
    return (np.arange(nr) + 0.5) * width, 0.5 * (pb[1:] + pb[:-1]), g, npair


def ring_gr(frames, R, prof_s, prof, width=1.0, dr=1.0, ntheta=720, raw=False):
    """g_ab(r | s): pair correlation vs distance r, for reference defects of type a in the ring s
    (width `width`) and partners of type b anywhere in the cavity.

    Expected count for no correlation = integral of the measured radial density n_b(s') over the
    annulus [r, r+dr] around the reference (restricted to the cavity), so the radial profile is
    divided out and only within-frame pairs are used (self excluded). a, b in {+, -}."""
    nr = int(np.ceil(R / width))
    rb = np.arange(0, 2 * R + dr, dr)
    rc = 0.5 * (rb[1:] + rb[:-1])
    th = (np.arange(ntheta) + 0.5) * 2 * np.pi / ntheta
    ds = prof_s[1] - prof_s[0]
    def nb_at(key, rad):
        sv = R - rad
        i = np.clip((sv / ds).astype(int), 0, len(prof_s) - 1)
        return np.where(rad < R, prof[key][i], 0.0)
    # expected partner number per reference in each (ring, r) bin, per frame
    expd = {}
    for key in ("+", "-"):
        E = np.zeros((nr, len(rc)))
        for i in range(nr):
            r0 = R - (i + 0.5) * width
            x = r0 + rc[:, None] * np.cos(th)[None, :]
            y = rc[:, None] * np.sin(th)[None, :]
            E[i] = nb_at(key, np.hypot(x, y)).mean(1) * 2 * np.pi * rc * dr
        expd[key] = E
    obs = {ab: np.zeros((nr, len(rc))) for ab in ("++", "+-", "-+", "--")}
    nref = {a: np.zeros(nr) for a in ("+", "-")}
    for f in frames:
        if len(f) == 0:
            continue
        ring = np.minimum(((R - np.hypot(f[:, 0], f[:, 1])) / width).astype(int).clip(0), nr - 1)
        sg = np.where(f[:, 2] > 0, "+", "-")
        D = np.hypot(*(f[:, None, :2] - f[None, :, :2]).transpose(2, 0, 1))
        np.fill_diagonal(D, np.inf)
        for a in ("+", "-"):
            ia = np.where(sg == a)[0]
            np.add.at(nref[a], ring[ia], 1)
            for b in ("+", "-"):
                ib = sg == b
                for i in ia:
                    obs[a + b][ring[i]] += np.histogram(D[i, ib], rb)[0]
    g = {}
    for ab in obs:
        e = expd[ab[1]] * nref[ab[0]][:, None]
        g[ab] = np.where((e > 3) & (nref[ab[0]][:, None] >= 30), obs[ab] / np.maximum(e, 1e-12), np.nan)
    if raw:
        return (np.arange(nr) + 0.5) * width, rc, obs, {ab: expd[ab[1]] * nref[ab[0]][:, None] for ab in obs}, nref
    return (np.arange(nr) + 0.5) * width, rc, g, nref


def analyse(d, sbin, rlink, dbulk):
    R, fr = d["R"], d["frames"]
    nf = len(fr)
    allf = np.vstack([f for f in fr if len(f)])
    # 2D maps
    edges = np.arange(-R - 0.5, R + 0.51, 0.5)
    maps = {}
    for sgn, key in ((1, "+"), (-1, "-")):
        sel = allf[np.sign(allf[:, 2]) == sgn]
        H, _, _ = np.histogram2d(sel[:, 0], sel[:, 1], [edges, edges])
        maps[key] = H / (nf * 0.25)
    # radial density vs s
    sb = np.arange(0, R + sbin, sbin)
    rin, rout = np.clip(R - sb[1:], 0, None), R - sb[:-1]
    area = np.pi * (rout ** 2 - rin ** 2)
    prof = {}
    for sgn, key in ((1, "+"), (-1, "-")):
        sel = allf[np.sign(allf[:, 2]) == sgn]
        s = R - np.hypot(sel[:, 0], sel[:, 1])
        prof[key] = np.histogram(s, sb)[0] / (area * nf)
    # speed vs s (nearest same-sign defect in the next frame)
    sp_s, sp_v, sp_q = [], [], []
    for f0, f1 in zip(fr[:-1], fr[1:]):
        for sgn in (1, -1):
            a, b = f0[np.sign(f0[:, 2]) == sgn, :2], f1[np.sign(f1[:, 2]) == sgn, :2]
            if len(a) == 0 or len(b) == 0:
                continue
            dd = np.hypot(*(a[:, None, :] - b[None, :, :]).transpose(2, 0, 1))
            j = dd.argmin(1)
            ok = dd[np.arange(len(a)), j] < rlink
            sp_s += list(R - np.hypot(a[ok, 0], a[ok, 1]))
            sp_v += list(dd[np.arange(len(a)), j][ok] / d["dt"])
            sp_q += [sgn] * int(ok.sum())
    sp_s, sp_v, sp_q = map(np.array, (sp_s, sp_v, sp_q))
    speed = {}
    for sgn, key in ((1, "+"), (-1, "-")):
        m = sp_q == sgn
        idx = np.digitize(sp_s[m], sb) - 1
        speed[key] = np.array([sp_v[m][idx == i].mean() if (idx == i).sum() >= 20 else np.nan for i in range(len(sb) - 1)])
    nsite, nlayer = layer_sites(fr, R, dbulk)
    ring_s, ring_phi, g_ring, ring_npair = ring_gphi(fr, R)
    gr_s, gr_r, g_rs, gr_nref = ring_gr(fr, R, 0.5 * (sb[1:] + sb[:-1]), prof)
    rbins = np.arange(0, 2 * R + 1.0, 1.0)
    rc, g_all = rdf(fr, R, rbins)
    _, g_bulk = rdf(fr, R, rbins, smin=dbulk)
    return dict(edges=edges, maps=maps, sc=0.5 * (sb[1:] + sb[:-1]), prof=prof, speed=speed,
                rc=rc, g_all=g_all, g_bulk=g_bulk, nf=nf, nsite=nsite, nlayer=nlayer, ring_s=ring_s, ring_phi=ring_phi, g_ring=g_ring, ring_npair=ring_npair, gr_s=gr_s, gr_r=gr_r, g_rs=g_rs, gr_nref=gr_nref)


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
    ap.add_argument("--labels", nargs="+")
    ap.add_argument("--tfrac", type=float, default=0.2)
    ap.add_argument("--sbin", type=float, default=0.25)
    ap.add_argument("--rlink", type=float, default=3.0)
    ap.add_argument("--dbulk", type=float, default=4.2)
    ap.add_argument("--out", default="defect_spatial")
    a = ap.parse_args()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, LogNorm

    labels = a.labels or [os.path.normpath(r) for r in a.runs]
    D = [load(r, a.tfrac) for r in a.runs]
    A = [analyse(d, a.sbin, a.rlink, a.dbulk) for d in D]
    n = len(D)
    cm = {"+": LinearSegmentedColormap.from_list("p", ["#ffffff", "#cde2fb", "#6da7ec", "#2a78d6", "#104281"]),
          "-": LinearSegmentedColormap.from_list("m", ["#ffffff", "#fde0d2", "#f4a27f", "#eb6834", "#8f3413"])}

    # maps
    fig, axs = plt.subplots(2, n, figsize=(2.6 * n, 5.6), squeeze=False)
    for c, (d, r, lab) in enumerate(zip(D, A, labels)):
        for row, key in enumerate(("+", "-")):
            ax = axs[row, c]
            H = r["maps"][key].T
            vmax = max(H.max(), 1e-6)
            ax.imshow(np.where(H > 0, H, np.nan), origin="lower", extent=[r["edges"][0], r["edges"][-1]] * 2,
                      cmap=cm[key], norm=LogNorm(vmin=vmax * 1e-3, vmax=vmax), interpolation="nearest")
            th = np.linspace(0, 2 * np.pi, 400)
            ax.plot(d["R"] * np.cos(th), d["R"] * np.sin(th), color=INK2, lw=0.8)
            ax.plot((d["R"] - a.dbulk) * np.cos(th), (d["R"] - a.dbulk) * np.sin(th), color=INK2, lw=0.6, ls=":")
            ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_visible(False)
            if row == 0:
                ax.set_title(lab, fontsize=9, color=INK)
        axs[0, 0].set_ylabel("+½ defects", color=INK, fontsize=9)
        axs[1, 0].set_ylabel("−½ defects", color=INK, fontsize=9)
    fig.suptitle(f"Time-averaged defect density (log scale, 3 decades; dotted: s = {a.dbulk:g})", fontsize=10, color=INK)
    fig.tight_layout(); fig.savefig(a.out + "_maps.png", dpi=150); plt.close(fig)

    # radial density + speed (two rows, shared x per column)
    fig, axs = plt.subplots(2, n, figsize=(2.8 * n, 5.4), squeeze=False, sharex="col")
    for c, (d, r, lab) in enumerate(zip(D, A, labels)):
        ax = axs[0, c]; style(ax)
        ax.plot(r["sc"], r["prof"]["+"], color=PLUS, lw=1.8, label="+½")
        ax.plot(r["sc"], r["prof"]["-"], color=MINUS, lw=1.8, label="−½")
        ax.set_yscale("log"); ax.set_ylim(1e-5, None)
        ax.set_title(lab, fontsize=9, color=INK)
        ax2 = axs[1, c]; style(ax2)
        ax2.plot(r["sc"], r["speed"]["+"], color=PLUS, lw=1.8)
        ax2.plot(r["sc"], r["speed"]["-"], color=MINUS, lw=1.8)
        ax2.set_xlabel("distance from wall s", fontsize=8, color=INK2)
        for x in (ax, ax2):
            x.axvline(a.dbulk, color=INK2, lw=0.8, ls=":")
    axs[0, 0].set_ylabel("number density (per cell², per frame)", fontsize=8, color=INK2)
    axs[1, 0].set_ylabel("mean speed (cells / time)", fontsize=8, color=INK2)
    axs[0, 0].legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(a.out + "_radial.png", dpi=150); plt.close(fig)

    # RDF: all defects (top), bulk only (bottom)
    fig, axs = plt.subplots(2, n, figsize=(2.8 * n, 5.4), squeeze=False)
    cols = {"++": PLUS, "--": MINUS, "+-": "#1baf7a"}
    for c, (d, r, lab) in enumerate(zip(D, A, labels)):
        for row, (g, name) in enumerate(((r["g_all"], "all"), (r["g_bulk"], f"bulk, s > {a.dbulk:g}"))):
            ax = axs[row, c]; style(ax)
            for key in ("++", "--", "+-"):
                ax.plot(r["rc"], g[key], color=cols[key], lw=1.6, label={"++": "+½ / +½", "--": "−½ / −½", "+-": "+½ / −½"}[key])
            ax.axhline(1, color=INK2, lw=0.6, ls="--")
            ax.set_title(f"{lab}\n{name}" if row == 0 else name, fontsize=8, color=INK)
            ax.set_xlim(0, 2 * (d["R"] if row == 0 else d["R"] - a.dbulk))
            if row == 1:
                ax.set_xlabel("pair distance r", fontsize=8, color=INK2)
    axs[0, 0].set_ylabel("g(r) (1 = uniform in disc)", fontsize=8, color=INK2)
    axs[1, 0].set_ylabel("g(r)", fontsize=8, color=INK2)
    axs[0, 0].legend(frameon=False, fontsize=7)
    fig.tight_layout(); fig.savefig(a.out + "_rdf.png", dpi=150); plt.close(fig)

    # same-ring angular correlation: rows ++, +-, --; x = dphi / pi, y = s (ring of width 1)
    fig, axs = plt.subplots(3, n, figsize=(2.6 * n, 7.4), squeeze=False)
    div = LinearSegmentedColormap.from_list("d", ["#2a78d6", "#f2f2f0", "#eb6834"])
    names = {"++": "+½ / +½", "+-": "+½ / −½", "--": "−½ / −½"}
    for c, (d, r, lab) in enumerate(zip(D, A, labels)):
        for row, key in enumerate(("++", "+-", "--")):
            ax = axs[row, c]
            G = r["g_ring"][key]
            im = ax.pcolormesh(r["ring_phi"] / np.pi, r["ring_s"], np.log2(G), cmap=div, vmin=-2, vmax=2, shading="nearest")
            ax.set_ylim(d["R"], 0)
            ax.tick_params(colors=INK2, labelsize=7)
            for sp in ax.spines.values():
                sp.set_color(INK2)
            if d["k"] and d["k"] > 1 and key == "++":
                for m in range(1, d["k"]):
                    ax.axvline(m / (d["k"] - 1), color=INK2, lw=0.4, ls=":")
            ax.set_title(f"{lab}\n{names[key]}" if row == 0 else names[key], fontsize=8, color=INK)
            if row == 2:
                ax.set_xlabel("Δφ / π", fontsize=8, color=INK2)
            if c == 0:
                ax.set_ylabel("ring: distance from wall s", fontsize=8, color=INK2)
    cb = fig.colorbar(im, ax=axs, shrink=0.5, pad=0.01)
    cb.set_label("log₂ g (0 = uncorrelated; blank: < 30 pairs)", fontsize=8, color=INK2)
    cb.ax.tick_params(labelsize=7, colors=INK2)
    fig.savefig(a.out + "_ring_gphi.png", dpi=150); plt.close(fig)

    # g(r | ring of the reference defect)
    fig, axs = plt.subplots(4, n, figsize=(2.6 * n, 9.6), squeeze=False)
    names = {"++": "ref +½, partner +½", "+-": "ref +½, partner −½", "-+": "ref −½, partner +½", "--": "ref −½, partner −½"}
    for c, (d, r, lab) in enumerate(zip(D, A, labels)):
        for row, key in enumerate(("++", "+-", "-+", "--")):
            ax = axs[row, c]
            im = ax.pcolormesh(r["gr_r"], r["gr_s"], np.log2(r["g_rs"][key]), cmap=div, vmin=-2, vmax=2, shading="nearest")
            ax.set_ylim(d["R"], 0); ax.set_xlim(0, 2 * d["R"])
            ax.tick_params(colors=INK2, labelsize=7)
            for sp in ax.spines.values():
                sp.set_color(INK2)
            ax.set_title(f"{lab}\n{names[key]}" if row == 0 else names[key], fontsize=8, color=INK)
            if row == 3:
                ax.set_xlabel("pair distance r", fontsize=8, color=INK2)
            if c == 0:
                ax.set_ylabel("reference ring s", fontsize=8, color=INK2)
    cb = fig.colorbar(im, ax=axs, shrink=0.4, pad=0.01)
    cb.set_label("log₂ g (0 = uncorrelated; blank: < 30 refs or < 3 expected)", fontsize=8, color=INK2)
    cb.ax.tick_params(labelsize=7, colors=INK2)
    fig.savefig(a.out + "_ring_gr.png", dpi=150); plt.close(fig)

    with open(a.out + "_profiles.csv", "w") as fo:
        fo.write("run,label,R,k,act,frames,s,n_plus,n_minus,v_plus,v_minus\n")
        for d, r, lab in zip(D, A, labels):
            for i, s in enumerate(r["sc"]):
                fo.write(f"{d['run']},{lab},{d['R']},{d['k']},{d['act']},{r['nf']},{s:.3f},"
                         f"{r['prof']['+'][i]:.6g},{r['prof']['-'][i]:.6g},{r['speed']['+'][i]:.4g},{r['speed']['-'][i]:.4g}\n")
    for d, r, lab in zip(D, A, labels):
        print(f"{lab:<16} k={d['k']}  layer +1/2 per frame (s < {a.dbulk:g}) = {r['nlayer']:.2f}   fixed angular sites = {r['nsite']}")
    print("written", ", ".join(a.out + x for x in ("_maps.png", "_radial.png", "_rdf.png", "_ring_gphi.png", "_ring_gr.png", "_profiles.csv")))


if __name__ == "__main__":
    main()
