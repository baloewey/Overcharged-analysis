#!/usr/bin/env python3
"""Vacancy vs full layer in active circles, and the duck-duck-goose exchange (user, 2026-10-14).
Field of the wall layer felt by one bulk +1/2 (analytic one-constant disk energy with same-sign wall images; layer alone; K = 32.5 kBT), for measured frames of the
k = 6, area-2500 circles with 9 (vacancy), 10 (full) and 11 (extra) wall +1/2. Per frame: range of the field inside the 0.6 level, and whether it points at the
vacancy: angle between the largest gap of the layer (its midpoint) and the direction of the interior field (minus the field's gradient at r = 0.4 R).
Duck-duck-goose (ddg.one on the shape-sweep circles): delay from an entry to the next exit vs random times, and the distance of the 'goose' from the entry site.
-> fig_vacancy_field.png (+pdf), vacancy_field.txt"""
import os, sys, glob, json, numpy as np
from multiprocessing import Pool
from scipy.spatial import cKDTree
R_ = "/home/baloewe/Physics/Overcharged/runs"; H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
sys.path.insert(0, R_ + "/cavity"); sys.path.insert(0, os.path.join(H, "..", "energetics"))
from edge_criteria import read_with_angle  # noqa
import energetics as En  # noqa
PIK = np.pi * 32.5; RUNS = [f"{R_}/fineflow/eps1.0_s1", f"{R_}/fineflow2/eps1.0_s3", f"{R_}/fineflow2/eps1.0_s4"]


def field(layer, R, X):
    """potential (pi K units) of a +1/2 at X from the layer (charges + same-sign images) and its own images removed (layer part only)"""
    zs, qs = layer[:, 0] + 1j * layer[:, 1], 0.5 * np.sign(layer[:, 2]); v = np.zeros(X.shape)
    for z, q in zip(zs, qs):
        v += -2 * 0.5 * q * (np.log(np.abs(X - z)) + np.log(np.abs(R - X * np.conj(z) / R)))
    return v


def frames():
    out = []
    for run in RUNS:
        bc = json.load(open(run + "/input.json"))["BC"][0]; R = bc["R"]; W, s, P = En.wall_polyline(bc); tree = cKDTree(W)
        t, fr = read_with_angle(run, np.array(bc["Q"][:2], float), 0.5)
        for f in fr[::3]:
            if not len(f): continue
            d = tree.query(f[:, :2])[0]; m = (d <= 4) & (f[:, 2] > 0)
            if m.sum() in (9, 10, 11): out.append((int(m.sum()), f[m], R))
    return out


if __name__ == "__main__":
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm
    F = frames(); R = F[0][2]
    xs = np.arange(-R, R + 0.01, 0.75); XX, YY = np.meshgrid(xs, xs, indexing="ij"); X = XX + 1j * YY; ins = np.abs(X) < 0.9 * R; inn = np.abs(X) < 0.6 * R
    rng = np.random.default_rng(1); stats = {9: [], 10: [], 11: []}; examples = {}
    for n, lay, _ in F:
        v = field(lay, R, X) * PIK; rr = np.ptp(v[inn])
        ph = np.sort(np.arctan2(lay[:, 1], lay[:, 0])); gaps = np.diff(np.r_[ph, ph[0] + 2 * np.pi]); g = np.argmax(gaps); gmid = ph[g] + gaps[g] / 2
        # direction of the interior force on a +1/2 (downhill): average of -grad V on the 0.4 R circle
        c = 0.4 * R * np.exp(1j * np.linspace(0, 2 * np.pi, 64, endpoint=False)); h = 0.1
        gx = (field(lay, R, c + h) - field(lay, R, c - h)) / (2 * h); gy = (field(lay, R, c + 1j * h) - field(lay, R, c - 1j * h)) / (2 * h)
        fdir = np.angle(-(gx + 1j * gy).mean()); cosang = np.cos(fdir - gmid)
        stats[n].append((rr, cosang, gaps.max() / (2 * np.pi / 10)))
        if n not in examples or abs(rr - np.median([s[0] for s in stats[n]] or [rr])) < 3: examples[n] = (lay, v)
    L = ["VACANCY vs FULL LAYER, k = 6 circles (area 2500, measured frames): range of the layer's field felt by a bulk +1/2 inside the 0.6 level, and alignment of the",
         "interior force with the largest gap (cos of the angle; 1 = pulls towards the gap). wall count | frames | range median [IQR] kBT | <cos> | fraction cos > 0.7 | largest gap (sites) median"]
    for n in (9, 10, 11):
        S = np.array(stats[n])
        if not len(S): continue
        L.append(f"  {n:2d} | {len(S):4d} | {np.median(S[:,0]):6.1f} [{np.percentile(S[:,0],25):.1f}, {np.percentile(S[:,0],75):.1f}] | {S[:,1].mean():+.2f} | {np.mean(S[:,1] > 0.7):.2f} | {np.median(S[:,2]):.2f}")
    # duck-duck-goose distances (shape-sweep circles)
    import ddg
    runs = sorted(glob.glob(f"{R_}/shapes_server/circle/circle/Req1.[05]/k*/s*"))
    with Pool(6) as pool: D = pool.map(ddg.one, runs)
    P = [p for x in D for p in x["pairs"]]; Nl = np.array([v for x in D for v in x["null"]])
    dl = np.array([p[0] for p in P]); oth = np.array([not p[1] for p in P]); dph = np.array([p[2] for p in P])
    nsite = np.array([(2 * abs(x["k"]) - 2) for x in D for _ in x["pairs"]])
    L.append(f"DUCK-DUCK-GOOSE (shape-sweep circles R_eq 1-1.5, {len(runs)} runs): next exit after an entry, median {np.median(dl):.0f} t.u. (random times {np.median(Nl):.0f}); by another defect {oth.mean():.2f}; "
             f"goose distance from the entry site, median {np.median(dph[oth]):.1f} sites")
    fig = plt.figure(figsize=(18, 9.4)); gs = fig.add_gridspec(2, 4, width_ratios=[1, 1, 1, 1.2]); th = np.linspace(0, 2 * np.pi, 300)
    for i, (n, lab) in enumerate(((10, "full layer (10)"), (9, "one vacancy (9)"), (11, "one extra (11)"))):
        if n not in examples: continue
        lay, v = examples[n]; ax = fig.add_subplot(gs[0, i]); vv = np.where(ins, v - np.nanmin(np.where(inn, v, np.nan)), np.nan)
        hi = max(np.nanpercentile(vv[inn], 98), 5); im = ax.pcolormesh(xs, xs, vv.T, cmap="RdBu_r", shading="nearest", norm=TwoSlopeNorm(0, -hi, hi))
        ax.contour(xs, xs, np.where(inn, vv, np.nan).T, levels=np.arange(5, hi, 5), colors="k", linewidths=0.35, alpha=0.6)
        ax.plot(R * np.cos(th), R * np.sin(th), "k-", lw=1); ax.plot(0.6 * R * np.cos(th), 0.6 * R * np.sin(th), "k--", lw=0.5)
        ax.plot(lay[:, 0], lay[:, 1], "o", ms=5, color="#1f5fa8", mec="w", mew=0.6); ax.set_aspect("equal"); ax.axis("off")
        ax.set_title(f"({'abc'[i]}) Circle, {lab}", fontsize=10); cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02, shrink=0.8); cb.set_label("$V-V_\\mathrm{min}$ in the inner 60% ($k_BT$)", fontsize=8)
    ax = fig.add_subplot(gs[0, 3]); data = [np.array(stats[n])[:, 0] for n in (9, 10, 11) if stats[n]]
    ax.boxplot(data, labels=[f"{n}" for n in (9, 10, 11) if stats[n]], showfliers=False); ax.set_yscale("log"); ax.set_xlabel("wall $+1/2$ in the layer"); ax.set_ylabel("field range, inner 60% ($k_BT$)")
    ax.set_title("(d) Leftover field vs wall count", fontsize=10); [ax.spines[s_].set_visible(False) for s_ in ("top", "right")]
    ax = fig.add_subplot(gs[1, 0]); S9 = np.array(stats[9]); ax.hist(S9[:, 1], bins=np.linspace(-1, 1, 21), color="#1f5fa8", alpha=0.85)
    ax.set_xlabel("cos(force direction − gap direction)"); ax.set_ylabel("frames with 9 wall $+1/2$"); ax.set_title("(e) The bulk is pulled towards the hole", fontsize=10)
    [ax.spines[s_].set_visible(False) for s_ in ("top", "right")]
    ax = fig.add_subplot(gs[1, 1]); b = np.arange(0, 400, 10); ax.hist(dl, bins=b, density=True, histtype="step", lw=1.6, color="#c2452d", label="after an entry")
    ax.hist(Nl, bins=b, density=True, histtype="step", lw=1.6, color="0.5", label="after a random time"); ax.set_xlabel("time to the next exit (t.u.)"); ax.set_ylabel("probability density")
    ax.set_title("(f) An entry triggers an exit", fontsize=10); ax.legend(fontsize=8, frameon=False); [ax.spines[s_].set_visible(False) for s_ in ("top", "right")]
    ax = fig.add_subplot(gs[1, 2]); dd = dph[oth] / (nsite[oth] / 2); ax.hist(dd, bins=np.linspace(0, 1, 11), density=True, color="#7aa36b", alpha=0.85, label="goose (another defect)")
    ax.axhline(1, color="k", ls=":", lw=1, label="uniform"); ax.set_xlabel("goose distance from the entry site\n(fraction of half the ring)"); ax.set_ylabel("probability density")
    ax.set_title("(g) The goose can be anywhere on the ring", fontsize=10); ax.legend(fontsize=8, frameon=False); [ax.spines[s_].set_visible(False) for s_ in ("top", "right")]
    ax = fig.add_subplot(gs[1, 3]); ax.axis("off")
    ax.text(0, 0.95, "Top: field of the wall layer felt\nby one bulk $+1/2$ (one-constant,\nsame-sign wall images, layer\nalone; measured frames, $k=6$).\n\n"
            "A vacancy or an extra wall\ndefect is an unscreened\nhalf-charge whose logarithmic\nfield spans the cavity.\n\n"
            "Bottom: tracked $+1/2$ in the\nshape-sweep circles.", va="top", fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(H, "fig_vacancy_field.png"), dpi=200)
    print("\n".join(L)); open(os.path.join(H, "vacancy_field.txt"), "w").write("\n".join(L) + "\n")
