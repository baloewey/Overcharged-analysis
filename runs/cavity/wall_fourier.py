#!/usr/bin/env python3
"""
Fourier analysis of the defects along the wall (user's idea, 2026-10-02): pinned vs sliding vs no layer.

Wall: circle, Cassini oval (shape 1) or ellipse (aInv), as a fine polyline with cumulative arclength; P = perimeter.
Near-wall +1/2 (distance from the wall < --band cells) are mapped to the arclength s of the nearest wall point.
For each frame n: rho_m(n) = sum_j exp(2 pi i m s_j / P), m = 1..mmax.
  S(m)    = < |rho_m|^2 > / < N >                instantaneous structure factor (peak = ordered spacing, even if sliding)
  f_m     = |< rho_m >|^2 / < |rho_m|^2 >        static fraction (1 = pinned lattice, ~0 = sliding or disordered)
  C_m(t)  = Re< rho_m(n+t) rho_m*(n) > - |< rho_m >|^2, normalised to 1 at t = 0; tau_m = first 1/e crossing
            (collective decay time of the fluctuating part; no tracking needed)
m* = argmax S(m) over m >= 2; reported with the expected m for a lattice of the mean near-wall count.
Also: time-averaged density map of +1/2 and -1/2 (PREFIX_maps.png) and S(m) curves (PREFIX_S.png).

Usage: wall_fourier.py RUN [RUN ...] [--labels ...] [--band 3] [--out PREFIX]
"""
import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from edge_criteria import read_with_angle  # noqa: E402


def wall_polyline(bc, n=6000):
    """Wall points (relative to the centre Q) and cumulative arclength."""
    th = np.linspace(0, 2 * np.pi, n, endpoint=False)
    if bc.get("shape", 0) == 1:
        a, b = bc.get("cassiniA", 0.0), bc["R"]
        r2 = a * a * np.cos(2 * th) + np.sqrt(b ** 4 - a ** 4 * np.sin(2 * th) ** 2)
        x, y = np.sqrt(r2) * np.cos(th), np.sqrt(r2) * np.sin(th)
    else:
        A, B = bc["aInv"][0] * bc["R"], bc["aInv"][1] * bc["R"]     # ellipse semi-axes (circle: A = B = R)
        x, y = A * np.cos(th), B * np.sin(th)
    W = np.column_stack([x, y])
    seg = np.hypot(*np.diff(np.vstack([W, W[:1]]), axis=0).T)
    s = np.concatenate([[0], np.cumsum(seg)[:-1]])
    return W, s, seg.sum()


def project(W, s, xy):
    """-> (distance to the wall, arclength of the nearest wall point) for points xy (n,2)."""
    d = np.empty(len(xy)); sa = np.empty(len(xy))
    for i in range(0, len(xy), 4000):
        blk = xy[i:i + 4000]
        D = ((blk[:, None, :] - W[None, :, :]) ** 2).sum(-1)
        j = D.argmin(1)
        d[i:i + 4000] = np.sqrt(D[np.arange(len(blk)), j]); sa[i:i + 4000] = s[j]
    return d, sa


def analyse(run, band, tfrac, mmax=60):
    j = json.load(open(os.path.join(run, "input.json"))); bc = j["BC"][0]
    t, fr = read_with_angle(run, np.array(bc["Q"][:2]), tfrac)
    dt = t[1] - t[0]
    W, s, P = wall_polyline(bc)
    m = np.arange(1, mmax + 1)
    rho = np.zeros((len(fr), mmax), complex); Nn = np.zeros(len(fr))
    allp, alln = [], []
    for n, f in enumerate(fr):
        p = f[f[:, 2] > 0, :2]
        allp.append(p); alln.append(f[f[:, 2] < 0, :2])
        if not len(p):
            continue
        d, sa = project(W, s, p)
        sel = sa[d < band]
        Nn[n] = len(sel)
        if len(sel):
            rho[n] = np.exp(2j * np.pi * np.outer(sel, m) / P).sum(0)
    mean = rho.mean(0); pw = (np.abs(rho) ** 2).mean(0)
    S = pw / max(Nn.mean(), 1e-9)
    fstat = np.abs(mean) ** 2 / np.maximum(pw, 1e-12)
    ms = 2 + int(np.argmax(S[1:]))                # m >= 2 (S[0] is m = 1)
    i = ms - 1                                    # index of m* in the arrays
    x = rho[:, i] - mean[i]
    var = np.mean(np.abs(x) ** 2)
    L = min(len(x) // 2, int(2000 / dt))
    C = np.array([np.mean((x[l:] * np.conj(x[:len(x) - l])).real) for l in range(L)]) / max(var, 1e-12)
    below = np.where(C < np.exp(-1))[0]
    tau = below[0] * dt if len(below) else np.inf
    return dict(run=run, k=bc.get("customAnch"), P=P, Nwall=Nn.mean(), m=m, S=S, fstat=fstat, mstar=int(ms),
                Sstar=S[i], fstar=fstat[i], tau=tau, C=C, dt=dt, allp=np.vstack(allp) if allp else np.zeros((0, 2)),
                alln=np.vstack([q for q in alln if len(q)]) if any(len(q) for q in alln) else np.zeros((0, 2)), W=W, nf=len(fr))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+"); ap.add_argument("--labels", nargs="+")
    ap.add_argument("--band", type=float, default=4.0); ap.add_argument("--tfrac", type=float, default=0.2)
    ap.add_argument("--out", default="wall_fourier")
    a = ap.parse_args()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    labels = a.labels or [os.path.basename(os.path.normpath(r)) for r in a.runs]
    R = []
    print(f"{'run':<30} {'k':>3} {'P':>6} {'N_wall':>6} {'m*':>4} {'S(m*)':>7} {'f_static':>8} {'tau':>7}")
    for run, lab in zip(a.runs, labels):
        r = analyse(run, a.band, a.tfrac); r["lab"] = lab; R.append(r)
        print(f"{lab:<30} {r['k']:>3} {r['P']:6.0f} {r['Nwall']:6.2f} {r['mstar']:>4d} {r['Sstar']:7.2f} {r['fstar']:8.3f} {r['tau']:7.0f}", flush=True)
    n = len(R)
    fig, axs = plt.subplots(1, n, figsize=(2.9 * n, 2.8), squeeze=False)
    for ax, r in zip(axs[0], R):
        ext = np.abs(r["W"]).max(0) + 1
        for pts, cm in ((r["allp"], "Blues"), (r["alln"], "Oranges")):
            if len(pts):
                H, xe, ye = np.histogram2d(pts[:, 0], pts[:, 1], [np.arange(-ext[0], ext[0] + 0.5, 0.5), np.arange(-ext[1], ext[1] + 0.5, 0.5)])
                H = H / (r["nf"] * 0.25)
                if H.max() > 0:
                    ax.imshow(np.where(H > 0, H, np.nan).T, origin="lower", extent=[xe[0], xe[-1], ye[0], ye[-1]], cmap=cm,
                              norm=LogNorm(vmin=H.max() * 1e-3, vmax=H.max()), interpolation="nearest", alpha=0.9)
        ax.plot(np.r_[r["W"][:, 0], r["W"][0, 0]], np.r_[r["W"][:, 1], r["W"][0, 1]], color="#2b2f36", lw=0.8)
        ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.set_title(f"{r['lab']}\nm*={r['mstar']}  S={r['Sstar']:.1f}  f={r['fstar']:.2f}", fontsize=7)
    fig.tight_layout(); fig.savefig(a.out + "_maps.png", dpi=150); plt.close(fig)
    fig, ax = plt.subplots(figsize=(6.5, 3.8))
    for r in R:
        ax.plot(r["m"], r["S"], lw=1.4, label=r["lab"])
    ax.axhline(1, color="#555a62", lw=0.7, ls=":"); ax.set_xlabel("mode m along the wall (wavelength P/m)"); ax.set_ylabel("S(m)")
    ax.set_yscale("log"); ax.legend(frameon=False, fontsize=6); ax.grid(color="#e9ebee", lw=0.6)
    fig.tight_layout(); fig.savefig(a.out + "_S.png", dpi=150)
    print("written", a.out + "_maps.png", a.out + "_S.png")


if __name__ == "__main__":
    main()
