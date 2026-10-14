#!/usr/bin/env python3
"""Lifetimes of dynamical defect pairs and who annihilates with whom (user, 2026-10-14).
Tracking: kickout/kick_identity.track (same-sign links < 2.5 cells per frame, gaps of up to 3 frames bridged). A PAIR birth = a new +1/2 and a new -1/2 within 3 cells
in the same frame, away from the wall (depth >= 4). Each pair's -1/2 is followed to its death (last frame seen; deaths in the final 50 frames are censored).
Death classes of the -1/2 (an opposite-sign track ending within 3 cells, within +-1 frame):
  PARTNER  - annihilates with its own birth partner;
  BULK+    - with another bulk +1/2 (an excess / older +1/2: the partner survives and takes its place, a 'charge relay');
  WALL+    - with a +1/2 in the wall band (depth < 4);
  LOST     - no opposite-sign death nearby (tracking loss; caveat).
Also: lifetime of the partner +1/2 when the -1/2 died with someone else (how long the relayed +1/2 lives on).
Geometries: periodic box (neutral control), cavities from the shape sweep, winding channels. Frames: second half, every frame (2 t.u.).
Usage: pair_lifetime.py -> pair_lifetime.txt, fig_pair_lifetime.png (+pdf)"""
import os, sys, glob, json, numpy as np
from multiprocessing import Pool
from scipy.spatial import cKDTree
R_ = "/home/baloewe/Physics/Overcharged/runs"; H = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, R_ + "/cavity"); sys.path.insert(0, os.path.join(H, "..", "kickout")); sys.path.insert(0, os.path.join(H, "..", "energetics"))
from edge_criteria import read_with_angle  # noqa
from kick_identity import track  # noqa
import energetics as En  # noqa


def depth_fn(run):
    J = json.load(open(run + "/input.json")); BC = J.get("BC", []); Lx = J["domain"][0]
    if not BC: return np.zeros(2), (lambda p: np.full(len(p), 99.0)), Lx, "box"
    if len(BC) >= 4 and BC[2].get("eMode", 0) == 4:
        y0, y1 = BC[2]["Q"][1], BC[3]["Q"][1]; return np.zeros(2), (lambda p: np.minimum(p[:, 1] - y0, y1 - p[:, 1])), Lx, "channel"
    W, s, P = En.wall_polyline(BC[0]); tree = cKDTree(W); return np.array(BC[0]["Q"][:2], float), (lambda p: tree.query(p)[0]), None, "cavity"


def one(run):
    Q, dep, px, kind = depth_fn(run); t, fr = read_with_angle(run, Q, 0.5); t = np.array(t); dt = float(np.median(np.diff(t)))
    if px:   # periodic in x: the tracker uses plain distances, so keep tracks across the seam by unwrapping is not needed for lifetimes > a few frames; ignore seam pairs
        pass
    D = [dep(f[:, :2]) if len(f) else np.zeros(0) for f in fr]
    ids, info = track(fr, D, t)
    last = {}; pos = {}
    for n, (f, idn) in enumerate(zip(fr, ids)):
        for i, k in enumerate(idn): last[k] = n; pos[(k, n)] = (f[i, 0], f[i, 1], D[n][i])
    births = {}                                            # minus id -> plus partner id
    for n, (f, idn) in enumerate(zip(fr, ids)):
        if n == 0: continue
        newm = [i for i, k in enumerate(idn) if info[k]["birth"] == t[n] and f[i, 2] < 0 and D[n][i] >= 4]
        newp = [i for i, k in enumerate(idn) if info[k]["birth"] == t[n] and f[i, 2] > 0]
        for i in newm:
            dd = [np.hypot(*(f[j, :2] - f[i, :2])) for j in newp]
            if dd and min(dd) < 3: births[idn[i]] = idn[newp[int(np.argmin(dd))]]
    deaths_at = {}
    for k, n in last.items(): deaths_at.setdefault(n, []).append(k)
    N = len(fr); out = []
    for m, p in births.items():
        nb = int(np.flatnonzero(t == info[m]["birth"])[0]); nd = last[m]
        if nd >= N - 50: out.append(("CENSORED", (nd - nb) * dt, np.nan)); continue
        xm, ym, _ = pos[(m, nd)]; cls, plife = "LOST", np.nan
        cands = [k for dn in (nd - 1, nd, nd + 1) for k in deaths_at.get(dn, []) if info[k]["sign"] > 0 and (k, last[k]) in pos]
        best = None
        for k in cands:
            xk, yk, dk = pos[(k, last[k])]; d = np.hypot(xk - xm, yk - ym)
            if px: d = np.hypot(((xk - xm) + px / 2) % px - px / 2, yk - ym)
            if d < 3 and (best is None or d < best[0]): best = (d, k, dk)
        if best is not None:
            k = best[1]; cls = "PARTNER" if k == p else ("WALL+" if best[2] < 4 else "BULK+")
            if cls != "PARTNER" and p in last: plife = (last[p] - nb) * dt if last[p] < N - 50 else np.nan
        out.append((cls, (nd - nb) * dt, plife))
    Ab = None
    return dict(run=os.path.relpath(run, R_), kind=kind, ev=out, T=N * dt)


GROUPS = {
    "periodic box, neutral": sorted(glob.glob(f"{R_}/activeLength/L200/act_0.070*")),
    "circles R_eq 1.5-2": sorted(glob.glob(f"{R_}/shapes_server/circle/circle/Req1.5/k*/s*") + glob.glob(f"{R_}/shapes_server/circle/circle/Req2.0/k*/s*")),
    "ellipses 2.0-2.5, R_eq 1.5-2": sorted(glob.glob(f"{R_}/shapes_server/ellipse/eps2.[05]0/Req1.5/k*/s*") + glob.glob(f"{R_}/shapes_server/ellipse/eps2.[05]0/Req2.0/k*/s*")),
    "channels W 40-60": sorted(glob.glob(f"{R_}/channel/W40_m6_L192_s*") + glob.glob(f"{R_}/channel/W60_m6_L192_s*")),
}

if __name__ == "__main__":
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    runs = [r for v in GROUPS.values() for r in v]
    with Pool(int(os.environ.get("NPROC", 6))) as pool: res = dict(zip(runs, pool.map(one, runs, chunksize=1)))
    L = ["PAIR LIFETIMES (tracked; pair births in the bulk; -1/2 followed to death; second half of each run). Statistics for ESTABLISHED pairs (-1/2 seen >= TMIN t.u.); shorter ones = flicker.",
         "group | runs | pairs | -1/2 death: partner / other bulk +1/2 / wall +1/2 / lost / censored (fractions) | median -1/2 lifetime by class (t.u.) partner / bulk+ / wall+ | mean (all, uncensored) | relayed +1/2 lives on (median t.u.)"]
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2)); cols = ["0.4", "#1f5fa8", "#c2452d", "#7aa36b"]
    for gi, (name, rr) in enumerate(GROUPS.items()):
        ev = [e for r in rr for e in res[r]["ev"]]
        if not ev: continue
        cl0 = np.array([e[0] for e in ev]); lt0 = np.array([e[1] for e in ev]); pl0 = np.array([e[2] for e in ev])
        TMIN = float(os.environ.get("TMIN", 10)); est = lt0 >= TMIN; flick = np.mean(~est)       # established pairs: -1/2 seen for >= TMIN t.u.
        cl, lt, pl = cl0[est], lt0[est], pl0[est]
        rate_all = len(ev) / sum(res[r]["T"] for r in rr); rate_est = est.sum() / sum(res[r]["T"] for r in rr)
        fr_ = [np.mean(cl == c) for c in ("PARTNER", "BULK+", "WALL+", "LOST", "CENSORED")]
        med = [np.median(lt[cl == c]) if (cl == c).any() else np.nan for c in ("PARTNER", "BULK+", "WALL+")]
        unc = cl != "CENSORED"
        L.append(f"{name:28s} | {len(rr):2d} | {len(ev):6d} | " + " / ".join(f"{x:.2f}" for x in fr_) + " | " + " / ".join(f"{x:.0f}" for x in med)
                 + f" | {lt[unc].mean():.0f} | {np.nanmedian(pl):.0f} | flicker (< {TMIN:.0f} t.u.) {flick:.2f}; established pairs per 1000 t.u. per run {1000*rate_est:.1f}")
        x = np.sort(lt[unc]); ax[0].plot(x, 1 - np.arange(len(x)) / len(x), color=cols[gi], label=name)
        ax[1].bar(np.arange(4) + 0.2 * gi - 0.3, fr_[:4], width=0.2, color=cols[gi], label=name)
    ax[0].set_yscale("log"); ax[0].set_xlim(0, 400); ax[0].set_xlabel("lifetime of a pair-born $-1/2$ (t.u., established pairs)"); ax[0].set_ylabel("fraction surviving longer")
    ax[0].set_title("(a) survival of dynamical $-1/2$", fontsize=10)
    ax[1].set_xticks(range(4)); ax[1].set_xticklabels(["own partner", "other bulk $+1/2$", "wall $+1/2$", "lost (tracking)"]); ax[1].set_ylabel("fraction of $-1/2$ deaths")
    ax[1].set_title("(b) who the $-1/2$ annihilates with", fontsize=10)
    for a_ in ax: [a_.spines[s_].set_visible(False) for s_ in ("top", "right")]
    fig.legend(*ax[0].get_legend_handles_labels(), fontsize=8, frameon=False, ncol=4, loc="lower center"); fig.tight_layout(rect=(0, 0.1, 1, 1))
    fig.savefig(os.path.join(H, "fig_pair_lifetime.png"), dpi=200)
    print("\n".join(L)); open(os.path.join(H, "pair_lifetime.txt"), "w").write("\n".join(L) + "\n")
