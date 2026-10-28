#!/usr/bin/env python3
"""Continuum K test analysis (2026-10-28; runs/continuum_Ktest, see run.sh there). Aspect-2 ellipse R_eq 31.83 (grid units = MPCD cells), zeta 1, lambda 2.
Defects every 2 t.u. (<run>_defects.npz), t >= TSKIP. +1/2 tracked with the validated tracker (trackers.track_new, gate 9, 2-frame gaps), the same
definitions as the MPCD shape sweep:
  p_esc   wall +1/2 (d < 4) at a start frame (every 10 t.u.) found deeper than 6 cells 20 t.u. later (per wall defect per 20 t.u.); p_esc40 the same at 40 t.u.
  S(x)    residents (d < 4 for the previous 20 t.u.): fraction whose deepest point in the next 20 t.u. exceeds x
  N_wall  mean wall +1/2; Q_b bulk charge (d > 6); coh / |circ| from the 100 t.u. observations (<run>.json)
Conditions: ref; K2 (K x2 alone); K2G05 (K and gamma_1 x2); N05 (noise variance /2 = K, zeta, a, c x2 up to a time rescale x2: its 20 t.u. window
corresponds to 10 t.u. of the rescaled system, so p_esc40 of N05 is the like-for-like comparison with... see the note in the output).
-> ktest_escape.csv, ktest_escape.txt"""
import csv, glob, json, os, sys, numpy as np
from scipy.spatial import cKDTree
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(H, "..", "crowd_protect", "tracker_validation"))
from trackers import track_new  # noqa
ROOT = "/home/baloewe/Physics/Overcharged/runs/continuum_Ktest"; REQ = 31.830989; EPS = 2.0; TSKIP = 500.0; DT = 2.0
A_, B_ = REQ * np.sqrt(EPS), REQ / np.sqrt(EPS); th = np.linspace(0, 2 * np.pi, 8001)[:-1]; TREE = cKDTree(np.column_stack([A_ * np.cos(th), B_ * np.sin(th)]))
XG = np.arange(1.5, 10.01, 0.5)
def dist(xy): return TREE.query(xy)[0]


def one(base):
    Z = np.load(base + "_defects.npz"); t, n, rows = Z["t"], Z["n"], Z["rows"]; off = np.r_[0, np.cumsum(n)]
    fr = [rows[off[i]:off[i + 1]] for i in range(len(t)) if t[i] >= TSKIP]; nf = len(fr)
    T = [tr for tr in track_new(fr, 9.0) if tr["sg"] > 0]; st = es = es40 = 0; xm = []
    for tr in T:
        nn = np.array(tr["n"]); d = dist(np.array(tr["xy"])); pos = dict(zip(nn, d))
        for m, dv in zip(nn, d):
            if m % 5 or dv >= 4 or m + 20 >= nf: continue
            st += 1; d10 = pos.get(m + 10); d20 = pos.get(m + 20)
            es += d10 is not None and d10 > 6; es40 += d20 is not None and d20 > 6
            hist = [pos.get(j) for j in range(m - 10, m + 1)]
            if m >= 10 and not any(h_ is None or h_ >= 4 for h_ in hist):
                fut = [pos[j] for j in range(m + 1, m + 11) if j in pos]
                if fut: xm.append(max(fut))
    nw, qb = [], []
    for f in fr:
        if not len(f): nw.append(0); qb.append(0); continue
        d = dist(f[:, :2]); nw.append(((d < 4) & (f[:, 2] > 0)).sum()); qb.append(f[d > 6, 2].sum())
    ob = [o for o in json.load(open(base + ".json"))["obs"] if o["t"] >= TSKIP]; xm = np.array(xm)
    name = os.path.basename(base); cond, k, s = name.split("_")
    return dict(run=name, cond=cond, k=int(k[1:]), seed=int(s[1:]), starts=st, p_esc=es / max(st, 1), p_esc40=es40 / max(st, 1), nres=len(xm),
                nwall=np.mean(nw), qb=np.mean(qb), coh=np.nanmean([o.get("coh", np.nan) for o in ob]), circ=np.mean([abs(o["circ"]) for o in ob]),
                vrms=np.mean([o["vrms"] for o in ob]), **{f"S{x:g}": (xm > x).mean() if len(xm) else np.nan for x in XG})


if __name__ == "__main__":
    bases = sorted(f[:-5] for f in glob.glob(ROOT + "/*_k*_s*.json"))
    R = [one(b) for b in bases]
    with open(os.path.join(H, "ktest_escape.csv"), "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=list(R[0])); w.writeheader(); w.writerows(R)
    L = ["CONTINUUM K TEST (aspect-2 ellipse, zeta 1, lambda 2; t >= 500; mean over 2 seeds, range in brackets)",
         "cond | k | wall +1/2 | Q_b | coherence | |circ| | v_rms | p_esc (20 t.u.) | p_esc (40 t.u.) | S(4) | S(6) | starts"]
    for cond in ("ref", "K2", "K2G05", "N05"):
        for k in (5, 14):
            g = [r for r in R if r["cond"] == cond and r["k"] == k]
            if not g: continue
            f = lambda c: f"{np.mean([r[c] for r in g]):.4g} [{min(r[c] for r in g):.3g}-{max(r[c] for r in g):.3g}]"
            L.append(f"{cond:5s} | {k:2d} | {f('nwall')} | {f('qb')} | {f('coh')} | {f('circ')} | {f('vrms')} | {f('p_esc')} | {f('p_esc40')} | {f('S4')} | {f('S6')} | {sum(r['starts'] for r in g)}")
    for cond in ("ref", "K2", "K2G05", "N05"):
        a = [r["p_esc"] for r in R if r["cond"] == cond and r["k"] == 5]; b = [r["p_esc"] for r in R if r["cond"] == cond and r["k"] == 14]
        if a and b: L.append(f"crowding suppression {cond}: p_esc(k5) / p_esc(k14) = {np.mean(a) / max(np.mean(b), 1e-9):.2f}")
    L.append("MPCD reference (aspect 2, R_eq 1.5 l_a, U 10): p_esc 0.140 (k5) vs 0.0132 (k14) -> 10.6x suppression")
    open(os.path.join(H, "ktest_escape.txt"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
