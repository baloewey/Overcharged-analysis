#!/usr/bin/env python3
"""Which factors explain the higher escape rate of wall +1/2 in elongated cavities? Shape sweep (675 runs), second half, frames every 10 t.u.
Per wall +1/2 sample: escape within 20 t.u. (as selfprop.py) and covariates: depth; cos(tau) (inward head component); distance to the nearest OTHER bulk +1/2 (depth > 4, capped 15)
and its heading cosine towards the wall defect; number of bulk +1/2 within 8 cells; region (tip/middle/flank); tilt coherence of the wall layer in that frame; R_eq.
Pooled logistic regression (scipy L-BFGS, binned covariates as one-hot) with a class indicator ELONG (aspect >= 1.25). The log-odds gap of ELONG with no covariates vs with each
group of covariates added (and all together) = how much of the excess each group explains. -> decomp.txt (features cached in decomp_features.npz)"""
import csv, glob, json, os, sys
import numpy as np
from multiprocessing import Pool
from scipy.optimize import minimize
from scipy.spatial import cKDTree
ROOT = "/home/baloewe/Physics/Overcharged/runs"; sys.path.insert(0, ROOT + "/cavity")
from edge_criteria import read_with_angle  # noqa
from wall_fourier import wall_polyline  # noqa
SW = ROOT + "/shapes_server"; H = os.path.dirname(os.path.abspath(__file__))
def one(rel):
    try:
        run = os.path.join(SW, rel); bc = json.load(open(run + "/input.json"))["BC"][0]; Q = bc["Q"]
        W, s, P = wall_polyline(bc); tree = cKDTree(W); A = np.abs(W[:, 0]).max(); Wc = W[:, 0] + 1j * W[:, 1]
        Tn = np.roll(Wc, -1) - np.roll(Wc, 1); Nn = -1j * Tn / np.abs(Tn)
        if np.mean((Nn * np.conj(Wc)).real) < 0: Nn = -Nn
        t, fr = read_with_angle(run, np.array(Q[:2], float), 0.5); fr = fr[::5]; rows = []
        PL = [f[f[:, 2] > 0] if len(f) else np.zeros((0, 4)) for f in fr]
        for a in range(len(PL) - 2):
            p = PL[a]
            if not len(p): continue
            d, j = tree.query(p[:, :2]); w = d < 4
            if not w.any(): continue
            tau = np.angle(np.exp(1j * (p[:, 3] - np.angle(Nn[j])))); coh = abs(np.mean(np.exp(1j * tau[w]))); bulk = np.flatnonzero(d > 4)
            k = json.load(open(run + "/input.json"))["BC"][0]["customAnch"]; nvac = (2 * k - 2) - int(w.sum()); sw = np.sort(s[j[w]]); spc = P / (2 * k - 2)
            pos = p[w, :2]; cur = pos.copy(); alive = np.ones(len(pos), bool); dfin = np.zeros(len(pos))
            for b in (1, 2):
                if not len(PL[a + b]): alive[:] = False; break
                q, k = cKDTree(PL[a + b][:, :2]).query(cur); alive &= q < 6; cur = PL[a + b][k, :2]; dfin = tree.query(cur)[0]
            esc = alive & (dfin > 6)
            for ii, i in enumerate(np.flatnonzero(w)):
                if len(bulk):
                    dv = p[i, :2] - p[bulk, :2]; r = np.hypot(*dv.T); k = np.argmin(r); nb = int((r < 8).sum())
                    head = np.array([-np.cos(p[bulk[k], 3]), -np.sin(p[bulk[k], 3])]); hc = float(np.dot(head, dv[k]) / max(r[k], 1e-9)); rn = min(r[k], 15.0)
                else: rn, hc, nb = 15.0, 0.0, 0
                reg = 0 if abs(p[i, 0]) / A > 0.5 else (2 if abs(p[i, 0]) / A < 0.25 else 1)
                si = s[j[i]]; dd_ = np.sort(np.abs(((sw - si + P / 2) % P) - P / 2)); gaps = dd_[1:3] / spc if len(dd_) >= 3 else np.array([2.0, 2.0])
                rows.append((d[i], np.cos(tau[i]), rn, hc, nb, reg, coh, float(esc[ii]), nvac, gaps.max(), gaps.min()))
        return rel, np.array(rows)
    except Exception:
        return rel, None
def fit(X, y):
    def f(b):
        z = X @ b; p = 1 / (1 + np.exp(-z)); eps = 1e-12
        return -np.sum(y * np.log(p + eps) + (1 - y) * np.log(1 - p + eps)), X.T @ (p - y)
    return minimize(f, np.zeros(X.shape[1]), jac=True, method="L-BFGS-B", options=dict(maxiter=500)).x
def onehot(v, edges):
    k = np.digitize(v, edges); return np.eye(len(edges) + 1)[k][:, 1:]
if __name__ == "__main__":
    man = {r["dir"]: r for r in csv.DictReader(open(SW + "/manifest_shape_sweep.tsv"), delimiter="\t")}
    cache = os.path.join(H, "decomp_features2.npz")
    if os.path.exists(cache): Z = dict(np.load(cache)); F, EL, RQ = Z["F"], Z["EL"], Z["RQ"]
    else:
        runs = sorted(os.path.relpath(os.path.dirname(f), SW) for f in glob.glob(SW + "/*/*/*/*/s*/DONE"))
        with Pool(8) as pool: res = pool.map(one, runs, chunksize=4)
        F, EL, RQ = [], [], []
        for rel, a in res:
            if a is None or not len(a): continue
            F.append(a); EL.append(np.full(len(a), float(man[rel]["aspect"]) >= 1.25)); RQ.append(np.full(len(a), float(man[rel]["Req_la"])))
        F, EL, RQ = np.concatenate(F), np.concatenate(EL).astype(float), np.concatenate(RQ); np.savez_compressed(cache, F=F, EL=EL, RQ=RQ)
    y = F[:, 7]; n = len(y)
    groups = {"depth": onehot(F[:, 0], [1, 1.5, 2, 2.5, 3]), "orientation": onehot(F[:, 1], [0, 0.25, 0.5, 0.75, 0.9]),
              "nearest bulk +1/2 (distance x heading)": np.column_stack([onehot(F[:, 2], [3, 5, 8, 12]), onehot(F[:, 2], [3, 5, 8, 12]) * F[:, 3:4]]),
              "bulk +1/2 count (r < 8)": onehot(F[:, 4], [1, 2, 3, 5]), "region": onehot(F[:, 5], [0.5, 1.5]), "layer coherence": onehot(F[:, 6], [0.3, 0.5, 0.7, 0.85]),
              "size": onehot(RQ, [1.25, 1.75]), "vacancies in the cavity": onehot(F[:, 8], [0.5, 1.5, 2.5, 3.5, 5]),
              "gaps to neighbouring wall defects": np.column_stack([onehot(F[:, 9], [1.2, 1.5, 2.0]), onehot(F[:, 10], [0.8, 1.0])])}
    base = np.column_stack([np.ones(n), EL])
    L = []; say = lambda x: (L.append(x), print(x, flush=True))
    b0 = fit(base, y); g0 = b0[1]; say(f"{n} samples; P(escape) circle-like {y[EL == 0].mean():.4f}, elongated {y[EL == 1].mean():.4f}; class log-odds gap, no covariates: {g0:.3f} (odds ratio {np.exp(g0):.2f})")
    for name, G in groups.items():
        b = fit(np.column_stack([base, G]), y); say(f"  + {name:40s}: gap {b[1]:.3f} -> explains {(g0 - b[1]) / g0:+.0%}")
    allG = np.column_stack([base] + list(groups.values())); b = fit(allG, y); say(f"  + ALL covariates together                 : gap {b[1]:.3f} -> explains {(g0 - b[1]) / g0:+.0%}")
    # effect sizes within the full model for the bulk +1/2 terms
    say("full model, nearest bulk +1/2 coefficients (distance bins <3, 3-5, 5-8, 8-12, >12; reference = first absent): see decomp_coef.txt")
    np.savetxt(os.path.join(H, "decomp2_coef.txt"), b)
    open(os.path.join(H, "decomp2.txt"), "w").write("\n".join(L) + "\n")
