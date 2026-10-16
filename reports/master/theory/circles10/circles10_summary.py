#!/usr/bin/env python3
"""Summary of circles10_dynamics.pkl (10-seed circle sweep, act 0.07, rho <= 0.15, k >= 2; +k and mirror -k). By R/l_a band (+k; mirror separately):
wall exchange (escape rate per wall defect, arrivals, residence), duck-duck-goose (entry -> next exit delay vs null; other-defect fraction), pair lifetimes
(flicker fraction < 10 t.u.; established -1/2 lifetime median; death partner/bulk+/wall+), nucleation (pair creations per bulk area per 1000 t.u.; bulk -1/2
density; within-run correlation of creations with the bulk charge fluctuation). -> circles10_summary.txt"""
import os, pickle, numpy as np
H = os.path.dirname(os.path.abspath(__file__)); LA = 22.21
D = pickle.load(open(os.path.join(H, "circles10_dynamics.pkl"), "rb"))
bands = [(0.3, 0.6), (0.6, 0.8), (0.8, 1.2), (1.2, 1.6), (1.6, 2.0), (2.0, 2.4), (2.4, 3.0)]
L = [f"10-SEED CIRCLE DYNAMICS ({len(D)} runs; act 0.07, rho <= 0.15). Rates per 1000 t.u.",
     "set | R/l_a | runs | escape /wall def | arrivals /cavity | DDG: delay median entry->exit vs null (ratio), other-defect frac | pairs: flicker frac, est. -1/2 life median, death partner/bulk+/wall+ | nucleation per 1000 cells^2 | corr(creations, dQ_b) within runs"]
for neg in (False, True):
    for lo, hi in bands:
        S = [x for x in D if (x["k"] < 0) == neg and lo <= x["R"] / LA < hi]
        if not S: continue
        esc = [x["exch"]["escape"] * 1000 for x in S if "escape" in x["exch"]]; arr = [x["exch"]["arrive"] * 1000 for x in S if "arrive" in x["exch"]]
        dl = [p[0] for x in S if "pairs" in x["ddg"] for p in x["ddg"]["pairs"]]; nl = [v for x in S if "null" in x["ddg"] for v in x["ddg"]["null"]]
        oth = [not p[1] for x in S if "pairs" in x["ddg"] for p in x["ddg"]["pairs"]]
        ev = [e for x in S if "ev" in x["life"] for e in x["life"]["ev"]]
        cl = np.array([e[0] for e in ev]); lt = np.array([e[1] for e in ev]); est = lt >= 10
        fr = [np.mean(cl[est] == c) for c in ("PARTNER", "BULK+", "WALL+")] if est.any() else [np.nan] * 3
        unc = est & (cl != "CENSORED")
        nuc, corr = [], []
        for x in S:
            n = x["nuc"]
            if "ev" not in n or n["bulk_area"] <= 0: continue
            nuc.append(n["ev"].sum() / (len(n["ev"]) * n["dt"]) * 1000 / n["bulk_area"] * 1000)
            q = n["qb"] - np.nanmean(n["qb"]); m = np.isfinite(q)
            if m.sum() > 50 and np.std(n["ev"][m]) > 0 and np.std(q[m]) > 0:
                k_ = np.ones(25) / 25; corr.append(np.corrcoef(np.convolve(n["ev"][m], k_, "same"), np.convolve(q[m], k_, "same"))[0, 1])
        L.append(f"{'mirror' if neg else '+k':6s} | {lo:.1f}-{hi:.1f} | {len(S):3d} | {np.mean(esc):.3f} | {np.mean(arr):.2f} | {np.median(dl) if dl else np.nan:5.0f} vs {np.median(nl) if nl else np.nan:5.0f} "
                 f"({np.median(nl)/np.median(dl) if dl else np.nan:.1f}x), {np.mean(oth) if oth else np.nan:.2f} | {np.mean(~est) if len(ev) else np.nan:.2f}, {np.median(lt[unc]) if unc.any() else np.nan:.0f} t.u., "
                 f"{fr[0]:.2f}/{fr[1]:.2f}/{fr[2]:.2f} | {np.mean(nuc) if nuc else np.nan:.2f} | {np.median(corr) if corr else np.nan:+.2f}")
print("\n".join(L)); open(os.path.join(H, "circles10_summary.txt"), "w").write("\n".join(L) + "\n")
