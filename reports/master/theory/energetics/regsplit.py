#!/usr/bin/env python3
"""Split the 36-run escape signatures by wall region (tip / mid / flank) of the escape site. Event q of each run = energetics.events(run, max_events=50)[q] (same list as
forcefit/vortdef/qcrit). Joins forcefit_ff36.npy (run, q), vortdef_ff36.npy + _runs.npy and qcrit_ff36s25.npy + _runs.npy. -> regsplit_ff36.txt"""
import os, sys, numpy as np
from multiprocessing import Pool
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
import energetics as En  # noqa
from tipcount import reg  # noqa
from forcefit import fit  # noqa
def regions(run):
    J, bc, k, t, fr, ev = En.events(run, max_events=50); B = float(bc["R"]); A = B * float(bc.get("aInv", [1])[0] or 1)
    return run, [reg(e["chain"][e["n"]], A, B) for e in ev]
if __name__ == "__main__":
    vr = np.load("vortdef_ff36_runs.npy", allow_pickle=True); runs = sorted(set(vr[:, 0]))
    with Pool(6) as pool: RG = dict(pool.map(regions, runs))
    L = []; RN = ("tip", "mid", "flank")
    # force fit by region
    a = np.load("forcefit_ff36.npy", allow_pickle=True); key = {os.path.basename(os.path.dirname(r)) + "/" + os.path.basename(r): r for r in runs}
    rg = np.array([RG.get(key.get(x[0], ""), [None] * 999)[int(x[2])] if x[1] in ("esc", "ctl") and key.get(x[0]) in RG and int(x[2]) < len(RG[key[x[0]]]) else -1 for x in a], dtype=object)
    R = dict(kind=a[:, 1].astype(str), ev=np.array([f"{x}|{y}" for x, y in zip(a[:, 0], a[:, 2])]), lag=a[:, 3].astype(int), d=a[:, 4].astype(float), ct=a[:, 5].astype(float),
             uc=a[:, 6].astype(float), ur=a[:, 7].astype(float), F0=a[:, 8].astype(float), v=a[:, 9].astype(float)); R["jet"] = R["uc"] - R["ur"]
    L.append("FORCE FIT by region (escapers): n events | alpha jet | beta ring | mu_n | const | jet at lag -10 / 0")
    for g in range(3):
        m = (R["kind"] == "esc") & (rg == g); S = {k_: v_[m] for k_, v_ in R.items()}
        if m.sum() < 100: continue
        b, e, r2, n = fit(S, ["jet", "ur", "F0"], nb=200); ne = len(set(S["ev"]))
        j10 = np.nanmean(S["jet"][S["lag"] == -5]); j0 = np.nanmean(S["jet"][S["lag"] == 0])
        L.append(f"  {RN[g]:5s} | {ne:3d} | {b[0]:+.2f}±{e[0]:.2f} | {b[1]:+.2f}±{e[1]:.2f} | {b[2]:+.4f}±{e[2]:.4f} | {b[3]:+.3f} | {j10:+.3f} / {j0:+.3f}")
    # vortdef / qcrit by region
    for nm, fn, cols in (("VORTICITY (escapers): wall-side lobe / free at lag", "vortdef_ff36", (12,)), ("Q CRITERION (escapers): Qn and P(d_Q<1) at lag", "qcrit_ff36s25", (5, 7))):
        A_ = np.load(fn + ".npy", allow_pickle=True); RR = np.load(fn + "_runs.npy", allow_pickle=True)
        kind = A_[:, 0].astype(str); q = A_[:, 1].astype(int); lag = A_[:, 2].astype(int)
        g = np.array([RG[RR[i, 0]][q[i]] if kind[i] in ("esc", "ctl") and q[i] >= 0 and q[i] < len(RG.get(RR[i, 0], [])) else -1 for i in range(len(A_))])
        L.append(f"\n{nm} -20 / -10 / 0 / +4 t.u.")
        for gg in range(3):
            row = []
            for lg in (-10, -5, 0, 2):
                m = (kind == "esc") & (g == gg) & (lag == lg)
                if fn.startswith("vortdef"):
                    wl, wr, sgm, lw, lat = (A_[m, i].astype(float) for i in (8, 9, 11, 12, 13)); O0 = 0.118
                    wl, wr = wl * sgm, wr * sgm; wallL = np.where(lw > 0, wl / O0, -wr / O0); wallL[lat < 0.5] = np.nan; row.append(f"{np.nanmean(wallL):+.2f}")
                else: row.append(f"{np.nanmean(A_[m, 5].astype(float)):+.2f},{np.mean(A_[m, 7].astype(float) < 1):.2f}")
            L.append(f"  {RN[gg]:5s} (n {int(((kind == 'esc') & (g == gg) & (lag == 0)).sum())}): " + "  ".join(row))
    print("\n".join(L)); open("regsplit_ff36.txt", "w").write("\n".join(L) + "\n")
