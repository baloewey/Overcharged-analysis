#!/usr/bin/env python3
"""SERVER COPY (paths -> ~/Overcharged/runs/cavity/main_act0.07; frame reader inlined) of curvature_profiles.py.
Tests of links A and B of curvature_theory.md on the main active circle sweep (lambda 2, R/l_a 0.8-2.2, 10 seeds; director and 100-t.u. flow
averages every 100 t.u., second half). Per run, azimuthally averaged profiles vs depth y = R - r (0.5-cell bins), circulation-signed:
  v(y)      azimuthal velocity;  Lambda(y) = |v_r + v/r| / |v_r - v/r| (rotation/strain, signed form); first order in 1/R: Lambda = 1 - 2 l_v/R with the flow length l_v = v/|v'|
            (= y for a linear profile); link A is an identity for azimuthal flow: its content is the profile shape l_v(y)
  alpha(y)  director angle from the local tangent over 1-cell windows [y, y+1), arg<S exp(2i(theta - phi - pi/2))>/2 (band_profile.py), uniformity |Z|/<S>
  alpha_J   = 0.5 arccos(Lambda/lambda) (link B: alpha -> alpha_J beyond l_E), from the run's own Lambda(y)
Fits: Lambda(y) - 1 vs y/R over y 1.5-4 cells, pooled over runs (prediction -2); alpha - alpha_J vs y.
-> curvature_profiles.txt, curvature_profiles.npz"""
import os, sys, glob, json, numpy as np
from multiprocessing import Pool
H = os.path.dirname(os.path.abspath(__file__)); import gzip
class BS:                                          # server copy: frame reader of tiltcrowd/bulk_spiral.py inlined
    @staticmethod
    def frames(fn, ncol): return np.loadtxt((l for l in gzip.open(fn, "rt") if len(l.split()) == ncol and l.split()[0][0].isdigit()))
LA, LAM = 22.21, 2.0; YE = np.arange(0, 12.01, 0.5); YC = 0.5 * (YE[1:] + YE[:-1])


def one(run):
    try:
        bc = json.load(open(run + "/input.json"))["BC"][0]; R = bc["R"]; k = bc["customAnch"]; Q = np.array(bc["Q"][:2], float)
        if not 0.8 <= R / LA <= 2.2: return None
        D = BS.frames(run + "/directorfield.dat.gz", 8); F = BS.frames(run + "/flowfield.dat.gz", 7)
        tf = np.unique(F[:, 0]); tf = tf[tf >= tf[0] + 0.5 * (tf[-1] - tf[0])]
        V = np.zeros(len(YC)); NV = np.zeros(len(YC)); Z = np.zeros(len(YC), complex); SM = np.zeros(len(YC)); ND = np.zeros(len(YC))
        for tt in tf:
            f = F[F[:, 0] == tt]; x = f[:, 1] + 0.5 - Q[0]; y = f[:, 2] + 0.5 - Q[1]; r = np.hypot(x, y); vt = (x * f[:, 5] - y * f[:, 4]) / np.maximum(r, 1e-9)
            circ = np.sign(vt[(r < R - 1) & (r > 0.3 * R)].mean())
            if circ == 0: continue
            dep = R - r; ok = (dep >= 0) & (dep < YE[-1]); b = np.digitize(dep[ok], YE) - 1
            np.add.at(V, b, circ * vt[ok]); np.add.at(NV, b, 1)
            d = D[np.abs(D[:, 0] - tt) < 1e-6]
            if not len(d): continue
            xd = d[:, 1] + 0.5 - Q[0]; yd = d[:, 2] + 0.5 - Q[1]; rd = np.hypot(xd, yd); ph = np.arctan2(yd, xd); dd = R - rd
            m = (dd >= 0) & (dd < YE[-1]) & (d[:, 7] > 0); bd = np.digitize(dd[m], YE) - 1
            np.add.at(Z, bd, d[m, 7] * np.exp(2j * circ * (np.arctan2(d[m, 5], d[m, 4]) - ph[m] - np.pi / 2))); np.add.at(SM, bd, d[m, 7]); np.add.at(ND, bd, 1)
        v = V / np.maximum(NV, 1); r = R - YC; dv = -np.gradient(v, YC)            # dv/dr = -dv/dy
        Lam = np.abs(dv + v / r) / np.maximum(np.abs(dv - v / r), 1e-12)         # signed form: valid past the velocity maximum too
        lv = v / np.maximum(np.abs(dv), 1e-12)                                     # flow length v/|v'|: first order Lambda = 1 - 2 lv/R
        Z2 = Z + np.r_[Z[1:], 0]; S2 = SM + np.r_[SM[1:], 0]                       # 1-cell windows [y, y+1): half-cell rings alias with the defect rows
        alpha = np.degrees(0.5 * np.angle(Z2)); unif = np.abs(Z2) / np.maximum(S2, 1e-9)
        aJ = np.where(Lam / LAM <= 1, np.degrees(0.5 * np.arccos(np.clip(Lam / LAM, -1, 1))), np.nan)   # Lambda > lambda: tumbling, no aligned angle
        rho = 2 * abs(k) / (2 * np.pi * R)
        return dict(run=run.split("main_act0.07/")[-1], k=k, R=R, rho=rho, v=v, Lam=Lam, lv=lv, alpha=alpha, unif=unif, aJ=aJ)
    except Exception as ex:
        return dict(run=run, err=str(ex))


if __name__ == "__main__":
    runs = sorted(os.path.dirname(p) for p in glob.glob(os.path.expanduser("~/Overcharged/runs/cavity/main_act0.07/rho*/k*/s[01][0-9]/directorfield.dat.gz")))
    with Pool(int(os.environ.get("NPROC", 12))) as p: res = p.map(one, runs, chunksize=2)
    err = [r for r in res if r and "err" in r]; res = [r for r in res if r and "err" not in r]
    np.savez_compressed(os.path.join(H, "curvature_profiles.npz"), YC=YC, runs=np.array([r["run"] for r in res]), k=np.array([r["k"] for r in res]),
                        R=np.array([r["R"] for r in res]), rho=np.array([r["rho"] for r in res]), v=np.array([r["v"] for r in res]),
                        Lam=np.array([r["Lam"] for r in res]), alpha=np.array([r["alpha"] for r in res]), unif=np.array([r["unif"] for r in res]),
                        aJ=np.array([r["aJ"] for r in res]), lv=np.array([r["lv"] for r in res]))
    L = [f"CURVATURE PROFILES (links A, B of curvature_theory.md): {len(res)} runs, {len(err)} errors {[e['run'][-30:] for e in err[:3]]}"]
    sel = (YC >= 1.5) & (YC <= 4.0)
    x = np.concatenate([YC[sel] / r["R"] for r in res]); yv = np.concatenate([r["Lam"][sel] - 1 for r in res]); ok = np.isfinite(yv)
    A = np.c_[np.ones(ok.sum()), x[ok]]; c, *_ = np.linalg.lstsq(A, yv[ok], rcond=None)
    L.append(f"LINK A pooled (y 1.5-4 cells): Lambda - 1 = {c[0]:+.3f} {c[1]:+.2f} y/R   (prediction 0 - 2 y/R)")
    for lo, hi in ((0.03, 0.055), (0.055, 0.085), (0.085, 0.16)):
        g = [r for r in res if lo <= r["rho"] < hi]
        if len(g) < 5: continue
        xg = np.concatenate([YC[sel] / r["R"] for r in g]); yg = np.concatenate([r["Lam"][sel] - 1 for r in g]); okg = np.isfinite(yg)
        cg, *_ = np.linalg.lstsq(np.c_[np.ones(okg.sum()), xg[okg]], yg[okg], rcond=None)
        L.append(f"   rho {lo:.3f}-{hi:.3f} ({len(g)} runs): Lambda - 1 = {cg[0]:+.3f} {cg[1]:+.2f} y/R")
    L.append("Profiles by size (sparse rho < 0.055 | crowded rho >= 0.085), mean over runs: y | Lambda | 1-2y/R(mean R) | flow length v/|v'| | alpha over [y, y+1) (unif) | alpha_J | alpha - alpha_J")
    for lab, cond in (("sparse", lambda r: r["rho"] < 0.055), ("crowded", lambda r: r["rho"] >= 0.085)):
        for Rlo, Rhi in ((0.8, 1.2), (1.2, 1.7), (1.7, 2.2)):
            g = [r for r in res if cond(r) and Rlo <= r["R"] / LA < Rhi]
            if len(g) < 3: continue
            Rm = np.mean([r["R"] for r in g]); L.append(f" {lab} R/l_a {Rlo}-{Rhi} ({len(g)} runs, mean R {Rm:.1f}):")
            for yy in (1.25, 2.25, 3.25, 4.25, 6.25, 8.25):
                i = int(np.argmin(np.abs(YC - yy))); f = lambda key: np.nanmean([r[key][i] for r in g])
                L.append(f"    {YC[i]:5.2f} | {f('Lam'):6.3f} | {1 - 2 * YC[i] / Rm:6.3f} | {f('lv'):5.2f} | {f('alpha'):+6.1f} ({f('unif'):.2f}) | {f('aJ'):5.1f} | {f('alpha') - f('aJ'):+6.1f}")
    open(os.path.join(H, "curvature_profiles.txt"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
