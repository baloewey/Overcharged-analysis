#!/usr/bin/env python3
"""Depth (finite-shear) correction of the tilt law as the cause of the crowding drop (user, 2026-10-16). Master Eq. (tilt):
    tau = arccos(1/lambda) + d sqrt(gamma1 s / K) (lambda^2 - 1)^(-1/4)
(d = defect depth below the wall, s = wall shear rate). The bare prefactor is known not to transfer to MPCD (gamma1 ~ 60, K ~ 32.5 give +16 deg; matching
needed gamma1/K 0.07-0.18), so test the FORM with one effective prefactor C: tau - 60 = C d sqrt(s) (lambda 2).
Data, per run: depth d (median R - r of wall +1/2) and tilt (coherent frames) from defect_flow; wall shear s = d<v_theta>/dr just inside the layer,
from the azimuthally averaged profile (wallshear.one: e_band + w_band in the band 1-4 cells from the wall).
Sets: active circles at R 31.83 (k 3-10, 10 seeds, main sweep) and the whole active 10-seed window; passive imposed-vortex runs (swirl 0.01 and 0.03).
-> depth_correction.txt, fig_depth_correction.png"""
import os, sys, glob, json, numpy as np
from multiprocessing import Pool
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H); ROOT = "/home/baloewe/Physics/Overcharged"
import defect_flow as DF, wallshear as WS


def one(run):
    try:
        a = DF.one(run); w = WS.one(run)
        return dict(run=run, a=a["a"], d=a["depth"], tilt=a["tilt"], s=w["e_band"] + w["w_band"], k=a["k"], R=a["R"])
    except Exception as ex:
        return None


if __name__ == "__main__":
    pas = sorted(glob.glob(ROOT + "/runs/swirl_crowd/k0[468]_s*") + glob.glob(ROOT + "/runs/swirl_crowd/A03_k*_s*")); pas = [r for r in pas if "A03_k16" not in r]
    act = sorted(os.path.dirname(p) for p in glob.glob(ROOT + "/runs/cavity/server_main/main_act0.07/rho*/k*/s[01][0-9]/flowfield.dat.gz"))
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 6) as pool: res = [x for x in pool.map(one, pas + act, chunksize=2) if x and np.isfinite(x["tilt"])]
    P = [x for x in res if "swirl_crowd" in x["run"]]; A = [x for x in res if "swirl_crowd" not in x["run"]]
    arr = lambda S, c: np.array([x[c] for x in S])
    L = ["DEPTH CORRECTION TEST: tau - 60 = C d sqrt(s)  (lambda 2; d depth, s wall shear rate; one effective C per set and jointly)"]
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
    for name, S, col in (("passive imposed vortex", P, "#1f5fa8"), ("active circles (10-seed window)", A, "#c2452d")):
        d, s, t, a = arr(S, "d"), arr(S, "s"), arr(S, "tilt"), arr(S, "a"); x = d * np.sqrt(s); C = np.sum(x * (t - 60)) / np.sum(x * x)
        r = t - 60 - C * x
        L.append(f"{name}: runs {len(S)}; depth {d.min():.2f}-{d.max():.2f}; shear s {s.min():.3f}-{s.max():.3f}; C = {C:.1f} deg/(cell s^0.5); rms residual {r.std():.1f} deg (raw sd {t.std():.1f}); "
                 f"corr(tilt, d sqrt s) {np.corrcoef(t, x)[0,1]:+.2f}; corr(tilt, d) {np.corrcoef(t, d)[0,1]:+.2f}; corr(tilt, s) {np.corrcoef(t, s)[0,1]:+.2f}; corr(residual, log a) {np.corrcoef(r, np.log(a))[0,1]:+.2f}")
        for lo, hi in ((5, 8), (8, 12), (12, 20), (20, 60)):
            m = (a >= lo) & (a < hi)
            if m.sum(): L.append(f"   spacing [{lo},{hi}): n {m.sum():3d}  depth {np.median(d[m]):.2f}  s {np.median(s[m]):.3f}  tilt {np.median(t[m]):5.1f}  predicted {60 + C * np.median(x[m]):5.1f}")
        ax[0].scatter(x, t, s=10, color=col, alpha=0.6, label=name); xx = np.linspace(0, x.max() * 1.05, 10); ax[0].plot(xx, 60 + C * xx, color=col, lw=1)
        ax[1].scatter(a, d, s=10, color=col, alpha=0.6, label=name); ax[2].scatter(a, r, s=10, color=col, alpha=0.6, label=name)
    ax[0].axhline(60, color="k", ls="--", lw=0.8); ax[0].set_xlabel(r"$d\sqrt{s}$ (cells $\times$ (1/t.u.)$^{1/2}$)"); ax[0].set_ylabel("tilt of the coherent layer (deg)"); ax[0].set_title("(a) Tilt vs the depth correction", fontsize=10)
    ax[1].set_xscale("log"); ax[1].set_xlabel("wall-defect spacing (cells)"); ax[1].set_ylabel("defect depth d (cells)"); ax[1].set_title("(b) Depth vs spacing", fontsize=10)
    ax[2].set_xscale("log"); ax[2].axhline(0, color="0.7"); ax[2].set_xlabel("wall-defect spacing (cells)"); ax[2].set_ylabel("residual (deg)"); ax[2].set_title("(c) What the correction leaves", fontsize=10)
    for a_ in ax: a_.legend(fontsize=7.5, frameon=False); [a_.spines[q].set_visible(False) for q in ("top", "right")]
    fig.tight_layout(); fig.savefig(os.path.join(H, "fig_depth_correction.png"), dpi=200)
    print("\n".join(L)); open(os.path.join(H, "depth_correction.txt"), "w").write("\n".join(L) + "\n")
