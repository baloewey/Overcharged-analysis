#!/usr/bin/env python3
"""Lock slope from FLUCTUATIONS (user, 2026-11-27). The scan in lock_slope_check.py cannot separate s = 2k-2 from s = 0 on an evenly spaced ring
(aliasing: (2k-2) phi_j is the same for all sites). Here, in frames with exactly N = 2k-2 wall +1/2: fit the regular ring phi_j^0 = alpha + 2 pi j / N
(alpha from the circular mean of N phi_j), deviations dphi_j = phi_j - phi_j^0 and dtau_j = tau_j - <tau>_frame (circular). Lock: dtau = (2k-2) dphi;
no lock (tilt independent of position): dtau = 0. OLS slope of dtau on dphi pooled per k; also the slope with the per-defect tracker noise in
dphi estimated (attenuation factor = var(true dphi) / var(dphi); lag-1 covariance at 2 t.u. used for the true variance). -> lock_slope_fluct.txt"""
import glob, json, os, sys, warnings, numpy as np
from multiprocessing import Pool
warnings.filterwarnings("ignore")
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
import static_band as SB  # noqa


def run_one(run):
    try:
        bc = json.load(open(run + "/input.json"))["BC"][0]; R = bc["R"]; k = bc["customAnch"]; Q = np.array(bc["Q"][:2], float); N = 2 * k - 2
        t, fr = SB.read_with_angle(run, Q, 0.5); X, Yv = [], []
        for f in fr:
            if not len(f): continue
            zs = f[:, 0] + 1j * f[:, 1]; w = (f[:, 2] > 0) & (R - np.abs(zs) < 4)
            if w.sum() != N: continue
            phi = np.sort(np.angle(zs[w])); o = np.argsort(np.angle(zs[w])); tau_all = np.angle(np.exp(1j * (f[w, 3] - np.angle(zs[w]))))[o]
            alpha = np.angle(np.mean(np.exp(1j * N * phi))) / N; j = np.round((phi - alpha) / (2 * np.pi / N)); dphi = np.angle(np.exp(1j * (phi - alpha - 2 * np.pi * j / N)))
            dtau = np.angle(np.exp(1j * (tau_all - np.angle(np.mean(np.exp(1j * tau_all))))))
            X += list(dphi); Yv += list(dtau)
        if len(X) < 100: return None
        return dict(k=k, x=np.array(X), y=np.array(Yv))
    except Exception as ex:
        return dict(err=repr(ex)[:150], run=run)


if __name__ == "__main__":
    runs = [os.path.dirname(p) for p in sorted(glob.glob(SB.ROOT + "/runs/cavity/server_main/main_act0.07/rho*/k*/s[01][0-9]/defects.dat.gz"))
            if abs(json.load(open(os.path.dirname(p) + "/input.json"))["BC"][0]["R"] - 31.83) < 0.5]
    with Pool(int(os.environ.get("NPROC", 12))) as p: res = p.map(run_one, runs, chunksize=1)
    ok = [r for r in res if r and "err" not in r]
    L = ["LOCK SLOPE FROM FLUCTUATIONS (active MPCD circles R 31.83; frames with exactly 2k-2 wall +1/2): dtau_j vs dphi_j about the frame's regular ring",
         "k | defect-frames | sd dphi (deg) = sd displacement (cells) | sd dtau (deg) | OLS slope dtau/dphi | lock 2k-2 | slope / (2k-2) | corr"]
    for k in sorted({r["k"] for r in ok}):
        x = np.concatenate([r["x"] for r in ok if r["k"] == k]); y = np.concatenate([r["y"] for r in ok if r["k"] == k])
        sl = np.cov(x, y)[0, 1] / np.var(x); cr = np.corrcoef(x, y)[0, 1]
        L.append(f"{k:2d} | {len(x):7d} | {np.degrees(x.std()):5.2f} = {x.std() * 31.83:4.2f} | {np.degrees(y.std()):5.1f} | {sl:6.2f} | {2 * k - 2:2d} | {sl / (2 * k - 2):.2f} | {cr:.2f}")
    open(os.path.join(H, "lock_slope_fluct.txt"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
