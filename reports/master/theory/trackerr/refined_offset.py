#!/usr/bin/env python3
"""Is the +4-5 deg offset between position shift and tilt a tracker artefact? (user, 2026-10-22). Real frames of the active main sweep
(server_main/main_act0.07, seeds 1-3, all k and densities; director frames of the second half). For every wall +1/2 found by the solver's tracker
(integer-cell centre of mass of the flagged core cells): (i) refined sub-cell position = centre of mass of (S_ref - S) over cells within 2.0 cells
of the tracker position (S cell order parameter, S_ref = 90th percentile of S in the frame; weights clipped at 0); (ii) orientation from div Q
averaged within 2 cells of the refined position. Per run (frames with exactly 2k-2 wall +1/2 and coherent tracker tilt, |<e^{i tau}>| > 0.9):
offset D = arg Z_N - tau with tracker positions + tracker tilt, refined positions + tracker tilt, refined positions + div-Q tilt; signed by the tilt.
Also the mean tangential displacement refined - tracker (cells, signed by the tilt) and its component along the comet axis.
-> refined_offset.txt"""
import os, sys, glob, json, warnings, numpy as np
from multiprocessing import Pool
warnings.filterwarnings("ignore")
ROOT = "/home/baloewe/Physics/Overcharged"; sys.path.insert(0, ROOT + "/reports/master/theory/tiltcrowd"); sys.path.insert(0, ROOT + "/runs/cavity")
import bulk_spiral as BS  # noqa
from edge_criteria import read_with_angle  # noqa
def one(run):
    try:
        bc = json.load(open(run + "/input.json"))["BC"][0]; R = bc["R"]; k = bc["customAnch"]; Q = np.array(bc["Q"][:2], float); N = int(2 * k - 2)
        D = BS.frames(run + "/directorfield.dat.gz", 8); td = np.unique(D[:, 0]); td = td[td >= td[0] + 0.5 * (td[-1] - td[0])]
        t, fr = read_with_angle(run, Q, 0.0); t = np.array(t); n = int(D[:, 1].max()) + 1; m = int(D[:, 2].max()) + 1
        X, Y = np.meshgrid(np.arange(n) + 0.5 - Q[0], np.arange(m) + 0.5 - Q[1], indexing="ij")
        A0, A1, T0, T1, DT, DA = [], [], [], [], [], []
        for tt in td:
            i = int(np.argmin(np.abs(t - tt)))
            if abs(t[i] - tt) > 1e-6: continue
            f = fr[i]; r = np.hypot(f[:, 0], f[:, 1]); w = (f[:, 2] > 0) & (R - r < 4)
            if w.sum() != N: continue
            tau_t = f[w, 3] - np.arctan2(f[w, 1], f[w, 0]); zt = np.mean(np.exp(1j * tau_t))
            if abs(zt) < 0.9: continue
            d = D[D[:, 0] == tt]; S = np.zeros((n, m)); q1 = np.zeros((n, m)); q2 = np.zeros((n, m)); th = np.arctan2(d[:, 5], d[:, 4]); ix, iy = d[:, 1].astype(int), d[:, 2].astype(int)
            S[ix, iy] = d[:, 7]; q1[ix, iy] = d[:, 7] * np.cos(2 * th); q2[ix, iy] = d[:, 7] * np.sin(2 * th)
            dvx = np.gradient(q1, axis=0) + np.gradient(q2, axis=1); dvy = np.gradient(q2, axis=0) - np.gradient(q1, axis=1)
            Sref = np.percentile(S[np.hypot(X, Y) < R - 1], 90); xr, yr, tr = [], [], []
            for x0, y0 in zip(f[w, 0], f[w, 1]):
                mk = (np.hypot(X - x0, Y - y0) < 2.0) & (np.hypot(X, Y) < R); wt = np.clip(Sref - S[mk], 0, None)
                if wt.sum() <= 0: xr.append(x0); yr.append(y0)
                else: xr.append(np.sum(wt * X[mk]) / wt.sum()); yr.append(np.sum(wt * Y[mk]) / wt.sum())
                mk2 = np.hypot(X - xr[-1], Y - yr[-1]) < 2.0; tr.append(np.arctan2(dvy[mk2].mean(), dvx[mk2].mean()))
            xr, yr, tr = map(np.array, (xr, yr, tr)); ph0 = np.arctan2(f[w, 1], f[w, 0]); ph1 = np.arctan2(yr, xr)
            sg = np.sign(np.angle(zt)) or 1.0
            A0.append(np.exp(1j * N * ph0).mean()); A1.append(np.exp(1j * N * ph1).mean()); T0.append(zt); T1.append(np.mean(np.exp(1j * (tr - ph1))))
            dphi = np.angle(np.exp(1j * (ph1 - ph0))) * np.hypot(xr, yr); DT.append(sg * np.mean(dphi))
            ax_ = f[w, 3]; DA.append(np.mean((xr - f[w, 0]) * np.cos(ax_) + (yr - f[w, 1]) * np.sin(ax_)))
        if len(A0) < 5: return None
        c = lambda z: np.angle(np.mean(np.array(z) / np.abs(np.array(z))))
        sg = np.sign(c(T0)) or 1.0
        out = dict(k=k, a=2 * np.pi * (R - 1.5) / N, n=len(A0), D_tt=sg * np.angle(np.exp(1j * (c(A0) - c(T0)))), D_rt=sg * np.angle(np.exp(1j * (c(A1) - c(T0)))),
                   D_rd=sg * np.angle(np.exp(1j * (c(A1) - c(T1)))), dtan=np.mean(DT), dax=np.mean(DA), tau=sg * c(T0))
        return out
    except Exception:
        return None
if __name__ == "__main__":
    runs = sorted(os.path.dirname(p) for p in glob.glob(ROOT + "/runs/cavity/server_main/main_act0.07/rho*/k*/s0[1-3]/directorfield.dat.gz"))
    with Pool(8) as p: res = [x for x in p.map(one, runs, chunksize=1) if x]
    L = [f"{len(res)} runs. Offset D = arg Z_N - tau (deg, signed by the tilt), median per spacing bin: tracker pos + tracker tilt | refined pos + tracker tilt | refined pos + div-Q tilt | mean tangential shift refined - tracker (cells, signed) | shift along comet axis (cells)"]
    a = np.array([x["a"] for x in res])
    for lo, hi in ((0, 8), (8, 12), (12, 20), (20, 80)):
        mm = (a >= lo) & (a < hi)
        if not mm.any(): continue
        g = lambda key: np.median([x[key] for x, m_ in zip(res, mm) if m_])
        L.append(f"spacing [{lo},{hi}): {mm.sum():3d} runs | {np.degrees(g('D_tt')):+5.1f} | {np.degrees(g('D_rt')):+5.1f} | {np.degrees(g('D_rd')):+5.1f} | {g('dtan'):+.3f} | {g('dax'):+.3f}   (tilt {np.degrees(g('tau')):.1f})")
    open("refined_offset.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
