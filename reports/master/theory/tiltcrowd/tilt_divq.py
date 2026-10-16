#!/usr/bin/env python3
"""Tilt measured two ways on the SAME frames (user: is the tilt measured the same way everywhere?). MPCD circles at R 31.83, k 3-20, second half.
At every director frame (every 100 t.u.): wall +1/2 positions from the tracker (defects.dat, nearest frame in time); tilt (a) from the tracker angle,
(b) from the direction of div Q of the cell director field (S-weighted Q, disc of radius 2 around the tracker position) — the definition used by the
continuum codes. Layer tilt per frame = |arg <e^{i tau}>|, coherent frames (C > 0.9 by the tracker); report medians per k. -> tilt_divq.txt"""
import os, sys, glob, json, numpy as np
from multiprocessing import Pool
ROOT = "/home/baloewe/Physics/Overcharged"; H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H); sys.path.insert(0, ROOT + "/runs/cavity")
import bulk_spiral as BS  # noqa
from edge_criteria import read_with_angle  # noqa


def one(run):
    bc = json.load(open(run + "/input.json"))["BC"][0]; R = bc["R"]; k = bc["customAnch"]; Q = np.array(bc["Q"][:2], float)
    D = BS.frames(run + "/directorfield.dat.gz", 8); td = np.unique(D[:, 0]); td = td[td >= td[0] + 0.5 * (td[-1] - td[0])]
    t, fr = read_with_angle(run, Q, 0.0); t = np.array(t); n = int(D[:, 1].max()) + 1; m = int(D[:, 2].max()) + 1
    X, Y = np.meshgrid(np.arange(n) + 0.5 - Q[0], np.arange(m) + 0.5 - Q[1], indexing="ij"); TT, TD, DIFF = [], [], []
    for tt in td:
        d = D[D[:, 0] == tt]; q1 = np.zeros((n, m)); q2 = np.zeros((n, m)); th = np.arctan2(d[:, 5], d[:, 4])
        q1[d[:, 1].astype(int), d[:, 2].astype(int)] = d[:, 7] * np.cos(2 * th); q2[d[:, 1].astype(int), d[:, 2].astype(int)] = d[:, 7] * np.sin(2 * th)
        dvx = np.gradient(q1, axis=0) + np.gradient(q2, axis=1); dvy = np.gradient(q2, axis=0) - np.gradient(q1, axis=1)
        f = fr[int(np.argmin(np.abs(t - tt)))]
        if not len(f): continue
        r = np.hypot(f[:, 0], f[:, 1]); w = (f[:, 2] > 0) & (R - r < 4)
        if w.sum() < 3: continue
        phi = np.arctan2(f[w, 1], f[w, 0]); tau_t = f[w, 3] - phi; tau_d = []
        for x0, y0, p in zip(f[w, 0], f[w, 1], phi):
            mk = np.hypot(X - x0, Y - y0) < 2.0; tau_d.append(np.arctan2(dvy[mk].mean(), dvx[mk].mean()) - p)
        tau_d = np.array(tau_d); zt = np.mean(np.exp(1j * tau_t)); zd = np.mean(np.exp(1j * tau_d))
        if abs(zt) > 0.9: TT.append(abs(np.angle(zt))); TD.append(abs(np.angle(zd))); DIFF += list(np.angle(np.exp(1j * (tau_t - tau_d))) * np.sign(np.angle(zt)))
    return dict(k=k, a=2 * np.pi * R / (2 * k - 2), tt=np.degrees(np.median(TT)) if TT else np.nan, td=np.degrees(np.median(TD)) if TD else np.nan,
                diff=np.degrees(np.median(DIFF)) if DIFF else np.nan, sd=np.degrees(np.std(DIFF)) if DIFF else np.nan, nf=len(TT))


if __name__ == "__main__":
    runs = sorted(os.path.dirname(p) for p in glob.glob(ROOT + "/runs/cavity/server_main/main_act0.07/rho*/k*/s[01][0-9]/directorfield.dat.gz")
                  if abs(json.load(open(os.path.dirname(p) + "/input.json"))["BC"][0]["R"] - 31.83) < 0.05)
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 6) as pool: res = pool.map(one, runs)
    L = ["TILT: tracker angle vs div Q of the director field, same frames (R 31.83, coherent frames by the tracker). Per k: medians over seeds.",
         "k | spacing | tracker tilt | div-Q tilt | per-defect (tracker - divQ), signed with the tilt: median, sd | coherent director frames"]
    for kk in sorted(set(x["k"] for x in res)):
        S_ = [x for x in res if x["k"] == kk]
        L.append(f"{kk:2d} | {S_[0]['a']:5.1f} | {np.nanmedian([x['tt'] for x in S_]):5.1f} | {np.nanmedian([x['td'] for x in S_]):5.1f} | {np.nanmedian([x['diff'] for x in S_]):+5.1f} {np.nanmedian([x['sd'] for x in S_]):5.1f} | {sum(x['nf'] for x in S_)}")
    print("\n".join(L)); open(os.path.join(H, "tilt_divq.txt"), "w").write("\n".join(L) + "\n")
