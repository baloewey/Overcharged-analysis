#!/usr/bin/env python3
"""Item 4 of the 27 Oct list: recheck the analyses built on old linkers with the validated tracker (NEW v2; tracker_validation/trackers.track_new, gate 9).
Per shape-sweep run (second half):
  KINKS (was kink_life.py, NN linking): per tracked wall +1/2 (d < 4), episodes with the tilt sign opposite to the frame majority: start rate per defect per
        1000 t.u. and mean duration (t.u.); wrong fraction (tracking-free).
  DEPTH MEMORY (was depth_memory.py, NN linking): normalised autocorrelation of resident depth and cos tilt along resident stretches (>= 60 frames), 1/e times.
  REFILL + ESCAPE REGIONS (was kinetics_carriers.py, 10 t.u. chain): refill probability per bulk +1/2 (d > 6) per 20 t.u. (in the wall, d < 4, at +20 t.u.);
        escape (d < 4 -> d > 6 at +20 t.u.) by start region |x|/A: tip > 1/2, flank < 1/4, middle; share of escapes and probability per region.
-> recheck4.csv"""
import csv, glob, json, os, sys, numpy as np
from multiprocessing import Pool
sys.path.insert(0, "../tracker_validation"); sys.path.insert(0, "../exchange")
from trackers import track_new  # noqa
import exchange_count as X  # noqa
from wall_fourier import wall_polyline  # noqa
from scipy.spatial import cKDTree
ROOT = X.ROOT
def acf_time(series_list, maxlag=30, dt=2.0):
    num = np.zeros(maxlag + 1); den = np.zeros(maxlag + 1)
    for a in series_list:
        a = a - a.mean()
        for L in range(maxlag + 1):
            if len(a) > L: num[L] += np.sum(a[:len(a) - L] * a[L:]); den[L] += len(a) - L
    c = num / np.maximum(den, 1); c = c / c[0] if c[0] > 0 else c; b = np.flatnonzero(c < np.exp(-1)); return (b[0] * dt) if len(b) else np.nan
def one(run):
    try:
        bc = json.load(open(run + "/input.json"))["BC"][0]; W, s, P = wall_polyline(bc); tree = cKDTree(W); A = np.abs(W[:, 0]).max()
        Wc = W[:, 0] + 1j * W[:, 1]; Tg = np.roll(Wc, -1) - np.roll(Wc, 1); Tg /= np.abs(Tg); Nn = -1j * Tg
        if np.mean((Nn * np.conj(Wc - Wc.mean())).real) < 0: Nn = -Nn
        t, fr = X.read_with_angle(run, np.array(bc["Q"][:2], float), 0.5); nf = len(fr); dt = t[1] - t[0]
        maj = []
        for f in fr:
            p = f[f[:, 2] > 0] if len(f) else np.zeros((0, 4)); d, j = tree.query(p[:, :2]) if len(p) else (np.zeros(0), np.zeros(0, int)); w = d < 4
            tl = np.angle(np.exp(1j * (p[w, 3] - np.angle(Nn[j[w]])))) if w.any() else np.zeros(0); maj.append(np.sign(np.sign(tl).sum()) or 1.0)
        T = [tr for tr in track_new(fr, 9.0) if tr["sg"] > 0]
        eps = []; nres = 0; dep_s = []; cos_s = []; st_b = rf = 0; reg = np.zeros((3, 2))
        for tr in T:
            nn = np.array(tr["n"]); xy = np.array(tr["xy"]); d, j = tree.query(xy); pos = dict(zip(nn, range(len(nn))))
            ang = []
            for n, (x, y) in zip(nn, xy):                              # orientation of this detection
                f = fr[n]; k = np.argmin(np.hypot(f[:, 0] - x, f[:, 1] - y)); ang.append(f[k, 3])
            tl = np.angle(np.exp(1j * (np.array(ang) - np.angle(Nn[j])))); wrong = (np.sign(tl) != np.array([maj[n] for n in nn])) & (d < 4)
            run_ = 0
            for i in range(len(nn)):
                if d[i] < 4: nres += 1
                if wrong[i] and (i == 0 or nn[i] == nn[i - 1] + 1): run_ += 1
                else:
                    if run_: eps.append(run_)
                    run_ = 1 if wrong[i] else 0
            if run_: eps.append(run_)
            res = d < 4; idx = np.flatnonzero(res); 
            for seg in np.split(idx, np.flatnonzero(np.diff(idx) > 1) + 1):
                if len(seg) >= 60 and np.all(np.diff(nn[seg]) == 1): dep_s.append(d[seg]); cos_s.append(np.cos(tl[seg]))
            for i, n in enumerate(nn):
                if n % 5 or n + 10 >= nf: continue
                k10 = pos.get(n + 10)
                if d[i] > 6: st_b += 1; rf += k10 is not None and d[k10] < 4
                elif d[i] < 4:
                    r_ = 0 if abs(xy[i, 0]) / A > 0.5 else (2 if abs(xy[i, 0]) / A < 0.25 else 1); reg[r_, 0] += 1; reg[r_, 1] += k10 is not None and d[k10] > 6
        e = np.array(eps) * dt; Tt = nf * dt
        return dict(run=os.path.relpath(run, ROOT + "/shapes_server"), kink_rate=1000 * len(eps) / max(nres * dt, 1), kink_life=e.mean() if len(e) else np.nan,
                    acf_depth=acf_time(dep_s, dt=dt), acf_cos=acf_time(cos_s, dt=dt), p_ref=rf / max(st_b, 1),
                    esc_tip=reg[0, 1] / max(reg[0, 0], 1), esc_mid=reg[1, 1] / max(reg[1, 0], 1), esc_flank=reg[2, 1] / max(reg[2, 0], 1), tip_share_esc=reg[0, 1] / max(reg[:, 1].sum(), 1))
    except Exception as ex:
        return dict(run=os.path.relpath(run, ROOT + "/shapes_server"), err=str(ex))
if __name__ == "__main__":
    runs = sorted(os.path.dirname(f) for f in glob.glob(ROOT + "/shapes_server/*/*/*/*/s*/DONE"))
    with Pool(int(os.environ.get("NPROC", 10))) as p: res = p.map(one, runs, chunksize=2)
    cols = ["run", "kink_rate", "kink_life", "acf_depth", "acf_cos", "p_ref", "esc_tip", "esc_mid", "esc_flank", "tip_share_esc", "err"]
    with open("recheck4.csv", "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=cols, extrasaction="ignore"); w.writeheader(); w.writerows(res)
    print("runs", len(res), "errors", sum(1 for r in res if "err" in r), [r.get("err") for r in res if "err" in r][:2])
