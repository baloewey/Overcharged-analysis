#!/usr/bin/env python3
"""
Flow and wall-defect tilt per circle run (main sweep, rho <= 0.15; flow reduced on the server by flow_reduce.py).

Flow (frames every 100 time units, first 20% dropped):
  vth_abs   time average of |<v_theta>|  (<.> = spatial mean over cells with r < R - 1; |.| removes the sign
            swaps of the circulation between runs and in time)
  vrms      time average of sqrt(<|v|^2>)
  ens       time average of <omega^2 / 2>
  nswap     number of sign changes of <v_theta> per 1000 time units
Tilt of the wall-sign defects (within 6 cells of the wall; defect frames every 2 time units, first 20% dropped):
  Delta = psi - phi, psi = defect orientation from defects.dat (Vromans-Giomi angle), phi = polar angle of the
  defect position (= outward wall normal). +1/2: wrapped to (-pi, pi]; Delta = 0: the radial arm of the comet
  (the line of the director field that leaves the core along psi) points at the wall. -1/2: three-fold, wrapped
  to (-pi/3, pi/3]; Delta = 0: one of the three radial arms points at the wall.
  tilt_abs   mean |Delta| over all wall-defect occurrences
  tilt_coh   time average of |frame circular mean of Delta|
  tilt_flow  mean over flow frames of sign(<v_theta>) x frame mean of Delta (same time stamps): the tilt in the
             frame where the circulation is counter-clockwise; agree = fraction with the same sign
Usage: flow_tilt.py [--workers 12]
"""
import argparse, csv, glob, json, os, sys
from multiprocessing import Pool
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
OUT = os.path.join(HERE, "analysis_flow")
COLS = ["set", "run", "k", "R", "rho", "vth_abs", "vth_signed", "vrms", "ens", "nswap", "nwall", "tilt_abs", "tilt_coh",
        "tilt_flow", "agree"]


def one(run):
    from edge_criteria import read_with_angle, wall
    try:
        full = os.path.join(HERE, "server_main", run)
        bc = json.load(open(os.path.join(full, "input.json")))["BC"][0]
        k = int(bc["customAnch"]); R = bc["R"]; sw = 1 if k > 0 else -1; m = 1 if k > 0 else 3
        row = dict(set="pos" if k > 0 else "neg", run=run, k=k, R=R, rho=abs(k) / (np.pi * R))
        fl = np.load(os.path.join(OUT, "flow_red", run.replace("/", "__") + ".npz"))
        n0 = int(0.2 * len(fl["t"])); t, vth = fl["t"][n0:], fl["vth"][n0:]
        row.update(vth_abs=np.mean(np.abs(vth)), vth_signed=np.mean(vth), vrms=np.mean(fl["vrms"][n0:]), ens=np.mean(fl["ens"][n0:]),
                   nswap=np.sum(np.diff(np.sign(vth)) != 0) / (t[-1] - t[0]) * 1000)
        td, fr = read_with_angle(full, np.array(bc["Q"][:2], float), 0.2)
        allD, fmean = [], {}
        hist = np.zeros(72)
        for ti, f in zip(td, fr):
            p = np.where(np.sign(f[:, 2]) == sw)[0]
            if not len(p):
                continue
            s, nang = wall(bc, f[p, 0], f[p, 1])
            q = s < 6
            if not q.any():
                continue
            z = np.exp(1j * m * (f[p[q], 3] - nang[q]))
            D = np.angle(z) / m; allD.append(D)
            zm = z.mean(); fmean[round(ti)] = np.angle(zm) / m if abs(zm) > 1e-6 else 0.0
        if allD:
            A = np.concatenate(allD)
            hist = np.histogram(A, 72, (-np.pi, np.pi))[0]
            F = np.array(list(fmean.values()))
            pairs = [(np.sign(v), fmean[round(tt)]) for tt, v in zip(t, vth) if round(tt) in fmean and v != 0]
            tf = np.mean([sg * d for sg, d in pairs]) if pairs else np.nan
            ag = np.mean([sg * d > 0 for sg, d in pairs]) if pairs else np.nan
            row.update(nwall=len(A) / len(fr), tilt_abs=np.mean(np.abs(A)), tilt_coh=np.mean(np.abs(F)), tilt_flow=tf, agree=ag)
        else:
            row.update(nwall=0, tilt_abs=np.nan, tilt_coh=np.nan, tilt_flow=np.nan, agree=np.nan)
        return row, hist
    except Exception as ex:
        return dict(set="err", run=run, k=str(ex)[:100]), None


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=12); a = ap.parse_args()
    runs = [l.strip() for l in open(os.path.join(OUT, "runlist.txt")) if l.strip()]
    H = {}
    with Pool(a.workers) as pool, open(os.path.join(OUT, "flow_tilt.csv"), "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=COLS, extrasaction="ignore"); w.writeheader()
        for i, (row, h) in enumerate(pool.imap_unordered(one, runs, chunksize=2)):
            if row["set"] == "err":
                print("FAILED", row, flush=True); continue
            w.writerow(row); H[row["run"]] = h
            if i % 100 == 0: print(i, flush=True)
    np.savez(os.path.join(OUT, "tilt_hist.npz"), runs=np.array(list(H)), hist=np.array(list(H.values())))
    print("written")


if __name__ == "__main__":
    main()
