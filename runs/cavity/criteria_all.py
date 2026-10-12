#!/usr/bin/env python3
"""
All edge/bulk criteria per circle run (main sweep, positive and mirror winding), rho = |k|/(pi R) <= --rhomax.

Sign convention: sw = sign(k) = sign of the imposed wall defects. Bulk charge q_b always in the original charges.
  knee      positive: knee of Q(r) (charge_knee.knee); mirror: plateau at the peak of Q(r) before the plunge to -|k|
            (peak value; "plateau" if the band within 0.1 of the peak is >= 0.04 in r/R wide and the plunge is >= 5x steeper)
  fixed     Q(r < 11 R / 12)
  site      pinned-site rule on the wall-sign defects (edge_sites.classify on sign-flipped frames for mirror runs)
  dens      minimum of the wall-sign defect density between its wall peak and the interior; bulk = Q(r < r_min)
  slow      wall-sign defects with displacement < 5 cells over 20 time units and within 6 cells of the wall are edge
  (slow for mirror runs: the wall's -1/2 followed with flipped signs)
  tilt      +k only: slow and |orientation - normal| < 1.8 rad (bulk +1/2 sit near 2.4 rad; the wall tilt drops from
            ~1.1 at rho 0.06 to ~0.85 at rho 0.15 and broadens, so no lower cut)
Also: r* = Q = +1 crossing (r/R), f* = 1 - r*^2, knee position r_knee/R.
Usage: criteria_all.py [--workers 12] [--out analysis_weekend/criteria_all.csv]
"""
import argparse, csv, glob, json, os, sys
from multiprocessing import Pool
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "cassini"))

COLS = ["set", "run", "k", "R", "rho", "q_knee", "knee_r", "has_knee", "q_fixed", "q_site", "nsite", "q_dens", "r_dens",
        "q_slow", "q_tilt", "rstar", "fstar"]


def one(run):
    from lobe_profiles import profile
    from charge_knee import knee
    from edge_sites import classify
    from defect_spatial import load
    from edge_criteria import analyse as ec
    try:
        bc = json.load(open(os.path.join(run, "input.json")))["BC"][0]
        k = int(bc["customAnch"]); R = bc["R"]; sw = 1 if k > 0 else -1
        p = profile(run, 200, 0.2); x = p["edges"][1:]; Q = p["Q"]
        row = dict(set="pos" if k > 0 else "neg", run=os.path.relpath(run, os.path.join(HERE, "server_main")), k=k, R=R, rho=abs(k) / (np.pi * R))
        # knee
        if k > 0:
            kn = knee(p)
            row.update(q_knee=kn[1] if kn else np.nan, knee_r=kn[0] if kn else np.nan, has_knee=int(kn is not None))
        else:
            ker = np.ones(9) / 9; Qs = np.convolve(np.r_[np.full(9, Q[0]), Q, np.full(9, Q[-1])], ker, "same")[9:-9]
            m = x < 0.97; i = int(np.argmax(np.where(m, Qs, -1e9))); qp = Qs[i]
            band = np.where(m & (Qs > qp - 0.1))[0]; width = x[band].max() - x[band].min() if len(band) else 0
            S = np.abs(np.gradient(Qs, x)); plunge = S[i:].max(); flat = S[band].mean() if len(band) else np.inf
            ok = width >= 0.04 and plunge >= 5 * max(flat, 1e-6) and qp > 0.05
            row.update(q_knee=qp if ok else np.nan, knee_r=x[i] if ok else np.nan, has_knee=int(ok))
        row["q_fixed"] = float(np.interp(11 / 12, x, Q))
        # r*, f*
        idx = np.where(Q < 1.0)[0]
        if k > 0 and len(idx) and idx[-1] + 1 < len(Q):
            j = idx[-1]; rs = x[j] + (1 - Q[j]) / (Q[j + 1] - Q[j]) * (x[j + 1] - x[j])
        else:
            above = np.where((Q >= 1.0) & (x < 0.97))[0]
            rs = x[above[0]] if (k < 0 and len(above)) else np.nan
        row.update(rstar=rs, fstar=1 - rs * rs if np.isfinite(rs) else np.nan)
        # density minimum of the wall-sign defects
        nw = p["npl"] if sw > 0 else p["nmi"]
        nws = np.convolve(nw, np.ones(5) / 5, "same")
        ip = int(np.argmax(np.where(x > 1 - 6 / R, nws, -1)))           # wall peak within 6 cells of the wall
        lo = np.searchsorted(x, max(1 - 10 / R, 0.3))                       # gap searched within 10 cells of the wall
        if ip > lo + 2:
            im = lo + int(np.argmin(nws[lo:ip])); row.update(q_dens=Q[im], r_dens=x[im])
        else:
            row.update(q_dens=np.nan, r_dens=np.nan)
        # site rule on the wall-sign defects
        d = load(run, 0.2)
        if sw < 0:
            for f in d["frames"]:
                f[:, 2] *= -1
        c = classify(d, 1.0, 0.02, 0.0, 4.2, 0.4)
        row.update(q_site=sw * c["qb"] if np.isfinite(c["qb"]) else np.nan, nsite=c["nsite"])
        # mobility / tilt (edge_criteria works on +1/2; flip for mirror via a temporary sign)
        if sw > 0:
            e = ec(run, 20, 8.0, 0.2)
            slow = (e["D"] < 5) & (e["S"] < 6); tilt = slow & (np.abs(e["Dl"]) < 1.8)          # upper cut only: the wall tilt shrinks and broadens at high rho
            n_slow = slow.sum() / e["nframes"]; n_tilt = tilt.sum() / e["nframes"]     # per frame used (unbiased)
            row.update(q_slow=k - 0.5 * n_slow, q_tilt=k - 0.5 * n_tilt)
        else:
            e = ec(run, 20, 8.0, 0.2, flip=True)
            slow = (e["D"] < 5) & (e["S"] < 6)
            row.update(q_slow=k + 0.5 * slow.sum() / e["nframes"], q_tilt=np.nan)   # bulk = k - (-1/2) x n_slow
        return row
    except Exception as ex:
        return dict(set="err", run=run, k=str(ex)[:80])


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=12); ap.add_argument("--rhomax", type=float, default=0.151)
    ap.add_argument("--out", default=os.path.join(HERE, "analysis_weekend", "criteria_all.csv"))
    a = ap.parse_args()
    runs = []
    for st in ("main_act0.07", "main_act0.07_neg"):
        for dn in sorted(glob.glob(os.path.join(HERE, "server_main", st, "rho*", "k*", "s0?", "DONE"))):
            r = os.path.dirname(dn); bc = json.load(open(os.path.join(r, "input.json")))["BC"][0]
            if abs(bc["customAnch"]) / (np.pi * bc["R"]) <= a.rhomax:
                runs.append(r)
    print(len(runs), "runs", flush=True)
    with Pool(a.workers) as pool, open(a.out, "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=COLS, extrasaction="ignore"); w.writeheader()
        for i, row in enumerate(pool.imap_unordered(one, runs, chunksize=2)):
            if row.get("set") == "err":
                print("FAILED", row, flush=True); continue
            w.writerow(row); fo.flush()
            if i % 50 == 0:
                print(i, row["run"], flush=True)
    print("written", a.out)


if __name__ == "__main__":
    main()
