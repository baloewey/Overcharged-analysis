#!/usr/bin/env python3
"""U-scan escape test (2026-10-27). Aspect-2 ellipses, R_eq 1.5 l_a, k 5 (rho 0.04) and 14 (rho 0.12), act 0.07, at mfpot 5 and 20 (runs/Uscan,
3 seeds) vs mfpot 10 (shape sweep, 10 seeds). Frank K (fluctuation splay, runs/Kcalib + Uscan/Kcalib_U*): 8, 43, 88 kBT at U 5, 10, 20.
Per run, all with the validated tracker (prediction + Hungarian, gate 9):
  p_esc   escape probability per wall +1/2 per 20 t.u. (shape_esc2.one: wall d < 4 at a start frame, deeper than 6 at +20 t.u.)
  S(x)    near-wall survival of resident excursions (drive_barrier.one: resident 20 t.u., deepest point in the next 20 t.u.)
  N_wall  mean number of wall +1/2 (d < 4); Q_b bulk charge (d > 6)
Question: does the elastic barrier set the escape rate (near-wall law S(x) ~ exp(-U(x)/T) with U proportional to K)?
-> uscan_escape.csv, printed summary"""
import csv, glob, os, sys, numpy as np
from multiprocessing import Pool
H = os.path.dirname(os.path.abspath(__file__)); CP = os.path.join(H, "..", "crowd_protect")
sys.path[:0] = [os.path.join(CP, "exchange"), os.path.join(CP, "drive_barrier"), os.path.join(CP, "tracker_validation")]
import exchange_count as X, shape_esc2 as SE, drive_barrier as DB  # noqa
ROOT = X.ROOT


def occupancy(run):
    kind, Q, dist, Lw, info = X.geometry(run); t, fr = X.read_with_angle(run, Q, 0.5); nw, qb = [], []
    for f in fr[len(fr) // 10:]:
        if not len(f): continue
        d = dist(f[:, :2]); nw.append(((d < 4) & (f[:, 2] > 0)).sum()); qb.append(0.5 * np.sign(f[d > 6, 2]).sum())
    return dict(k=info["k"], nwall=np.mean(nw), qb=np.mean(qb))


def one(run):
    r = dict(run=run.replace(ROOT + "/", ""))
    for f in (SE.one, DB.one, occupancy):
        o = f(run)
        if "err" in o: return dict(r, err=o["err"])
        o.pop("run", None); r.update(o)
    return r


if __name__ == "__main__":
    runs = [(5, r) for r in sorted(glob.glob(ROOT + "/Uscan/ell2_k*_U5_s*"))] + [(20, r) for r in sorted(glob.glob(ROOT + "/Uscan/ell2_k*_U20_s*"))]
    runs += [(10, os.path.dirname(f)) for f in sorted(glob.glob(ROOT + "/shapes_server/ellipse/eps2.00/Req1.5/k0[5]/s*/DONE") + glob.glob(ROOT + "/shapes_server/ellipse/eps2.00/Req1.5/k14/s*/DONE"))]
    with Pool(int(os.environ.get("NPROC", 12))) as p: res = p.map(one, [r for _, r in runs])
    for (U, _), r in zip(runs, res): r["U"] = U
    cols = ["run", "U", "k", "starts", "esc", "p_esc", "p_late", "n", "d0"] + [f"S{x:g}" for x in DB.XG] + ["nwall", "qb", "err"]
    with open(os.path.join(H, "uscan_escape.csv"), "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=cols, extrasaction="ignore"); w.writeheader(); w.writerows(res)
    print("runs", len(res), "errors", [r["run"] for r in res if "err" in r])
