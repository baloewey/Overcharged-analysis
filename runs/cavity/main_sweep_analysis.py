#!/usr/bin/env python3
"""
Main sweep (server circles, runs/cavity/server_main/): per-run layer diagnostics, written row by row to a CSV.

Per run (first 20% of frames dropped): R, k, rho = k/(pi R), seed, act;
  knee of Q(r) (cassini/charge_knee.knee): Q_knee (bulk charge) or none;
  Q = 1 crossing: c*/R, d* = R - c*, f* = 1 - (c*/R)^2;
  site rule (edge_sites.classify): number of sites, q_b, P(q_b = +1);
  wall Fourier (wall_fourier.analyse, band 4): m*, S(m*), f_static, N_wall.
Restartable: runs already in the CSV are skipped.
Usage: main_sweep_analysis.py [--root server_main] [--out analysis_weekend/main_sweep.csv] [--sets main_act0.07 main_passive]
"""
import argparse, csv, glob, json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "cassini"))
from lobe_profiles import profile          # noqa: E402
from charge_knee import knee               # noqa: E402
from edge_sites import classify            # noqa: E402
from defect_spatial import load            # noqa: E402
from wall_fourier import analyse as wf     # noqa: E402

COLS = ["set", "run", "R", "k", "rho", "seed", "act", "frames", "Q_knee", "knee_c", "knee_contrast", "cstar", "dstar", "fstar",
        "nsite", "qb_site", "P1_site", "m_star", "S_star", "f_static", "N_wall"]


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--root", default=os.path.join(HERE, "server_main"))
    ap.add_argument("--out", default=os.path.join(HERE, "analysis_weekend", "main_sweep.csv"))
    ap.add_argument("--sets", nargs="+", default=["main_act0.07", "main_passive"])
    a = ap.parse_args()
    done = set()
    if os.path.exists(a.out):
        done = {r["run"] for r in csv.DictReader(open(a.out))}
    new = not os.path.exists(a.out)
    fo = open(a.out, "a", newline=""); w = csv.DictWriter(fo, fieldnames=COLS)
    if new:
        w.writeheader()
    for st in a.sets:
        runs = sorted(os.path.dirname(p) for p in glob.glob(os.path.join(a.root, st, "rho*", "k*", "s0?", "DONE")))
        for run in runs:
            rel = os.path.relpath(run, a.root)
            if rel in done:
                continue
            try:
                j = json.load(open(os.path.join(run, "input.json"))); bc = j["BC"][0]
                R = bc["R"]; k = bc.get("customAnch"); act = j["species"][0].get("act", 0.0)
                p = profile(run, 200, 0.2)
                kn = knee(p)
                x = p["edges"][1:]; Q = p["Q"]
                idx = np.where(Q < 1.0)[0]
                cs = np.nan
                if len(idx) and idx[-1] + 1 < len(Q):
                    i = idx[-1]; cs = x[i] + (1 - Q[i]) / (Q[i + 1] - Q[i]) * (x[i + 1] - x[i])
                d = load(run, 0.2)
                c = classify(d, 1.0, 0.02, 0.0, 4.2, 0.4)
                f = wf(run, 4.0, 0.2) if k and abs(k) > 0 else None
                row = dict(set=st, run=rel, R=R, k=k, rho=abs(k) / (np.pi * R), seed=j.get("seed"), act=act, frames=len(d["frames"]),
                           Q_knee=kn[1] if kn else "", knee_c=kn[0] if kn else "", knee_contrast=kn[2] if kn else "",
                           cstar=cs, dstar=R * (1 - cs) if np.isfinite(cs) else "", fstar=1 - cs * cs if np.isfinite(cs) else "",
                           nsite=c["nsite"], qb_site=c["qb"], P1_site=c["P1"],
                           m_star=f["mstar"] if f else "", S_star=f["Sstar"] if f else "", f_static=f["fstar"] if f else "", N_wall=f["Nwall"] if f else "")
                w.writerow(row); fo.flush()
                print(f"{rel}: k {k} R {R:.1f} knee {row['Q_knee']} d* {row['dstar']} sites {c['nsite']} f_static {row['f_static']}", flush=True)
            except Exception as e:
                print(f"{rel}: FAILED {e}", flush=True)
    fo.close()


if __name__ == "__main__":
    main()
