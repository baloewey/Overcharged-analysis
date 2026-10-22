#!/usr/bin/env python3
"""Fig. 1(c) check (user, 2026-10-22: 'bulk charge stays 1 at large R but vacancies grow; charge does not add up'). 10-seed circle sweep (runs/cavity/server_main,
main_act0.07, +k, rho <= 0.15, k 2-10), seeds s01-s03, second half. Time-averaged +1/2 and -1/2 counts and charge in bands of wall distance d: 0-4 (the 'wall'
of Fig. 1c), 4-8, 8-12 and > 12 cells; total charge vs k; and the knee radius / Q_knee from criteria_all.csv for the same runs. -> charge_bands.txt"""
import csv, glob, json, os, sys, numpy as np
from multiprocessing import Pool
ROOT = "/home/baloewe/Physics/Overcharged/runs/cavity"; sys.path.insert(0, ROOT)
from edge_criteria import read_with_angle  # noqa
LA = 22.21; B = [0, 4, 8, 12, 1e9]
def one(run):
    bc = json.load(open(run + "/input.json"))["BC"][0]; R = bc["R"]; k = bc["customAnch"]; Q = np.array(bc["Q"][:2], float)
    t, fr = read_with_angle(run, Q, 0.5); acc = np.zeros((len(B) - 1, 2))
    for f in fr:
        if not len(f): continue
        d = R - np.hypot(f[:, 0], f[:, 1]); b = np.clip(np.digitize(d, B) - 1, 0, len(B) - 2)
        acc[:, 0] += np.bincount(b[f[:, 2] > 0], minlength=len(B) - 1); acc[:, 1] += np.bincount(b[f[:, 2] < 0], minlength=len(B) - 1)
    return run, k, R / LA, acc / max(len(fr), 1)
if __name__ == "__main__":
    C = {r["run"]: r for r in csv.DictReader(open(ROOT + "/analysis_weekend/criteria_all.csv")) if r["set"] == "pos"}
    runs = [os.path.dirname(f) for f in glob.glob(ROOT + "/server_main/main_act0.07/rho*/k*/s0[123]/DONE")]
    runs = [r for r in runs if 2 <= json.load(open(r + "/input.json"))["BC"][0]["customAnch"] <= 10]
    with Pool(8) as p: res = p.map(one, runs, chunksize=2)
    G = {}
    for run, k, R, a in res: G.setdefault((k, round(R, 2)), []).append((a, run))
    L = ["CHARGE BY WALL-DISTANCE BAND, circles (main_act0.07, seeds 1-3): k | R/l_a | rho | wall +1/2 (0-4) - (2k-2) | +1/2, -1/2 at 4-8 | charge 4-8 | charge 8-12 | charge > 12 | "
         "bulk charge outside 4 cells (k - wall/2 check) | total charge | Q_knee (mean of seeds with a knee) | knee radius from wall"]
    for (k, R), v in sorted(G.items()):
        if R < 1.9 and not (1.40 < R < 1.46): continue
        a = np.mean([x for x, _ in v], 0); q = 0.5 * (a[:, 0] - a[:, 1]); P = 2 * np.pi * R * LA
        kn = [C.get(os.path.relpath(r, ROOT + "/server_main")) for _, r in v]; kn = [c for c in kn if c and c["has_knee"] in ("1", "True", "true", "1.0")]
        qk = np.mean([float(c["q_knee"]) for c in kn]) if kn else np.nan; kr = np.mean([R * LA - float(c["knee_r"]) for c in kn]) if kn else np.nan
        L.append(f"{k:2d} | {R:.2f} | {2 * k / P:.3f} | {a[0, 0] - (2 * k - 2):+5.2f} | {a[1, 0]:.2f}, {a[1, 1]:.2f} | {q[1]:+.2f} | {q[2]:+.2f} | {q[3]:+.2f} | {q[1:].sum():+.2f} | {q.sum():+.2f} | {qk:.2f} ({len(kn)}/{len(v)}) | {kr:.1f}")
    open("charge_bands.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
