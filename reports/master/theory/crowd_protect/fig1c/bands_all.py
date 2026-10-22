#!/usr/bin/env python3
"""Per-run wall-distance profiles for the Fig. 1(c) alternatives (2026-10-22). 10-seed circle sweep (runs/cavity/server_main/main_act0.07, +k, k 2-10,
rho <= 0.15), second half: time-averaged number of +1/2 and -1/2 per frame in 0.5-cell bins of the wall distance d (0-40 cells, last bin open).
-> bands_all.npz (run, k, R/l_a, rho, Hp, Hm, edges) + the knee and N_wall from runs/cavity/analysis_weekend."""
import glob, json, os, sys, numpy as np
from multiprocessing import Pool
ROOT = "/home/baloewe/Physics/Overcharged/runs/cavity"; sys.path.insert(0, ROOT)
from edge_criteria import read_with_angle  # noqa
LA = 22.21; E = np.r_[np.arange(0, 40.01, 0.5), 1e9]
def one(run):
    bc = json.load(open(run + "/input.json"))["BC"][0]; R = bc["R"]; k = bc["customAnch"]; Q = np.array(bc["Q"][:2], float)
    t, fr = read_with_angle(run, Q, 0.5); hp = np.zeros(len(E) - 1); hm = np.zeros(len(E) - 1)
    for f in fr:
        if not len(f): continue
        d = R - np.hypot(f[:, 0], f[:, 1]); hp += np.histogram(d[f[:, 2] > 0], E)[0]; hm += np.histogram(d[f[:, 2] < 0], E)[0]
    n = max(len(fr), 1); return os.path.relpath(run, ROOT + "/server_main"), k, R / LA, 2 * k / (2 * np.pi * R), hp / n, hm / n
if __name__ == "__main__":
    runs = [os.path.dirname(f) for f in glob.glob(ROOT + "/server_main/main_act0.07/rho*/k*/s*/DONE")]
    runs = [r for r in runs if 2 <= json.load(open(r + "/input.json"))["BC"][0]["customAnch"] <= 10]
    with Pool(8) as p: res = p.map(one, runs, chunksize=4)
    np.savez_compressed("bands_all.npz", run=np.array([r[0] for r in res]), k=np.array([r[1] for r in res]), R=np.array([r[2] for r in res]), rho=np.array([r[3] for r in res]),
                        Hp=np.array([r[4] for r in res]), Hm=np.array([r[5] for r in res]), edges=E)
    print(len(res), "runs")
