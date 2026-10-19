#!/usr/bin/env python3
"""Wall -1/2 tilt, lock and charge budget in the mirror shape sweep (user, 2026-10-19). Per run (second half, every 5th frame):
  charge budget: wall band (<4 cells) -1/2 and +1/2, mid band (4-6), bulk (>6); total vs k
  -1/2 tilt: raw angle tau_raw = solver angle - outward normal angle; candidate lock relations m*tau_raw - c*psi_rel (psi_rel = anchoring relative to the
  normal, as in shape_sweep_analysis); the one with concentration near 1 fixes the convention. Prediction (charge s, director theta = s alpha + theta0, arm
  where theta = alpha): tau = psi_rel / (1 - s) mod pi/(1-s); +1/2: tau = 2 psi_rel; -1/2: 3 tau = 2 psi_rel mod 2pi -> commensurate ring 2|k-1| for both signs.
  For -1/2 the tilt is defined mod 120 deg: per frame tau_c = arg<exp(3 i tau)>/3 in (-60, 60] (common tilt), |tau_c| averaged; coherence |<exp(3 i tau)>|.
  tau = 0: one arm into the wall (homeotropic point, rest); +-60: one arm into the bulk. lock_phase = arg<exp(i(3 tau - 2 psi_rel))> (convention offset, expect 0).
Usage: mirror_tilt.py [test] -> mirror_tilt.csv"""
import csv, glob, json, os, sys
import numpy as np
from multiprocessing import Pool
from scipy.spatial import cKDTree
ROOT = "/home/baloewe/Physics/Overcharged/runs"; sys.path.insert(0, ROOT + "/cavity")
from edge_criteria import read_with_angle  # noqa
from wall_fourier import wall_polyline  # noqa
SW = ROOT + "/shapes_mirror"; H = os.path.dirname(os.path.abspath(__file__))
CANDS = [(m, c) for m in (1, -1, 3, -3) for c in (2, -2, 1, -1, 6, -6)]
def one(run):
    try:
        bc = json.load(open(run + "/input.json"))["BC"][0]; k = bc["customAnch"]; Q = bc["Q"]
        W, s, P = wall_polyline(bc); Wc = W[:, 0] + 1j * W[:, 1]; T = np.roll(Wc, -1) - np.roll(Wc, 1); Nn = -1j * T / np.abs(T)
        phin = np.unwrap(np.angle(Nn)); psi = k * 2 * np.pi * s / P - phin; tree = cKDTree(W)
        t, fr = read_with_angle(run, np.array(Q[:2], float), 0.5); fr = fr[::5]
        bud = np.zeros(5); cand = np.zeros(len(CANDS), complex); nc = 0; coh, tilt = [], []
        for f in fr:
            if not len(f): continue
            d, i = tree.query(f[:, :2]); q = f[:, 2]; wl = d < 4
            bud += [np.sum(wl & (q < 0)), np.sum(wl & (q > 0)), 0.5 * np.sum(np.sign(q[(d >= 4) & (d <= 6)])), 0.5 * np.sum(np.sign(q[d > 6])), 0.5 * np.sum(np.sign(q))]
            m = wl & (q < 0)
            if not m.any(): continue
            tr = f[m, 3] - np.angle(Nn[i[m]]); pr = psi[i[m]]
            for j, (a, c) in enumerate(CANDS): cand[j] += np.sum(np.exp(1j * (a * tr - c * pr)))
            nc += m.sum(); z = np.exp(3j * tr); coh.append(abs(z.mean())); tilt.append(np.angle(z.mean()) / 3)
        n = len(fr); cand_raw = cand.copy(); cand = np.abs(cand) / max(nc, 1); jb = int(np.argmax(cand))
        return dict(run=os.path.relpath(run, SW), k=k, n_wall_minus=bud[0] / n, n_wall_plus=bud[1] / n, q_mid=bud[2] / n, q_bulk=bud[3] / n, q_total=bud[4] / n,
                    best_m=CANDS[jb][0], best_c=CANDS[jb][1], best_conc=cand[jb], conc_3_2=cand[CANDS.index((3, 2))], conc_3_m2=cand[CANDS.index((3, -2))],
                    conc_m3_2=cand[CANDS.index((-3, 2))], tilt_abs_deg=np.degrees(np.mean(np.abs(tilt))) if tilt else np.nan, lock_phase_deg=np.degrees(np.angle(cand_raw[CANDS.index((3, 2))])), coh3=np.mean(coh) if coh else np.nan)
    except Exception as e:
        return dict(run=os.path.relpath(run, SW), error=str(e))
if __name__ == "__main__":
    runs = sorted(os.path.dirname(f) for f in glob.glob(SW + "/*/*/*/*/s*/DONE")) + sorted(os.path.dirname(f) for f in glob.glob(ROOT + "/neglayer/*/DONE"))
    if len(sys.argv) > 1: runs = [r for r in runs if "Req1.5/km09/s01" in r or ("eps2.00/Req1.5" in r and "s01" in r)][:3]
    with Pool(16) as p: res = p.map(one, runs, chunksize=2)
    if len(sys.argv) > 1:
        for r in res: print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items()})
        sys.exit()
    man = {r["dir"]: r for r in csv.DictReader(open(SW + "/manifest_shape_mirror.tsv"), delimiter="\t")}
    cols = ["run", "family", "aspect", "Req_la", "k", "rho"] + [c for c in res[0] if c not in ("run", "k")]
    with open(H + "/mirror_tilt.csv", "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=cols, extrasaction="ignore"); w.writeheader()
        for o in res:
            m = man.get(o["run"], {}); o.update(family=m.get("family"), aspect=m.get("aspect"), Req_la=m.get("Req_la"), rho=m.get("rho")); w.writerow(o)
    print("errors", sum("error" in o for o in res))
