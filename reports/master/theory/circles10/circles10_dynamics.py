#!/usr/bin/env python3
"""Defect dynamics on the 10-seed circle sweep (user, 2026-10-15: extend the circle analyses to 10 seeds). Main sweep, act 0.07, rho <= 0.15, +k (k >= 2) and
mirror -k, seeds s01-s10 (runs/cavity/server_main). Per run (second half of the run unless stated), reusing the original per-run functions unchanged:
  ddg        bulkdyn/ddg.one             entries into / exits from the wall band (+1/2, hysteresis 4/6 cells), delay entry -> next exit vs null, same/other defect
  life       bulkdyn/pair_lifetime.one   pair-born -1/2 lifetimes and what they annihilate with
  nuc        energetics/nucleation.one   pair creations per bulk area, bulk -1/2 density, bulk charge Q_b
  exch       as runs/cavity/exchange_rates.py (first 20% dropped; wall < 3, bulk > 6 cells): wall +1/2 escape rate per wall defect, residence, arrivals
Writes circles10_dynamics.pkl (one dict per run; incremental: runs already in the pickle are skipped). Aggregation: circles10_summary.py.
Usage: circles10_dynamics.py [NPROC=8]"""
import os, sys, glob, json, pickle, numpy as np
from multiprocessing import Pool
R_ = "/home/baloewe/Physics/Overcharged/runs"; H = os.path.dirname(os.path.abspath(__file__)); T_ = os.path.dirname(H)
for p in (R_ + "/cavity", T_ + "/bulkdyn", T_ + "/energetics", T_ + "/kickout"): sys.path.insert(0, p)
import ddg, pair_lifetime, nucleation  # noqa
from defect_tracks import build_tracks  # noqa
from edge_criteria import read_with_angle  # noqa
from wall_fourier import wall_polyline, project  # noqa
OUT = os.path.join(H, "circles10_dynamics.pkl")


def exch(run, swall=3.0, sbulk=6.0, tfrac=0.2):
    bc = json.load(open(os.path.join(run, "input.json")))["BC"][0]
    t, fr = read_with_angle(run, np.array(bc["Q"][:2]), tfrac); dt = t[1] - t[0]; T = len(fr) * dt
    W, s, P = wall_polyline(bc); tracks, info = build_tracks(fr, dt, 8.0, 4)
    res, nleave, narrive, wall_occ = [], 0, 0, 0
    for tr in tracks.values():
        if tr["q"] < 0: continue
        d, _ = project(W, s, tr["xy"]); state = np.empty(len(d), int); cur = 1 if d[0] < swall else (0 if d[0] > sbulk else -1)
        for i, di in enumerate(d):
            cur = 1 if di < swall else (0 if di > sbulk else cur); state[i] = cur
        wall_occ += np.sum(state == 1); i = 0; n = len(state)
        while i < n:
            if state[i] == 1:
                j = i
                while j + 1 < n and state[j + 1] == 1: j += 1
                a_in = i > 0 and state[i - 1] == 0; a_out = j + 1 < n and state[j + 1] == 0
                narrive += a_in; nleave += a_out
                if a_in and a_out: res.append((tr["n"][j] - tr["n"][i] + 1) * dt)
                i = j + 1
            else: i += 1
    nwall = wall_occ / len(fr)
    return dict(nwall=nwall, escape=nleave / max(nwall, 1e-9) / T, arrive=narrive / T, res_med=np.median(res) if res else np.nan, nvis=len(res), T=T)


def one(run):
    J = json.load(open(run + "/input.json")); bc = J["BC"][0]
    r = dict(run=os.path.relpath(run, R_), k=bc["customAnch"], R=bc["R"], seed=os.path.basename(run))
    for name, f in (("ddg", ddg.one), ("life", pair_lifetime.one), ("nuc", nucleation.one), ("exch", exch)):
        try: r[name] = f(run)
        except Exception as ex: r[name] = dict(err=repr(ex))
    return r


if __name__ == "__main__":
    runs = []
    for st in ("main_act0.07", "main_act0.07_neg"):
        for d in sorted(glob.glob(f"{R_}/cavity/server_main/{st}/rho*/k*/s[01][0-9]/DONE")):
            d = os.path.dirname(d); rho = float(d.split("/rho")[1].split("/")[0]); k = abs(json.load(open(d + "/input.json"))["BC"][0]["customAnch"])
            if rho <= 0.151 and k >= 2: runs.append(d)
    done = {}
    if os.path.exists(OUT): done = {x["run"]: x for x in pickle.load(open(OUT, "rb"))}
    todo = [r for r in runs if os.path.relpath(r, R_) not in done]
    print(len(runs), "runs,", len(todo), "to do", flush=True)
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 8) as pool:
        for i, x in enumerate(pool.imap_unordered(one, todo, chunksize=2)):
            done[x["run"]] = x
            if i % 50 == 0:
                pickle.dump(list(done.values()), open(OUT, "wb")); print(i, x["run"], flush=True)
    pickle.dump(list(done.values()), open(OUT, "wb")); print("done", len(done))
