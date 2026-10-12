#!/usr/bin/env python3
"""
Wall <-> bulk exchange of +1/2 defects (H7): residence times at the wall and arrival rates, from gap-closed tracks.

Track +1/2 (defect_tracks.build_tracks, rlink 8, gap 4); state per frame: wall (distance to the wall < --swall) or
bulk (> --sbulk); in between keeps the previous state (hysteresis). A wall visit = a maximal run of wall states.
  residence: duration of visits that start and end inside the track (left- and right-censored ones excluded)
  leave rate  = visits ending with a move to the bulk, per wall defect per time unit
  arrive rate = bulk -> wall transitions per time unit (whole cavity)
Also: fraction of wall visits that end by annihilation / tracker loss (track end) rather than leaving.
Usage: exchange_rates.py RUN [RUN ...] [--labels ...] [--swall 3] [--sbulk 6]
"""
import argparse, json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from defect_tracks import build_tracks          # noqa: E402
from edge_criteria import read_with_angle       # noqa: E402
from wall_fourier import wall_polyline, project  # noqa: E402


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("runs", nargs="+"); ap.add_argument("--labels", nargs="+")
    ap.add_argument("--swall", type=float, default=3.0); ap.add_argument("--sbulk", type=float, default=6.0)
    ap.add_argument("--tfrac", type=float, default=0.2)
    a = ap.parse_args()
    labels = a.labels or [os.path.basename(os.path.normpath(r)) for r in a.runs]
    print(f"{'run':<30} {'k':>3} {'wall +1/2':>9} {'visits':>6} {'median res':>10} {'mean res':>9} {'leave/def/t':>11} {'arrive/t':>8} {'end by death':>12}")
    for run, lab in zip(a.runs, labels):
        bc = json.load(open(os.path.join(run, "input.json")))["BC"][0]
        t, fr = read_with_angle(run, np.array(bc["Q"][:2]), a.tfrac); dt = t[1] - t[0]; T = len(fr) * dt
        W, s, P = wall_polyline(bc)
        tracks, info = build_tracks(fr, dt, 8.0, 4)
        res, nleave, narrive, nend, wall_occ = [], 0, 0, 0, 0
        for tr in tracks.values():
            if tr["q"] < 0:
                continue
            d, _ = project(W, s, tr["xy"])
            state = np.empty(len(d), int); cur = 1 if d[0] < a.swall else (0 if d[0] > a.sbulk else -1)
            for i, di in enumerate(d):
                if di < a.swall: cur = 1
                elif di > a.sbulk: cur = 0
                state[i] = cur
            wall_occ += np.sum(state == 1)
            # visits
            i = 0; n = len(state)
            while i < n:
                if state[i] == 1:
                    j = i
                    while j + 1 < n and state[j + 1] == 1: j += 1
                    started_inside = i > 0 and state[i - 1] == 0
                    ended_leave = j + 1 < n and state[j + 1] == 0
                    if started_inside: narrive += 1
                    if ended_leave: nleave += 1
                    if j == n - 1 and tr["n"][-1] < len(fr) - 1: nend += 1
                    if started_inside and ended_leave:
                        res.append((tr["n"][j] - tr["n"][i] + 1) * dt)
                    i = j + 1
                else:
                    i += 1
        nwall = wall_occ / len(fr)
        tot_visits_end = nleave + nend
        print(f"{lab:<30} {bc.get('customAnch'):>3} {nwall:9.2f} {len(res):6d} {np.median(res) if res else np.nan:10.0f} {np.mean(res) if res else np.nan:9.0f} "
              f"{nleave / max(nwall, 1e-9) / T:11.5f} {narrive / T:8.4f} {nend / max(tot_visits_end, 1):12.2f}", flush=True)


if __name__ == "__main__":
    main()
