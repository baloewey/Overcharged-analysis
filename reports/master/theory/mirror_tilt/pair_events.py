#!/usr/bin/env python3
"""Pair-creation-cost test, event side (user, 2026-10-19). Mirror shape sweep (wall -1/2), second half, every defect frame (2 t.u.).
Tracking: -1/2 and +1/2 linked frame to frame separately (mutual nearest neighbour, < 3 cells). Wall band: depth < 4 cells.
Per -1/2 track, wall stays (consecutive frames in the band, >= 5 frames = 10 t.u.; shorter = flicker, ignored):
  entry: 'arrive' (track was deeper than 6 cells before the stay), 'excursion_in' (track came from 4-6 only), or 'born' (track starts in the band);
  exit:  'depart' (track later goes deeper than 6), 'excursion_out' (stays within 6) or 'die' (track ends in the band; 'annih' if a +1/2 track also ends within 3 cells in the same/next frame).
Smoothed excess n(t) = N_wall(t) - 2|k| (running median, 11 frames); time spent at each n; entries and exits counted against n at the event.
Writes pair_events.pkl (per run: k, R_eq, rho, family, aspect, T, occupancy, entries/exits by type and n)."""
import csv, glob, json, os, sys, pickle
import numpy as np
from multiprocessing import Pool
from scipy.spatial import cKDTree
from scipy.ndimage import median_filter
ROOT = "/home/baloewe/Physics/Overcharged/runs"; sys.path.insert(0, ROOT + "/cavity")
from edge_criteria import read_with_angle  # noqa
from wall_fourier import wall_polyline  # noqa
SW = ROOT + "/shapes_mirror"; H = os.path.dirname(os.path.abspath(__file__)); NMIN = 5; NS = range(-3, 5)
def link(prev, cur):
    """mutual nearest neighbours within 3 cells -> dict cur_index -> prev_index"""
    if not len(prev) or not len(cur): return {}
    d1, j1 = cKDTree(prev).query(cur); d2, j2 = cKDTree(cur).query(prev)
    return {i: j1[i] for i in range(len(cur)) if d1[i] < 3 and j2[j1[i]] == i}
def tracks(frames):
    """frames: list of (N,2) arrays -> list of tracks [(frame, index)]"""
    tid_prev = []; T = []
    for f, p in enumerate(frames):
        lk = link(frames[f - 1], p) if f else {}; tid = []
        for i in range(len(p)):
            if i in lk: t = tid_prev[lk[i]]; T[t].append((f, i))
            else: t = len(T); T.append([(f, i)])
            tid.append(t)
        tid_prev = tid
    return T
def one(run):
    try:
        bc = json.load(open(run + "/input.json"))["BC"][0]; k = bc["customAnch"]; Q = bc["Q"]
        W, s, P = wall_polyline(bc); tree = cKDTree(W)
        t, fr = read_with_angle(run, np.array(Q[:2], float), 0.5); dt = float(np.median(np.diff(t)))
        M = [f[f[:, 2] < 0, :2] for f in fr]; Pp = [f[f[:, 2] > 0, :2] for f in fr]
        DM = [tree.query(m)[0] if len(m) else np.zeros(0) for m in M]
        N = np.array([np.sum(d < 4) for d in DM]); n = np.round(median_filter(N.astype(float), 11, mode="nearest")).astype(int) - 2 * abs(k)
        occ = {v: float(np.sum(n == v)) * dt for v in NS}
        TM = tracks(M); TP = tracks(Pp); pend = {}
        for tr in TP:
            f, i = tr[-1]; pend.setdefault(f, []).append(Pp[f][i])
        ev = {(typ, v): 0 for typ in ("arrive", "born", "excursion_in", "depart", "die", "annih", "excursion_out") for v in NS}
        nf = len(M); stays = []
        for tr in TM:
            dep = np.array([DM[f][i] for f, i in tr]); inb = dep < 4; fs = np.array([f for f, _ in tr])
            # runs of consecutive in-band frames
            j = 0
            while j < len(tr):
                if not inb[j]: j += 1; continue
                j0 = j
                while j < len(tr) and inb[j]: j += 1
                if j - j0 < NMIN: continue
                fa, fb = fs[j0], fs[j - 1]
                stays.append((fb - fa + 1) * dt)
                if fa > 5:
                    typ = ("arrive" if (dep[:j0] > 6).any() else "excursion_in") if j0 > 0 else "born"; v = int(n[fa])
                    if v in NS: ev[(typ, v)] += 1
                if fb < nf - 6:
                    v = int(n[fb])
                    if j < len(tr): typ = "depart" if (dep[j:] > 6).any() else "excursion_out"
                    else:
                        z = M[fb][tr[j - 1][1]]; near = [p for ff in (fb, fb + 1) for p in pend.get(ff, [])]
                        typ = "annih" if near and np.min(np.hypot(*(np.array(near) - z).T)) < 3 else "die"
                    if v in NS: ev[(typ, v)] += 1
        return dict(run=os.path.relpath(run, SW), k=k, T=nf * dt, P=P, occ=occ, ev=ev, Nmean=N.mean(), stays=stays)
    except Exception as e:
        return dict(run=os.path.relpath(run, SW), error=repr(e))
if __name__ == "__main__":
    runs = sorted(os.path.dirname(f) for f in glob.glob(SW + "/*/*/*/*/s*/DONE"))
    if len(sys.argv) > 1: runs = [r for r in runs if sys.argv[1] in r]
    man = {r["dir"]: r for r in csv.DictReader(open(SW + "/manifest_shape_mirror.tsv"), delimiter="\t")}
    with Pool(6) as p: res = p.map(one, runs, chunksize=2)
    for o in res: o.update({c: man.get(o["run"], {}).get(c) for c in ("family", "aspect", "Req_la", "rho")})
    pickle.dump(res, open(H + ("/pair_events_test.pkl" if len(sys.argv) > 1 else "/pair_events.pkl"), "wb"))
    print(len(res), "runs, errors", [o for o in res if "error" in o][:3])
