#!/usr/bin/env python3
"""
Defect tracks with gap closing, for defects.dat (position, charge, orientation per frame).

1. Frame-to-frame linking: same-sign defects, optimal assignment (Hungarian), max jump --rlink cells.
2. Gap closing: a segment ending at frame n is joined to a same-sign segment starting at n + g
   (2 <= g <= --gap + 1, i.e. up to --gap missing frames) if the distance is below
   rlink + vgap * (g - 1) * dt; joins are chosen by a global optimal assignment (cost = distance,
   plus a small penalty per missing frame).
3. Quality: real creation/annihilation events are +-1/2 pairs. A track start (end) away from the first
   (last) frame is "paired" if an opposite-sign track starts (ends) within --rpair cells and --tpair frames.
   The paired fraction measures how much tracker flicker remains.

Library use:  tracks, info = build_tracks(frames, dt, ...)   with frames = list of (n,4) [x, y, q, angle]
  tracks: dict id -> dict(n=frame indices, xy=(m,2), q=charge sign, ang=(m,))   (frames missing inside a
  closed gap are absent from n; positions are not interpolated)
CLI:  defect_tracks.py RUN [RUN ...] [--gap 4] [--vgap 1.5]   prints track statistics before/after gap closing.
"""
import argparse
import json
import os
import sys

import numpy as np
from scipy.optimize import linear_sum_assignment

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from edge_criteria import read_with_angle  # noqa: E402


def link_frames(frames, rlink):
    """-> list of segments: dict(q, n=[...], idx=[...]) using frame-to-frame optimal assignment per sign."""
    segs = []
    open_ = {1: {}, -1: {}}            # sign -> {defect index in previous frame: segment id}
    for n, f in enumerate(frames):
        new_open = {1: {}, -1: {}}
        for sg in (1, -1):
            cur = np.where(np.sign(f[:, 2]) == sg)[0]
            prev = open_[sg]
            matched = set()
            if len(prev) and len(cur):
                pi = np.array(list(prev.keys()))
                P = frames[n - 1][pi, :2]
                D = np.hypot(*(P[:, None, :] - f[cur, None, :2].transpose(1, 0, 2)).transpose(2, 0, 1))
                r, c = linear_sum_assignment(np.where(D < rlink, D, 1e6))
                for i, j in zip(r, c):
                    if D[i, j] < rlink:
                        sid = prev[pi[i]]
                        segs[sid]["n"].append(n); segs[sid]["idx"].append(cur[j])
                        new_open[sg][cur[j]] = sid
                        matched.add(cur[j])
            for j in cur:
                if j not in matched:
                    segs.append(dict(q=sg, n=[n], idx=[j]))
                    new_open[sg][j] = len(segs) - 1
        open_ = new_open
    return segs


def close_gaps(segs, frames, dt, rlink, gap, vgap, penalty=0.5):
    """Join segment ends to later segment starts across up to `gap` missing frames."""
    nfr = len(frames)
    pos = lambda s, which: frames[s["n"][which]][s["idx"][which], :2]
    ends = [i for i, s in enumerate(segs) if s["n"][-1] < nfr - 1]
    starts = [i for i, s in enumerate(segs) if s["n"][0] > 0]
    by_start = {}
    for j in starts:
        by_start.setdefault(segs[j]["n"][0], []).append(j)
    cand = []
    for i in ends:
        ne = segs[i]["n"][-1]; pe = pos(segs[i], -1)
        for g in range(2, gap + 2):
            for j in by_start.get(ne + g, []):
                if segs[j]["q"] != segs[i]["q"] or i == j:
                    continue
                d = np.hypot(*(pos(segs[j], 0) - pe))
                if d < rlink + vgap * (g - 1) * dt:
                    cand.append((i, j, d + penalty * (g - 1)))
    if not cand:
        return segs, 0
    I = sorted({c[0] for c in cand}); J = sorted({c[1] for c in cand})
    ii = {v: k for k, v in enumerate(I)}; jj = {v: k for k, v in enumerate(J)}
    C = np.full((len(I), len(J)), 1e6)
    for i, j, cst in cand:
        C[ii[i], jj[j]] = min(C[ii[i], jj[j]], cst)
    r, c = linear_sum_assignment(C)
    nxt = {}
    for a, b in zip(r, c):
        if C[a, b] < 1e5:
            nxt[I[a]] = J[b]
    has_prev = set(nxt.values())
    merged = []
    for i in range(len(segs)):
        if i in has_prev:
            continue
        s = dict(q=segs[i]["q"], n=list(segs[i]["n"]), idx=list(segs[i]["idx"]))
        k = i
        while k in nxt:
            k = nxt[k]
            s["n"] += segs[k]["n"]; s["idx"] += segs[k]["idx"]
        merged.append(s)
    return merged, len(nxt)


def to_tracks(segs, frames):
    tracks = {}
    for t, s in enumerate(segs):
        n = np.array(s["n"]); idx = np.array(s["idx"])
        arr = np.array([frames[a][b] for a, b in zip(n, idx)])
        tracks[t] = dict(n=n, xy=arr[:, :2], q=s["q"], ang=arr[:, 3])
    return tracks


def pairing_quality(tracks, nfr, rpair, tpair):
    """Fraction of track starts/ends (away from the run edges) with an opposite-sign partner nearby."""
    st = [(tr["n"][0], tr["xy"][0], tr["q"]) for tr in tracks.values() if tr["n"][0] > tpair]
    en = [(tr["n"][-1], tr["xy"][-1], tr["q"]) for tr in tracks.values() if tr["n"][-1] < nfr - 1 - tpair]

    def frac(ev):
        if not ev:
            return np.nan, 0
        n = np.array([e[0] for e in ev]); xy = np.array([e[1] for e in ev]); q = np.array([e[2] for e in ev])
        ok = 0
        order = np.argsort(n)
        n, xy, q = n[order], xy[order], q[order]
        for a in range(len(n)):
            lo, hi = np.searchsorted(n, n[a] - tpair), np.searchsorted(n, n[a] + tpair, side="right")
            m = np.arange(lo, hi)
            m = m[(q[m] != q[a])]
            if len(m) and np.min(np.hypot(*(xy[m] - xy[a]).T)) < rpair:
                ok += 1
        return ok / len(n), len(n)
    return frac(st), frac(en)


def build_tracks(frames, dt, rlink=3.0, gap=4, vgap=1.5, rpair=4.0, tpair=2):
    segs = link_frames(frames, rlink)
    n0 = len(segs)
    if gap > 0:
        segs, njoin = close_gaps(segs, frames, dt, rlink, gap, vgap)
    else:
        njoin = 0
    tracks = to_tracks(segs, frames)
    (ps, ns), (pe, ne) = pairing_quality(tracks, len(frames), rpair, tpair)
    life = np.array([(tr["n"][-1] - tr["n"][0] + 1) * dt for tr in tracks.values()])
    return tracks, dict(segments=n0, joins=njoin, tracks=len(tracks), paired_starts=ps, n_starts=ns,
                        paired_ends=pe, n_ends=ne, life_median=float(np.median(life)),
                        life_p90=float(np.percentile(life, 90)), life_max=float(life.max()))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--gap", type=int, default=4)
    ap.add_argument("--vgap", type=float, default=1.5)
    ap.add_argument("--rlink", type=float, default=3.0)
    ap.add_argument("--tfrac", type=float, default=0.0)
    a = ap.parse_args()
    for run in a.runs:
        bc = json.load(open(os.path.join(run, "input.json")))["BC"][0]
        t, fr = read_with_angle(run, np.array(bc["Q"][:2]), a.tfrac)
        dt = t[1] - t[0]
        for g in (0, a.gap):
            _, info = build_tracks(fr, dt, a.rlink, g, a.vgap)
            print(f"{run:<36} gap {g}: segments {info['segments']:6d} joins {info['joins']:6d} tracks {info['tracks']:6d}  "
                  f"lifetime median {info['life_median']:6.0f} p90 {info['life_p90']:6.0f} max {info['life_max']:6.0f}  "
                  f"paired starts {info['paired_starts']:.2f} ({info['n_starts']})  paired ends {info['paired_ends']:.2f} ({info['n_ends']})")


if __name__ == "__main__":
    main()
