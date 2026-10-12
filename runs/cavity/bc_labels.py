#!/usr/bin/env python3
"""
Follow the wall ("BC") +1/2 defects in time, circles and Cassini ovals.

Tracks: defect_tracks.build_tracks (rlink 8, gap 4). Birth pairs: each track start is paired with the nearest
opposite-sign track start within --rpair cells and --tpair frames (greedy, one-to-one); same for deaths.
Labels at reference times t0 (every --every time units after the first --tfrac of the run):
  circle (a = 0): +1/2 on the pinned sites (edge_sites.classify; position criterion)
  oval:           slow (displacement < 5 cells over 20 time units around t0) and tilted, |orientation - wall
                  normal| within this run's [P10, P90] of slow defects within 3 cells of the wall, and s < 6
Following a label: along its track; when the track ends in a paired annihilation with a -1/2, the label passes
to the +1/2 born together with that -1/2 (birth pair), if that track is alive at that frame; otherwise the
label is lost (tracker loss, unpaired event, or partner already gone).
Output per run, vs lag tau: fraction of labels alive, fraction alive and within 3 cells of the wall, median
wall distance of alive labels, fraction that left (s > 6) at least once and came back; positions at a few lags.

Usage: bc_labels.py RUN [RUN ...] [--labels ...] [--every 100] [--tmax 1000] [--out PREFIX]
"""
import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from defect_tracks import build_tracks  # noqa: E402
from edge_criteria import read_with_angle, wall  # noqa: E402
from edge_sites import classify  # noqa: E402
from defect_spatial import load  # noqa: E402


def event_pairs(tracks, which, rpair, tpair, nfr):
    """Pair track starts (which='start') or ends ('end') of opposite sign: dict track -> partner track."""
    ev = []
    for k, tr in tracks.items():
        n = tr["n"][0] if which == "start" else tr["n"][-1]
        if (which == "start" and n == 0) or (which == "end" and n == nfr - 1):
            continue
        xy = tr["xy"][0] if which == "start" else tr["xy"][-1]
        ev.append((n, xy, tr["q"], k))
    ev.sort(key=lambda e: e[0])
    plus = [e for e in ev if e[2] > 0]; minus = [e for e in ev if e[2] < 0]
    mt = np.array([m[0] for m in minus]); mxy = np.array([m[1] for m in minus]).reshape(-1, 2)
    cand = []
    for p in plus:
        lo, hi = np.searchsorted(mt, p[0] - tpair), np.searchsorted(mt, p[0] + tpair, side="right")
        if hi <= lo:
            continue
        d = np.hypot(*(mxy[lo:hi] - p[1]).T)
        for o in np.where(d < rpair)[0]:
            m = minus[lo + o]
            cand.append((d[o] + 0.5 * abs(m[0] - p[0]), p[3], m[3]))
    cand.sort()
    used, pair = set(), {}
    for c, a, b in cand:
        if a in used or b in used:
            continue
        used.add(a); used.add(b); pair[a] = b; pair[b] = a
    return pair


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--labels", nargs="+")
    ap.add_argument("--every", type=float, default=100.0)
    ap.add_argument("--tmax", type=float, default=1000.0)
    ap.add_argument("--tfrac", type=float, default=0.2)
    ap.add_argument("--rpair", type=float, default=6.0)
    ap.add_argument("--tpair", type=int, default=5)
    ap.add_argument("--out", default="bc_labels")
    a = ap.parse_args()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    labels = a.labels or [os.path.basename(os.path.normpath(r)) for r in a.runs]
    results = []
    for run, lab in zip(a.runs, labels):
        bc = json.load(open(os.path.join(run, "input.json")))["BC"][0]
        is_oval = bc.get("shape", 0) == 1
        t, fr = read_with_angle(run, np.array(bc["Q"][:2]), 0.0)
        dt = t[1] - t[0]; nfr = len(fr)
        tracks, info = build_tracks(fr, dt, 8.0, 4)
        births = event_pairs(tracks, "start", a.rpair, a.tpair, nfr)
        deaths = event_pairs(tracks, "end", a.rpair, a.tpair, nfr)
        # lookup: frame -> list of (track, row)
        at = {}
        for k, tr in tracks.items():
            for r_, n in enumerate(tr["n"]):
                at.setdefault(n, []).append((k, r_))
        h = int(round(10 / dt))                       # half window of the 20-time-unit mobility test
        if is_oval:
            sl, tl = [], []
            for n in range(int(a.tfrac * nfr) + h, nfr - h, max(1, int(round(10 / dt)))):
                for k, r_ in at.get(n, []):
                    tr = tracks[k]
                    if tr["q"] < 0 or r_ - h < 0 or r_ + h >= len(tr["n"]):
                        continue
                    s, na = wall(bc, tr["xy"][r_, 0:1], tr["xy"][r_, 1:2])
                    D = np.hypot(*(tr["xy"][r_ + h] - tr["xy"][r_ - h]))
                    if D < 5 and s[0] < 3:
                        tl.append(abs(np.angle(np.exp(1j * (tr["ang"][r_] - na[0])))))
            win = (np.percentile(tl, 10), np.percentile(tl, 90))
        else:
            d = load(run, a.tfrac)
            c = classify(d, 1.0, 0.02, 0.0, 4.2, 0.4)
            sites, dsite = c["sites"] + 0.0, c["dsite"]       # sites relative to the centre (load() subtracts Q)
            win = None
        lags = np.arange(0, a.tmax + dt, dt * max(1, int(round(10 / dt))))
        L = len(lags)
        alive = np.zeros(L); atwall = np.zeros(L); came_back = np.zeros(L); s_med = [[] for _ in range(L)]
        snap = {0: [], int(L * 0.1): [], int(L * 0.5): [], L - 1: []}
        nlab = 0; nref = 0; ncycle = [0]
        step = max(1, int(round(a.every / dt)))
        for n0 in range(int(a.tfrac * nfr) + h, nfr - int(a.tmax / dt) - 1, step):
            nref += 1
            labs = []
            for k, r_ in at.get(n0, []):
                tr = tracks[k]
                if tr["q"] < 0:
                    continue
                x, y = tr["xy"][r_]
                s, na = wall(bc, np.array([x]), np.array([y]))
                if is_oval:
                    if r_ - h < 0 or r_ + h >= len(tr["n"]) or s[0] > 6:
                        continue
                    D = np.hypot(*(tr["xy"][r_ + h] - tr["xy"][r_ - h]))
                    tilt = abs(np.angle(np.exp(1j * (tr["ang"][r_] - na[0]))))
                    if D < 5 and win[0] <= tilt <= win[1]:
                        labs.append(k)
                else:
                    if len(sites) and np.min(np.hypot(*(sites - [x, y]).T)) < dsite:
                        labs.append(k)
            nlab += len(labs)
            for k in labs:
                cur = k; left = False; seen = {k}
                for li, lag in enumerate(lags):
                    n = n0 + int(round(lag / dt))
                    # pass the label on while the current track has ended before n
                    while cur is not None and tracks[cur]["n"][-1] < n:
                        m = deaths.get(cur)
                        if m is None or tracks[m]["q"] > 0:
                            cur = None; break
                        p = births.get(m)
                        ne = tracks[cur]["n"][-1]
                        if p is None or tracks[p]["q"] < 0 or tracks[p]["n"][0] > ne or tracks[p]["n"][-1] < ne:
                            cur = None; break
                        if p in seen:              # cycle in the partner chain: stop
                            ncycle[0] += 1; cur = None; break
                        seen.add(p); cur = p
                    if cur is None:
                        break
                    tr = tracks[cur]
                    r_ = np.searchsorted(tr["n"], n)
                    if r_ >= len(tr["n"]):
                        r_ = len(tr["n"]) - 1
                    x, y = tr["xy"][r_]
                    s = wall(bc, np.array([x]), np.array([y]))[0][0]
                    alive[li] += 1
                    if s < 3:
                        atwall[li] += 1
                        if left:
                            came_back[li] += 1
                    if s > 6:
                        left = True
                    s_med[li].append(s)
                    if li in snap:
                        snap[li].append((x, y))
        results.append(dict(lab=lab, bc=bc, lags=lags, alive=alive / max(nlab, 1), atwall=atwall / max(nlab, 1),
                            came_back=came_back / max(nlab, 1), smed=np.array([np.median(v) if v else np.nan for v in s_med]),
                            snap=snap, nlab=nlab, nref=nref, win=win, info=info, ncycle=ncycle[0]))
        r = results[-1]
        def at_lag(arr, T):
            return arr[min(np.searchsorted(r["lags"], T), len(r["lags"]) - 1)]
        print(f"{lab:<26} labels/ref {nlab / max(nref, 1):5.1f} ({nref} refs)  tilt window {('%.2f-%.2f' % win) if win else 'site rule'}  "
              f"tracks {info['tracks']}  | at wall after 50/200/1000: {at_lag(r['atwall'], 50):.2f} {at_lag(r['atwall'], 200):.2f} {at_lag(r['atwall'], 1000):.2f}  "
              f"cycles {r['ncycle']}  alive {at_lag(r['alive'], 1000):.2f}  median s alive @1000 {at_lag(r['smed'], 1000):.1f}  came back by 1000 {at_lag(r['came_back'], 1000):.2f}", flush=True)
    # figure: runs with the same label are averaged (seeds)
    groups = {}
    for r in results:
        groups.setdefault(r["lab"], []).append(r)
    fig, axs = plt.subplots(1, 2, figsize=(10.5, 3.9))
    blues = ["#86b6ef", "#2a78d6", "#104281"]; oranges = ["#f4a27f", "#eb6834", "#c24d1c", "#8f3413", "#5c210c"]
    ib = io = 0
    for lab, rr in groups.items():
        oval = rr[0]["bc"].get("shape", 0) == 1
        col = oranges[io % len(oranges)] if oval else blues[ib % len(blues)]
        io += oval; ib += (not oval)
        L = min(len(r["lags"]) for r in rr); lags = rr[0]["lags"][:L]
        aw = np.mean([r["atwall"][:L] for r in rr], 0); al = np.mean([r["alive"][:L] for r in rr], 0)
        sm = np.nanmean([r["smed"][:L] for r in rr], 0)
        axs[0].plot(lags, aw, color=col, lw=2.0, label=f"{lab} ({len(rr)} seed{'s' if len(rr) > 1 else ''})")
        axs[0].plot(lags, al, color=col, lw=0.9, ls=":")
        axs[1].plot(lags, sm, color=col, lw=1.8)
    axs[0].set_xlabel("lag τ (time units)"); axs[0].set_ylabel("labelled wall +½ still within 3 cells of the wall\n(dotted: label still tracked)")
    axs[1].set_xlabel("lag τ (time units)"); axs[1].set_ylabel("median wall distance of tracked labels")
    axs[0].set_ylim(0, 1.02)
    for ax in axs:
        ax.grid(color="#e9ebee", lw=0.6); ax.spines[["top", "right"]].set_visible(False)
    axs[0].legend(frameon=False, fontsize=7)
    fig.tight_layout(); fig.savefig(a.out + ".png", dpi=150)
    print("written", a.out + ".png")


if __name__ == "__main__":
    main()
