#!/usr/bin/env python3
"""
Shielding test (user's argument, 2026-10-05): is the attraction between a mobile bulk defect and the opposite-sign wall defects screened?

Gap-closed tracks (build_tracks, rlink 8 cells, gap 4 frames), frames after 20%. Wall-sign = sign of the winding (-1/2 for -k, +1/2 for +k);
the "mobile" defect is the opposite sign (+1/2 for -k, -1/2 for +k).
Encounter: the first frame in which a mobile track comes within RE = 4 cells of a given opposite-sign track (one encounter per pair of
tracks). The target is classified by its distance s to the wall at that frame: wall (s < 4 cells) or bulk (s > 6 cells); others skipped.
Outcome: annihilation if both tracks end within TP = 5 frames (10 time units) of each other, at most 4 cells apart, and before the end of
the run; otherwise the pair separated (escape). Reported per class: encounters, annihilation probability, and the same restricted to
encounters at a minimum distance <= 2 cells (closest approach over the next 10 frames).
Usage (from runs/): shielding_test.py [--workers 10] [--set neg|pos|both] [--rhos ...] [--rmin 1.0 --rmax 2.3] [--runs DIR ...]
Writes ../reports/master/data/shielding_<set>.tsv (per run) and prints pooled numbers.
"""
import argparse, glob, json, os, sys
from multiprocessing import Pool
import numpy as np
from scipy.spatial import cKDTree
sys.path.insert(0, "cavity")
from edge_criteria import read_with_angle      # noqa: E402
from wall_fourier import wall_polyline, project  # noqa: E402
from defect_tracks import build_tracks         # noqa: E402

RE, TP, RP, SW, SB = 4.0, 5, 4.0, 4.0, 6.0
LA = 22.21
DATA = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data"))


def one(run):
    try:
        bc = json.load(open(os.path.join(run, "input.json")))["BC"][0]
        k = int(bc["customAnch"]); wsign = 1 if k > 0 else -1
        W, s, P = wall_polyline(bc)
        t, fr = read_with_angle(run, np.array(bc["Q"][:2], float), 0.2); dt = t[1] - t[0]; nfr = len(fr)
        tracks, _ = build_tracks(fr, dt, 8.0, 4)
        # per-frame index of track positions
        per = [[] for _ in range(nfr)]
        end = {}
        for tid, T in tracks.items():
            d, _ = project(W, s, T["xy"])
            for n, xy, dd in zip(T["n"], T["xy"], d):
                per[n].append((tid, T["q"], xy[0], xy[1], dd))
            end[tid] = (T["n"][-1], T["xy"][-1], T["q"])
        seen = set(); stats = {"wall": [0, 0, 0, 0], "bulk": [0, 0, 0, 0]}  # encounters, annihilations, close encounters, close annihilations
        for n in range(nfr):
            if not per[n]: continue
            arr = per[n]
            mob = [a for a in arr if np.sign(a[1]) == -wsign]; tgt = [a for a in arr if np.sign(a[1]) == wsign]
            if not mob or not tgt: continue
            tree = cKDTree([(a[2], a[3]) for a in tgt])
            for m in mob:
                for j in tree.query_ball_point((m[2], m[3]), RE):
                    g = tgt[j]; key = (m[0], g[0])
                    if key in seen: continue
                    seen.add(key)
                    cls = "wall" if g[4] < SW else ("bulk" if g[4] > SB else None)
                    if cls is None: continue
                    # closest approach over the next 10 frames
                    Tm, Tg = tracks[m[0]], tracks[g[0]]
                    im = {v: i for i, v in enumerate(Tm["n"])}; ig = {v: i for i, v in enumerate(Tg["n"])}
                    dmin = min((np.hypot(*(Tm["xy"][im[v]] - Tg["xy"][ig[v]])) for v in range(n, n + 11) if v in im and v in ig), default=RE)
                    em, eg = end[m[0]], end[g[0]]
                    ann = (em[0] < nfr - 1 and eg[0] < nfr - 1 and abs(em[0] - eg[0]) <= TP and em[0] >= n
                           and np.hypot(*(em[1] - eg[1])) <= RP)
                    st = stats[cls]; st[0] += 1; st[1] += ann
                    if dmin <= 2.0: st[2] += 1; st[3] += ann
        return dict(run=run, k=k, R=bc["R"], P=P, **{f"{c}_{i}": stats[c][j] for c in stats for j, i in enumerate(("enc", "ann", "cenc", "cann"))})
    except Exception as ex:
        import traceback
        return dict(run=run, err=traceback.format_exc()[-300:])


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=10); ap.add_argument("--set", default="both")
    ap.add_argument("--rhos", nargs="+", default=["0.05", "0.06", "0.08", "0.10", "0.12"])
    ap.add_argument("--rmin", type=float, default=1.0); ap.add_argument("--rmax", type=float, default=2.3)
    ap.add_argument("--runs", nargs="+"); ap.add_argument("--tag", default=None)
    a = ap.parse_args()
    if a.runs:
        runs = a.runs; tag = a.tag or "custom"
    else:
        sets = {"neg": ["main_act0.07_neg"], "pos": ["main_act0.07"], "both": ["main_act0.07_neg", "main_act0.07"]}[a.set]
        runs = []
        for st in sets:
            for rho in a.rhos:
                for d in sorted(glob.glob(f"cavity/server_main/{st}/rho{rho}/k*/s0*")):
                    bc = json.load(open(d + "/input.json"))["BC"][0]
                    if a.rmin <= bc["R"] / LA <= a.rmax and abs(bc["customAnch"]) >= 2: runs.append(d)
        tag = a.tag or a.set
    print(len(runs), "runs", flush=True)
    with Pool(a.workers) as p: res = p.map(one, runs)
    bad = [r for r in res if "err" in r]
    for b in bad: print("FAILED", b["run"], b["err"])
    res = [r for r in res if "err" not in r]
    os.makedirs(DATA, exist_ok=True)
    cols = ["run", "k", "R", "P", "wall_enc", "wall_ann", "wall_cenc", "wall_cann", "bulk_enc", "bulk_ann", "bulk_cenc", "bulk_cann"]
    with open(os.path.join(DATA, f"shielding_{tag}.tsv"), "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in res: f.write("\t".join(str(r[c]) for c in cols) + "\n")
    for sign in (-1, 1):
        rs = [r for r in res if np.sign(r["k"]) == sign]
        if not rs: continue
        S = {c: sum(r[c] for r in rs) for c in cols[4:]}
        mob = "+1/2" if sign < 0 else "-1/2"; wl = "-1/2" if sign < 0 else "+1/2"
        print(f"\nwinding {'-k' if sign < 0 else '+k'} ({len(rs)} runs): mobile {mob} meeting a {wl}")
        for c in ("wall", "bulk"):
            e, an, ce, ca = S[f"{c}_enc"], S[f"{c}_ann"], S[f"{c}_cenc"], S[f"{c}_cann"]
            pe = an / e if e else np.nan; pc = ca / ce if ce else np.nan
            se = np.sqrt(pe * (1 - pe) / e) if e else np.nan; sc = np.sqrt(pc * (1 - pc) / ce) if ce else np.nan
            print(f"  target in the {c:4s}: {e:7d} encounters, P(annihilation) = {pe:.3f} ± {se:.3f};  closest approach <= 2 cells: {ce:6d}, P = {pc:.3f} ± {sc:.3f}")
        # per-run spread
        pw = [r["wall_ann"] / r["wall_enc"] for r in rs if r["wall_enc"] >= 20]; pb = [r["bulk_ann"] / r["bulk_enc"] for r in rs if r["bulk_enc"] >= 20]
        if pw and pb:
            print(f"  per run (>= 20 encounters): wall median {np.median(pw):.3f} (IQR {np.percentile(pw, 25):.3f}-{np.percentile(pw, 75):.3f}), "
                  f"bulk median {np.median(pb):.3f} (IQR {np.percentile(pb, 25):.3f}-{np.percentile(pb, 75):.3f})")


if __name__ == "__main__":
    main()
