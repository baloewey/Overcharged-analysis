#!/usr/bin/env python3
"""Per-lobe analysis of Cassini ovals (user, 2026-10-12): do the lobes of pinched ovals carry their own vortex and their own bulk charge, and do defects transfer between lobes?
Shape sweep, second half. Ovals b/a 1.10, 1.30 (pinched, concave for b/a < sqrt 2) and 1.45, 1.60 (convex references); ellipses of comparable aspect as controls for the transfer rate.
Lobes: x < 0 (left) and x > 0 (right) in cavity-centred coordinates (the neck is at x = 0).
Per run: (1) bulk charge per lobe Q_L, Q_R (defects deeper than 4 cells), wall +1/2 per lobe; distribution of (Q_L, Q_R): P(Q_L = Q_R), corr(Q_L, Q_R);
(2) lobe circulation from the flow frames (every 100 t.u.): mean azimuthal velocity about each lobe's centroid within the lobe, signs -> co-rotating fraction, reversals; whole-cavity circulation;
(3) transfers: bulk defects linked frame to frame (same sign, < 2.5 cells) that change the sign of x while deeper than 4 cells -> crossings per bulk defect per 1000 t.u., by sign;
    wall +1/2 crossing x = 0 along the wall (< 4 cells) likewise.
Usage: lobes.py [NPROC] -> lobes.npy, lobes.txt"""
import os, sys, csv, gzip, numpy as np
from multiprocessing import Pool
from scipy.spatial import cKDTree
from matplotlib.path import Path
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
import energetics as En  # noqa
from selfflow import load_flow  # noqa
SW = "/home/baloewe/Physics/Overcharged/runs/shapes_server"
def one(rel):
    try:
        run = os.path.join(SW, rel); import json; bc = json.load(open(run + "/input.json"))["BC"][0]; Q = bc["Q"]; k = bc["customAnch"]; sg = 1 if k > 0 else -1
        W, s, P = En.wall_polyline(bc); tree = cKDTree(W); t, fr = En.read_with_angle(run, np.array(Q[:2], float), 0.5); t = np.array(t)
        QL, QR, WL, WR = [], [], [], []; cross_b = {1: 0, -1: 0}; cross_w = 0; nb_samples = 0
        prev = None
        for n, f in enumerate(fr):
            if not len(f): prev = None; continue
            d = tree.query(f[:, :2])[0]; b = d >= 4; w = (np.sign(f[:, 2]) == sg) & (d < 4)
            QL.append(0.5 * np.sign(f[b & (f[:, 0] < 0), 2]).sum()); QR.append(0.5 * np.sign(f[b & (f[:, 0] >= 0), 2]).sum())
            WL.append((w & (f[:, 0] < 0)).sum()); WR.append((w & (f[:, 0] >= 0)).sum()); nb_samples += b.sum()
            if prev is not None:
                pf, pd = prev
                for i in range(len(f)):
                    same = np.flatnonzero(np.sign(pf[:, 2]) == np.sign(f[i, 2]))
                    if not len(same): continue
                    dd = np.hypot(*(pf[same, :2] - f[i, :2]).T); kk = same[np.argmin(dd)]
                    if dd.min() > 2.5 or np.sign(pf[kk, 0]) == np.sign(f[i, 0]): continue
                    if d[i] >= 4 and pd[kk] >= 4: cross_b[int(np.sign(f[i, 2]))] += 1
                    elif d[i] < 4 and pd[kk] < 4 and np.sign(f[i, 2]) == sg: cross_w += 1
            prev = (f, d)
        dt = float(np.median(np.diff(t))); T = len(fr) * dt
        # lobe circulation
        rows, ft, idx = load_flow(run); X = rows[:, 1] + 0.5 - Q[0]; Y = rows[:, 2] + 0.5 - Q[1]; ins = Path(W).contains_points(np.c_[X, Y])
        L_ = ins & (X < 0); R_ = ins & (X >= 0); cL = np.array([X[L_].mean(), Y[L_].mean()]); cR = np.array([X[R_].mean(), Y[R_].mean()])
        circ = []
        for tt in ft[len(ft) // 2:]:
            ii = idx[tt]; x, y, vx, vy = X[ii], Y[ii], rows[ii, 4], rows[ii, 5]; m = ins[ii]; out = []
            for c, sel in ((cL, m & (x < 0)), (cR, m & (x >= 0)), (np.zeros(2), m)):
                rx, ry = x[sel] - c[0], y[sel] - c[1]; r = np.hypot(rx, ry) + 1e-9; out.append(np.mean((rx * vy[sel] - ry * vx[sel]) / r))
            circ.append(out)
        return dict(rel=rel, QL=np.array(QL), QR=np.array(QR), WL=np.array(WL), WR=np.array(WR), cross_b=cross_b, cross_w=cross_w, nb=nb_samples / max(len(fr), 1), T=T, circ=np.array(circ))
    except Exception as ex:
        return dict(rel=rel, err=str(ex))
if __name__ == "__main__":
    man = {r["dir"]: r for r in csv.DictReader(open(SW + "/manifest_shape_sweep.tsv"), delimiter="\t")}
    rels = [r for r in man if os.path.exists(os.path.join(SW, r, "DONE")) and ((man[r]["family"] == "oval" and man[r]["shape"] in ("1.1", "1.3", "1.45", "1.6")) or (man[r]["family"] == "ellipse" and float(man[r]["aspect"]) >= 1.95))]
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 3) as pool: res = pool.map(one, rels, chunksize=2)
    ok = [r for r in res if "err" not in r]; np.save(os.path.join(H, "lobes.npy"), np.array(ok, dtype=object), allow_pickle=True)
    L = [f"{len(ok)} runs ({len(res) - len(ok)} errors). Groups: oval b/a (pinched < 1.414) and ellipses aspect >= 2 (controls). Means over runs; ranges over sizes in brackets where useful."]
    groups = {}
    for r in ok:
        m = man[r["rel"]]; key = (m["family"], m["shape"] if m["family"] == "oval" else f"asp{float(m['aspect']):.1f}", m["Req_la"]); groups.setdefault(key, []).append(r)
    L.append("family shape R_eq | runs | Q_L Q_R (bulk charge per lobe) | P(Q_L = Q_R) | corr(Q_L,Q_R) | wall L R | co-rotating fraction | |lobe circ| | |whole circ| | bulk crossings per bulk defect per 1000 t.u. (+1/2, -1/2) | wall crossings per 1000 t.u.")
    for key in sorted(groups):
        S = groups[key]
        ql = np.concatenate([r["QL"] for r in S]); qr = np.concatenate([r["QR"] for r in S]); co = np.concatenate([np.sign(r["circ"][:, 0]) == np.sign(r["circ"][:, 1]) for r in S if len(r["circ"])])
        lc = np.concatenate([np.abs(r["circ"][:, :2]).ravel() for r in S if len(r["circ"])]); wc = np.concatenate([np.abs(r["circ"][:, 2]) for r in S if len(r["circ"])])
        cbp = np.mean([1000 * r["cross_b"][1] / max(r["nb"] * r["T"], 1e-9) for r in S]); cbm = np.mean([1000 * r["cross_b"][-1] / max(r["nb"] * r["T"], 1e-9) for r in S])
        cw = np.mean([1000 * r["cross_w"] / r["T"] for r in S])
        L.append(f"{key[0]:7s} {key[1]:6s} R{key[2]} | {len(S)} | {ql.mean():+.2f} {qr.mean():+.2f} | {np.mean(ql == qr):.2f} | {np.corrcoef(ql, qr)[0, 1]:+.2f} | {np.mean([r['WL'].mean() for r in S]):.1f} {np.mean([r['WR'].mean() for r in S]):.1f} | {co.mean():.2f} | {lc.mean():.3f} | {wc.mean():.3f} | {cbp:.1f} {cbm:.1f} | {cw:.2f}")
    print("\n".join(L)); open(os.path.join(H, "lobes.txt"), "w").write("\n".join(L) + "\n")
