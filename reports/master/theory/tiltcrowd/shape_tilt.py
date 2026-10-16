#!/usr/bin/env python3
"""Tilt vs spacing in elongated cavities (user, 2026-10-15): does the circle collapse tau(a) hold in ellipses and Cassini ovals?
TILT DEFINITION: tau_j = psi_j - theta_n(s_j), from the local outward normal at the nearest wall point (as in circles). By the lock this equals the
registry: tau = 2 (theta_w - theta_n) + c with a constant c (checked per frame: 'lock' = |<exp(i c_j)>|, lock_c = its angle). Derivation of c: Arclength-uniform anchoring (eMode 2/3): director at the wall theta_w(s) = 2 pi k s/P
(lab frame; s from the +x axis, ccw). A +1/2 sitting at a homeotropic point (theta_w = theta_n mod pi) with its tail along the normal is untilted;
the lock (zero mode) rotates the tail by 2 (theta_w - theta_n) when the defect slides along the wall. So the common tilt of the layer is
    c_j = psi_j - 2 theta_w(s_j) + theta_n(s_j)        (psi = defect angle, theta_n = outward normal angle at the nearest wall point)
In a circle (theta_w = k phi, theta_n = phi) c = psi - (2k-1) phi, the constant of the lock tau = (2k-2) phi + c, equal to the tilt from the normal at a site.
In an ellipse the homeotropic sites are unevenly spaced (theta_n does not advance uniformly), so the tilt from the local normal mixes in the registry;
c removes it. LOCAL SITE SPACING a(s) = pi / (2 pi k/P - kappa(s)) (kappa = wall curvature; circle: 2 pi R/(2k - 2)); also the actual arc distance to the
nearest wall neighbour.
Per frame (every 2 t.u., second half), wall +1/2 (within 4 cells): c_j; local coherence of each defect with its two wall neighbours
|<exp(i c)>|_{j-1,j,j+1}; LOCALLY COHERENT triplets (> 0.9) give |arg<exp(i c)>| = the local tilt. Also the whole-layer coherence |<exp(i c)>| (lock check).
Runs: runs/shapes_server (675: circle, ellipse eps 1.2-2.5, Cassini ovals b/a 1.1-2.5; R_eq 1/1.5/2 l_a; 3 seeds).  -> shape_tilt.npz, shape_tilt.txt"""
import os, sys, glob, json, numpy as np
from multiprocessing import Pool
from scipy.spatial import cKDTree
ROOT = "/home/baloewe/Physics/Overcharged"; H = os.path.dirname(os.path.abspath(__file__)); SW = ROOT + "/runs/shapes_server"
sys.path.insert(0, ROOT + "/runs/cavity")
from edge_criteria import read_with_angle  # noqa
from wall_fourier import wall_polyline  # noqa


def one(run):
    J = json.load(open(run + "/input.json")); bc = J["BC"][0]; k = bc["customAnch"]
    if "aInv" not in bc: bc["aInv"] = [1, 1, 1]
    W, s, P = wall_polyline(bc); n = len(W); tree = cKDTree(W)
    dW = (np.roll(W, -1, 0) - np.roll(W, 1, 0)); tn = np.arctan2(dW[:, 1], dW[:, 0]); thn = tn - np.pi / 2           # ccw wall: outward normal = tangent - 90 deg
    ds = np.hypot(*(np.roll(W, -1, 0) - W).T); kap = np.angle(np.exp(1j * (np.roll(tn, -1) - tn))) / ds                    # curvature (> 0 convex)
    asite = np.pi / np.maximum(2 * np.pi * k / P - kap, 1e-6)
    t, fr = read_with_angle(run, np.array(bc["Q"][:2], float), 0.5)
    out, Cl, Lk, Lc = [], [], [], []
    for f in fr:
        if not len(f): continue
        d, j = tree.query(f[:, :2]); w = (f[:, 2] > 0) & (d < 4)
        if w.sum() < 3: continue
        jj = j[w]; o = np.argsort(s[jj]); jj = jj[o]; psi = f[w, 3][o]
        lock = np.exp(1j * (psi - 2 * (2 * np.pi * k * s[jj] / P) + thn[jj])); Lk.append(abs(np.mean(lock))); Lc.append(np.angle(np.mean(lock)))
        c = np.angle(np.exp(1j * (psi - thn[jj]))); Cl.append(abs(np.mean(np.exp(1j * c))))       # tilt from the local outward normal
        sj = s[jj]; gap = np.diff(np.r_[sj, sj[0] + P]); nn = np.minimum(gap, np.roll(gap, 1))
        z = (np.exp(1j * c) + np.exp(1j * np.roll(c, 1)) + np.exp(1j * np.roll(c, -1))) / 3
        ok = np.abs(z) > 0.9
        for i in np.flatnonzero(ok): out.append((asite[jj[i]], nn[i], abs(np.angle(z[i])), kap[jj[i]]))
    o = np.array(out, np.float32) if out else np.zeros((0, 4), np.float32)
    return dict(run=os.path.relpath(run, SW), k=k, P=P, C=np.median(Cl) if Cl else np.nan, lock=np.median(Lk) if Lk else np.nan,
                lock_c=np.degrees(np.angle(np.mean(np.exp(1j * np.array(Lc))))) if Lc else np.nan, trip=o)


if __name__ == "__main__":
    import csv
    man = {r["dir"]: r for r in csv.DictReader(open(SW + "/manifest_shape_sweep.tsv"), delimiter="\t")}
    runs = [os.path.join(SW, r) for r in man if os.path.exists(os.path.join(SW, r, "DONE"))]
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 6) as pool: res = pool.map(one, runs, chunksize=4)
    fam = np.array([man[x["run"]]["family"] for x in res]); asp = np.array([float(man[x["run"]]["aspect"]) for x in res]); req = np.array([float(man[x["run"]]["Req_la"]) for x in res])
    T = [x["trip"] for x in res]; idx = np.concatenate([np.full(len(tt), i) for i, tt in enumerate(T)]); A = np.concatenate(T)
    np.savez(os.path.join(H, "shape_tilt.npz"), trip=A, run_idx=idx, fam=fam, asp=asp, req=req, k=np.array([x["k"] for x in res]), P=np.array([x["P"] for x in res]),
             C=np.array([x["C"] for x in res]), lock=np.array([x["lock"] for x in res]), lock_c=np.array([x["lock_c"] for x in res]), run=np.array([x["run"] for x in res]))
    L = ["SHAPE TILT (tilt from the local normal; lock |c| checked; locally coherent triplets > 0.9; second half; 3 seeds). Local tilt (deg, median) by local site spacing a(s) (cells)",
         "group | runs | median layer tilt coherence | median lock |c-resultant| (angle) | triplets | tilt for a in [6,9) [9,12) [12,16) [16,22) [22,40)"]
    bins = [(6, 9), (9, 12), (12, 16), (16, 22), (22, 40)]
    for g, m in (("circle", fam == "circle"), ("aspect < 1.3", (fam != "circle") & (asp < 1.3)), ("aspect 1.3-1.7", (asp >= 1.3) & (asp < 1.7)), ("aspect >= 1.7", asp >= 1.7),
                 ("ellipses >= 1.5", (fam == "ellipse") & (asp >= 1.5)), ("ovals >= 1.5", (fam == "oval") & (asp >= 1.5))):
        sel = np.isin(idx, np.flatnonzero(m)); a_ = A[sel]
        L.append(f"{g:16s} | {m.sum():3d} | {np.nanmedian([res[i]['C'] for i in np.flatnonzero(m)]):.2f} | {np.nanmedian([res[i]['lock'] for i in np.flatnonzero(m)]):.2f} ({np.nanmedian([res[i]['lock_c'] for i in np.flatnonzero(m)]):+.0f}) | {sel.sum():7d} | " +
                 "  ".join(f"{np.degrees(np.median(a_[(a_[:, 0] >= lo) & (a_[:, 0] < hi), 2])):5.1f} ({((a_[:, 0] >= lo) & (a_[:, 0] < hi)).sum()})" if ((a_[:, 0] >= lo) & (a_[:, 0] < hi)).sum() > 100 else "  -  " for lo, hi in bins))
    print("\n".join(L)); open(os.path.join(H, "shape_tilt.txt"), "w").write("\n".join(L) + "\n")
