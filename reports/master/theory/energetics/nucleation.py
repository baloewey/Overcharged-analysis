#!/usr/bin/env python3
"""Does the bulk overcharge affect dynamical defect nucleation? (user, 2026-10-13). Shape sweep, second half, every frame (2 t.u.).
Per frame: bulk net charge Q_b (defects >= 4 cells from the wall), bulk -1/2 count (each is half of a dynamical pair, since wall -1/2 are rare for k > 0).
Pair creation events: a new -1/2 (no -1/2 within 2.5 cells in the previous frame) with a new +1/2 within 3 cells (no +1/2 within 2.5 cells of it in the previous frame), in the bulk.
Per run: creation rate per bulk area per 1000 t.u., mean bulk -1/2 density, mean Q_b and its density. Across runs: by family/aspect/size; within runs: creation rate vs dQ_b = Q_b - run mean.
Usage: nucleation.py [NPROC] -> nucleation.npy, nucleation.txt"""
import os, sys, csv, json, numpy as np
from multiprocessing import Pool
from scipy.spatial import cKDTree
from matplotlib.path import Path
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
import energetics as En  # noqa
SW = "/home/baloewe/Physics/Overcharged/runs/shapes_server"
def one(rel):
    try:
        run = os.path.join(SW, rel); bc = json.load(open(run + "/input.json"))["BC"][0]; k = bc["customAnch"]; W, s, P = En.wall_polyline(bc); tree = cKDTree(W)
        xs = np.arange(W[:, 0].min(), W[:, 0].max(), 0.5); ys = np.arange(W[:, 1].min(), W[:, 1].max(), 0.5); XX, YY = np.meshgrid(xs, ys); pts = np.c_[XX.ravel(), YY.ravel()]
        ins = Path(W).contains_points(pts); area = ins.sum() * 0.25; bulk_area = (ins & (tree.query(pts)[0] >= 4)).sum() * 0.25
        t, fr = En.read_with_angle(run, np.array(bc["Q"][:2], float), 0.5); t = np.array(t); dt = float(np.median(np.diff(t)))
        qb, nm, ev = [], [], []; prev = None
        for f in fr:
            if not len(f): qb.append(np.nan); nm.append(np.nan); ev.append(0); prev = None; continue
            d = tree.query(f[:, :2])[0]; b = d >= 4; qb.append(0.5 * np.sign(f[b, 2]).sum()); nm.append(((f[:, 2] < 0) & b).sum()); e = 0
            if prev is not None:
                pm = prev[prev[:, 2] < 0, :2]; pp = prev[prev[:, 2] > 0, :2]
                newm = [i for i in np.flatnonzero((f[:, 2] < 0) & b) if not len(pm) or np.hypot(*(pm - f[i, :2]).T).min() > 2.5]
                newp = [i for i in np.flatnonzero((f[:, 2] > 0) & b) if not len(pp) or np.hypot(*(pp - f[i, :2]).T).min() > 2.5]
                used = set()
                for i in newm:
                    for j in newp:
                        if j not in used and np.hypot(*(f[i, :2] - f[j, :2])) < 3: e += 1; used.add(j); break
            ev.append(e); prev = f
        return dict(rel=rel, k=k, P=P, area=area, bulk_area=bulk_area, dt=dt, qb=np.array(qb, float), nm=np.array(nm, float), ev=np.array(ev, float))
    except Exception as ex:
        return dict(rel=rel, err=str(ex))
if __name__ == "__main__":
    man = {r["dir"]: r for r in csv.DictReader(open(SW + "/manifest_shape_sweep.tsv"), delimiter="\t")}
    rels = [r for r in man if os.path.exists(os.path.join(SW, r, "DONE"))]
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 8) as pool: res = [r for r in pool.map(one, rels, chunksize=4) if "err" not in r]
    np.save(os.path.join(H, "nucleation.npy"), np.array(res, dtype=object), allow_pickle=True)
    rows = []
    for r in res:
        m = man[r["rel"]]; T = np.isfinite(r["qb"]).sum() * r["dt"]
        rows.append((m["family"], float(m["aspect"]), float(m["Req_la"]), r["k"], np.nanmean(r["qb"]), np.nanmean(r["qb"]) / r["bulk_area"], 1000 * r["ev"].sum() / (r["bulk_area"] * T), np.nanmean(r["nm"]) / r["bulk_area"]))
    A = np.array(rows, dtype=object); fam = A[:, 0].astype(str); asp, req, k, Qb, rhoQ, rate, rhom = (A[:, i].astype(float) for i in range(1, 8))
    L = ["ACROSS RUNS (means): group | n | bulk charge | bulk charge density (x1e-4) | pair creation per bulk area per 1000 t.u. (x1e-3) | bulk -1/2 density (x1e-4)"]
    for R in (1.0, 1.5, 2.0):
        for lab, sel in (("circle", fam == "circle"), ("aspect 1.1-1.3", (fam != "circle") & (asp < 1.3)), ("aspect >= 1.4", asp >= 1.38)):
            s = sel & (req == R)
            if s.sum(): L.append(f"  R_eq {R} {lab:15s} | {s.sum():3d} | {Qb[s].mean():+.2f} | {1e4 * rhoQ[s].mean():5.2f} | {1e3 * rate[s].mean():6.2f} | {1e4 * rhom[s].mean():5.2f}")
    # partial correlation of the creation rate with the bulk charge density, controlling for size (R_eq) and k
    X = np.c_[np.ones(len(A)), rhoQ, req, k, asp]; b = np.linalg.lstsq(X, rate, rcond=None)[0]
    L.append(f"\nregression across runs: rate = a + b1*rhoQ + b2*R_eq + b3*k + b4*aspect: b1 = {b[1]:+.4f} (rate units per unit charge density), b3 (k) = {b[3]:+.2e}, b4 (aspect) = {b[4]:+.2e}")
    # within runs: creation rate in frames binned by dQ_b
    L.append("\nWITHIN RUNS (elongated aspect >= 1.4, and circles): pair creation per frame relative to the run mean, binned by dQ_b = Q_b - run mean")
    for lab, selfn in (("elongated", lambda m: float(m["aspect"]) >= 1.38), ("circle", lambda m: m["family"] == "circle")):
        Bn = np.arange(-2.25, 2.3, 0.5); num = np.zeros(len(Bn) - 1); den = np.zeros(len(Bn) - 1); nmn = np.zeros(len(Bn) - 1)
        for r in res:
            if not selfn(man[r["rel"]]): continue
            ok = np.isfinite(r["qb"]); q = r["qb"][ok]; e = r["ev"][ok]; nmm = r["nm"][ok]
            if e.mean() == 0: continue
            bi = np.digitize(q - q.mean(), Bn) - 1
            for j in range(len(Bn) - 1):
                s = bi == j; num[j] += (e[s] / e.mean()).sum(); den[j] += s.sum(); nmn[j] += (nmm[s] / max(nmm.mean(), 1e-9)).sum()
        L.append(f"  {lab}: dQ_b " + "  ".join(f"{0.5 * (Bn[j] + Bn[j + 1]):+.1f}: rate {num[j] / den[j]:.2f}, -1/2 {nmn[j] / den[j]:.2f} (n {int(den[j])})" for j in range(len(Bn) - 1) if den[j] > 500))
    print("\n".join(L)); open(os.path.join(H, "nucleation.txt"), "w").write("\n".join(L) + "\n")
