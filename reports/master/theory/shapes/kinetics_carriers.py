#!/usr/bin/env python3
"""(A) escape/refill kinetics of the wall layer and (B) carriers of the bulk charge, shape sweep (675 runs; second half; frames every 10 t.u.).
(A) every wall +1/2 (depth < 4) and every bulk +1/2 (depth > 6) is followed 20 t.u. (nearest-neighbour chain, < 6 cells per 10 t.u.): escape = wall -> depth > 6,
    refill = bulk -> depth < 4. Per run: escape probability per wall defect, refill probability per bulk +1/2 (per 20 t.u.), event rates per 1000 t.u., and the region
    (|x|/A: tip > 1/2, flank < 1/4, middle) where escapes start and refills arrive.
(B) bulk +1/2 and -1/2 counts; a bulk +1/2 is 'paired' if a -1/2 lies within 4 cells; unpaired +1/2 = carriers of the excess charge; their mean |x|/A and their mean
    'elliptic radius' rho (0 centre, 1 wall; from the wall distance: rho = 1 - depth/Req approx).
Writes kinetics_carriers.csv (one row per run) and prints a summary by family/aspect/R_eq."""
import csv, glob, json, os, sys
import numpy as np
from multiprocessing import Pool
from scipy.spatial import cKDTree
ROOT = "/home/baloewe/Physics/Overcharged/runs"; sys.path.insert(0, ROOT + "/cavity")
from edge_criteria import read_with_angle  # noqa
from wall_fourier import wall_polyline  # noqa
H = os.path.dirname(os.path.abspath(__file__)); SW = ROOT + "/shapes_server"; LA = 22.21
def region(x, A): ax = abs(x) / A; return 0 if ax > 0.5 else (2 if ax < 0.25 else 1)     # 0 tip, 1 middle, 2 flank
def one(rel):
    run = os.path.join(SW, rel)
    try:
        bc = json.load(open(run + "/input.json"))["BC"][0]; k = bc["customAnch"]; Q = bc["Q"]
        W, s, P = wall_polyline(bc); tree = cKDTree(W); A = np.abs(W[:, 0]).max()
        Req = np.sqrt(abs(0.5 * np.sum(W[:, 0] * np.roll(W[:, 1], -1) - np.roll(W[:, 0], -1) * W[:, 1])) / np.pi)
        t, fr = read_with_angle(run, np.array(Q[:2], float), 0.5); fr = fr[::5]
        PL, DEP = [], []
        nw, nbp, nbm, unp, unp_x, unp_rho, paired = [], [], [], [], [], [], []
        for f in fr:
            if not len(f): PL.append(np.zeros((0, 4))); DEP.append(np.zeros(0)); continue
            d, _ = tree.query(f[:, :2]); pl = f[:, 2] > 0
            PL.append(f[pl]); DEP.append(d[pl])
            b = d > 6; bp = f[b & pl]; bm = f[b & ~pl]
            nw.append(((d < 4) & pl).sum()); nbp.append(len(bp)); nbm.append(len(bm))
            if len(bp):
                if len(bm): dd, _ = cKDTree(bm[:, :2]).query(bp[:, :2]); isp = dd < 4
                else: isp = np.zeros(len(bp), bool)
                paired.append(isp.mean()); u = bp[~isp]; unp.append(len(u))
                if len(u):
                    du, _ = tree.query(u[:, :2]); unp_x.extend(np.abs(u[:, 0]) / A); unp_rho.extend(1 - du / Req)
            else: unp.append(0)
        esc = np.zeros(3); esc_start = np.zeros(3); ref = np.zeros(3); ref_start = 0; NF = 2
        for a in range(len(PL) - NF):
            p0, d0 = PL[a], DEP[a]
            for cls in ("wall", "bulk"):
                sel = d0 < 4 if cls == "wall" else d0 > 6
                if not sel.any(): continue
                pos = p0[sel, :2]; alive = np.ones(len(pos), bool); x0 = pos[:, 0].copy()
                for b in range(1, NF + 1):
                    if not len(PL[a + b]): alive[:] = False; break
                    dd, jj = cKDTree(PL[a + b][:, :2]).query(pos); alive &= dd < 6; pos = PL[a + b][jj, :2]; dfin = DEP[a + b][jj]
                if cls == "wall":
                    for x, ok, df in zip(x0, alive, dfin):
                        r_ = region(x, A); esc_start[r_] += 1
                        if ok and df > 6: esc[r_] += 1
                else:
                    ref_start += len(pos)
                    for x, ok, df in zip(pos[:, 0], alive, dfin):
                        if ok and df < 4: ref[region(x, A)] += 1
        T = len(PL) * 10.0
        return dict(run=rel, k=k, nwall=np.mean(nw), vac=(2 * k - 2) - np.mean(nw), nbp=np.mean(nbp), nbm=np.mean(nbm), unpaired=np.mean(unp),
                    paired_frac=np.mean(paired) if paired else np.nan, unp_x=np.mean(unp_x) if unp_x else np.nan, unp_rho=np.mean(unp_rho) if unp_rho else np.nan,
                    p_esc=esc.sum() / max(esc_start.sum(), 1), p_ref=ref.sum() / max(ref_start, 1), r_esc=esc.sum() / T * 1000, r_ref=ref.sum() / T * 1000,
                    esc_tip=esc[0] / max(esc_start[0], 1), esc_mid=esc[1] / max(esc_start[1], 1), esc_flank=esc[2] / max(esc_start[2], 1),
                    ref_tip=ref[0] / max(ref.sum(), 1), ref_mid=ref[1] / max(ref.sum(), 1), ref_flank=ref[2] / max(ref.sum(), 1),
                    tip_share=esc_start[0] / max(esc_start.sum(), 1))
    except Exception as ex:
        return dict(run=rel, error=str(ex))
if __name__ == "__main__":
    man = {r["dir"]: r for r in csv.DictReader(open(SW + "/manifest_shape_sweep.tsv"), delimiter="\t")}
    runs = sorted(os.path.relpath(os.path.dirname(f), SW) for f in glob.glob(SW + "/*/*/*/*/s*/DONE"))
    with Pool(8) as pool: res = pool.map(one, runs, chunksize=4)
    cols = ["run", "family", "aspect", "Req_la", "k", "nwall", "vac", "nbp", "nbm", "unpaired", "paired_frac", "unp_x", "unp_rho", "p_esc", "p_ref", "r_esc", "r_ref",
            "esc_tip", "esc_mid", "esc_flank", "ref_tip", "ref_mid", "ref_flank", "tip_share"]
    with open(os.path.join(H, "kinetics_carriers.csv"), "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=cols, extrasaction="ignore"); w.writeheader()
        for o in res:
            m = man.get(o["run"], {}); o.update(family=m.get("family"), aspect=m.get("aspect"), Req_la=m.get("Req_la")); w.writerow(o)
    print("errors:", sum(1 for o in res if "error" in o), flush=True)
