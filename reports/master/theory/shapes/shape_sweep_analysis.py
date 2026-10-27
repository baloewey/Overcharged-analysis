#!/usr/bin/env python3
"""Shape sweep (server, 675 runs: circles, Cassini ovals b/a 1.1-3.0 (eMode 3), ellipses eps 1.2-2.5 (eMode 2); R_eq = 1.0, 1.5, 2.0 l_a; rho 0.04-0.14;
3 seeds, extended to 10 seeds (s04-s10, server 16-26 Oct; labels from the s01 manifest row); act 0.07), analysed with the circle framework. Per run (second half; defects every 2 t.u., subsampled every 5th frame; flow every 100 t.u.):
  topology  wall +1/2 (depth < 4): mean, mode; vacancies 2k - 2 - mean; bulk charge (depth > 6)
  lock      anchoring relative to the normal psi_rel(s) = k 2 pi s/P - phi_n(s) (arclength-uniform anchoring, s from the +x axis); concentration of
            tail tilt - 2 psi_rel(s) over wall defects (1 = perfect local lock)
  registry  distance of each wall defect to the nearest homeotropic point (psi_rel = 0 mod pi), rms; |tilt|, tilt coherence (|mean exp(i tilt)|)
  escapes   wall +1/2 followed 20 t.u. (nearest-neighbour chain, < 3 cells per frame) ending deeper than 6 cells: probability per wall defect per 20 t.u.,
            binned by the local wall curvature kappa R_eq at the defect's starting position (pooled table written separately)
  flow      |circulation|, fraction of intervals with two vortices (halves along the long axis counter-rotating)
Writes shape_runs.csv (one row per run) and shape_escape_curvature.csv (pooled escapes by curvature bin, per family)."""
import csv, glob, gzip, json, os, sys
import numpy as np
from multiprocessing import Pool
from scipy.spatial import cKDTree
ROOT = "/home/baloewe/Physics/Overcharged/runs"; sys.path.insert(0, ROOT + "/cavity")
from edge_criteria import read_with_angle  # noqa
from wall_fourier import wall_polyline  # noqa
H = os.path.dirname(os.path.abspath(__file__)); SW = ROOT + "/shapes_server"; LA = 22.21
KBINS = np.array([0, 0.5, 0.8, 1.0, 1.25, 1.6, 2.0, 3.0, 5.0, 100.0])
def geometry(bc):
    W, s, P = wall_polyline(bc); n = len(W)
    Wc = W[:, 0] + 1j * W[:, 1]; T = np.roll(Wc, -1) - np.roll(Wc, 1); T /= np.abs(T); Nn = -1j * T          # outward normal (ccw wall)
    phin = np.unwrap(np.angle(Nn)); ds = P / n
    kap = np.gradient(phin, ds)
    return Wc, s, P, Nn, phin, kap
def one(run):
    try:
        J = json.load(open(run + "/input.json")); bc = J["BC"][0]; k = bc["customAnch"]; Q = bc["Q"]
        Wc, s, P, Nn, phin, kap = geometry(bc); psi = k * 2 * np.pi * s / P - phin
        Req = np.sqrt(abs(0.5 * np.sum((Wc.real * np.roll(Wc.imag, -1) - np.roll(Wc.real, -1) * Wc.imag))) / np.pi)
        tree = cKDTree(np.column_stack([Wc.real, Wc.imag]))
        t, fr = read_with_angle(run, np.array(Q[:2], float), 0.5); fr = fr[::5]; t = np.array(t)[::5]
        nw, qb, lockv, dist, tilts, esc_bins, start_bins = [], [], [], [], [], np.zeros(len(KBINS) - 1), np.zeros(len(KBINS) - 1)
        PL = []; DEP = []; IDX = []
        for f in fr:
            p = f[f[:, 2] > 0]; d, i = tree.query(p[:, :2]) if len(p) else (np.array([]), np.array([], int))
            inside_sign = 1.0
            PL.append(p); DEP.append(d); IDX.append(i)
            all_d, all_i = tree.query(f[:, :2]) if len(f) else (np.array([]), np.array([], int))
            qb.append(0.5 * np.sum(np.sign(f[all_d > 6, 2])) if len(f) else 0.0)
            w = d < 4; nw.append(w.sum())
            if w.any():
                tl = np.angle(np.exp(1j * (p[w, 3] - np.angle(Nn[i[w]])))); tilts.append(np.mean(np.exp(1j * tl)))
                lockv.extend(np.exp(1j * (tl - 2 * psi[i[w]])))
                dpsi = np.angle(np.exp(2j * psi[i[w]])) / 2; slope = np.gradient(psi, P / len(s))[i[w]]; dist.extend(np.abs(dpsi / np.where(np.abs(slope) > 1e-6, slope, 1e-6)))
        # escapes: 20 t.u. = 2 subsampled frames (dt 10)
        NF = 2
        for a in range(len(fr) - NF):
            w = DEP[a] < 4
            if not w.any(): continue
            pos = PL[a][w, :2]; alive = np.ones(len(pos), bool); k0 = np.abs(kap[IDX[a][w]]) * Req
            for b in range(1, NF + 1):
                if not len(PL[a + b]): alive[:] = False; break
                dd, jj = cKDTree(PL[a + b][:, :2]).query(pos); alive &= dd < 6; pos = PL[a + b][jj, :2]; dfin = DEP[a + b][jj]
            e = alive & (dfin > 6); bins = np.clip(np.digitize(k0, KBINS) - 1, 0, len(KBINS) - 2)
            np.add.at(start_bins, bins, 1); np.add.at(esc_bins, bins[e], 1)
        nw = np.array(nw); out = dict(run=os.path.relpath(run, SW), k=k, Req=Req / LA, nwall=nw.mean(), mode=int(np.bincount(nw).argmax()) if len(nw) else -1,
                                       vac=(2 * k - 2) - nw.mean(), qbulk=np.mean(qb), lock=abs(np.mean(lockv)) if lockv else np.nan,
                                       dist=np.sqrt(np.mean(np.square(dist))) if dist else np.nan, tilt=np.degrees(np.mean(np.abs(np.angle(tilts)))) if tilts else np.nan,
                                       coh=np.mean(np.abs(tilts)) if tilts else np.nan, esc=esc_bins.sum() / max(start_bins.sum(), 1))
        ff = run + "/flowfield.dat.gz"
        if os.path.exists(ff):
            rows = np.loadtxt((l for l in gzip.open(ff, "rt") if len(l.split()) == 7 and l.split()[0][0].isdigit()))
            ts = np.unique(rows[:, 0]); ts = ts[ts > ts.max() / 2]; xs = rows[:, 1] + 0.5 - Q[0]; ys = rows[:, 2] + 0.5 - Q[1]
            path = __import__("matplotlib.path", fromlist=["Path"]).Path(np.column_stack([Wc.real, Wc.imag]))
            circ, two = [], []
            for tt in ts[::2]:
                m = rows[:, 0] == tt; x, y = xs[m], ys[m]; vx, vy = rows[m, 4], rows[m, 5]
                ins = path.contains_points(np.column_stack([x, y])); r = np.hypot(x, y)
                circ.append(np.mean(((x * vy - y * vx) / np.maximum(r, 0.5))[ins]))
                def hr(sel):
                    cx, cy = x[sel].mean(), y[sel].mean(); dx, dy = x[sel] - cx, y[sel] - cy
                    return np.mean((dx * vy[sel] - dy * vx[sel]) / np.maximum(np.hypot(dx, dy), 0.5))
                two.append(hr(ins & (x < 0)) * hr(ins & (x > 0)) < 0)
            out.update(circ=np.mean(np.abs(circ)), two=np.mean(two))
        return out, esc_bins, start_bins
    except Exception as ex:
        return dict(run=os.path.relpath(run, SW), error=str(ex)), None, None
if __name__ == "__main__":
    man = {r["dir"]: r for r in csv.DictReader(open(SW + "/manifest_shape_sweep.tsv"), delimiter="\t")}
    runs = sorted(os.path.dirname(f) for f in glob.glob(SW + "/*/*/*/*/s*/DONE"))
    for rr in runs:                                                       # seeds s04-s10 (manifest_seeds10.tsv has no shape columns): labels from s01
        rel = os.path.relpath(rr, SW)
        if rel not in man: man[rel] = man.get(os.path.dirname(rel) + "/s01", {})
    print(len(runs), "runs", flush=True)
    with Pool(12) as pool: res = pool.map(one, runs, chunksize=4)
    cols = ["run", "family", "shape", "aspect", "Req_la", "k", "rho", "nwall", "mode", "vac", "qbulk", "lock", "dist", "tilt", "coh", "esc", "circ", "two"]
    with open(os.path.join(H, "shape_runs.csv"), "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=cols, extrasaction="ignore"); w.writeheader()
        for o, _, _ in res:
            m = man.get(o["run"], {}); o.update(family=m.get("family"), shape=m.get("shape"), aspect=m.get("aspect"), Req_la=m.get("Req_la"), rho=m.get("rho")); w.writerow(o)
    pooled = {}
    for o, e, s_ in res:
        if e is None: continue
        fam = man.get(o["run"], {}).get("family", "?"); pe, ps = pooled.get(fam, (np.zeros_like(e), np.zeros_like(s_))); pooled[fam] = (pe + e, ps + s_)
    with open(os.path.join(H, "shape_escape_curvature.csv"), "w") as fo:
        fo.write("family,kappaR_lo,kappaR_hi,starts,escapes,prob\n")
        for fam, (e, s_) in pooled.items():
            for i in range(len(e)): fo.write(f"{fam},{KBINS[i]},{KBINS[i + 1]},{int(s_[i])},{int(e[i])},{e[i] / max(s_[i], 1):.5f}\n")
    print("errors:", sum(1 for o, _, _ in res if "error" in o))
