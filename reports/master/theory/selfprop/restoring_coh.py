#!/usr/bin/env python3
"""Does the restoring force on wall +1/2 follow the flow coherence or the cavity class? Shape sweep, second half, frames every 10 t.u.
Per frame: tilt coherence of the wall layer |<exp(i tau)>| (wall +1/2, depth < 4), smoothed over +-5 frames (100 t.u.). Per wall +1/2: depth d and its depth change
over 10 t.u. (nearest-neighbour step < 6). Pooled per (class, coherence bin): effective force -U'/kT = 2<dd>/var(dd) per depth bin, and the well U(3.75) - U_min.
Also the imposed-swirl ellipses (runs/swirl_ellipse*, active vs passive) per run. Writes restoring_coh.txt."""
import csv, glob, json, os, sys
import numpy as np
from multiprocessing import Pool
from scipy.spatial import cKDTree
ROOT = "/home/baloewe/Physics/Overcharged/runs"; sys.path.insert(0, ROOT + "/cavity")
from edge_criteria import read_with_angle  # noqa
from wall_fourier import wall_polyline  # noqa
H = os.path.dirname(os.path.abspath(__file__)); SW = ROOT + "/shapes_server"
def samples(run, tmin_frac=0.5):
    J = json.load(open(run + "/input.json")); bc = J["BC"][0]; Q = bc["Q"]
    if "aInv" not in bc: bc["aInv"] = [1, 1, 1]
    W, s, P = wall_polyline(bc); tree = cKDTree(W); Wc = W[:, 0] + 1j * W[:, 1]; Tn = np.roll(Wc, -1) - np.roll(Wc, 1); Nn = -1j * Tn / np.abs(Tn)
    if np.mean((Nn * np.conj(Wc)).real) < 0: Nn = -Nn
    t, fr = read_with_angle(run, np.array(Q[:2], float), tmin_frac); fr = fr[::5]
    PL = [f[f[:, 2] > 0] if len(f) else np.zeros((0, 4)) for f in fr]; coh = np.full(len(PL), np.nan)
    for a, p in enumerate(PL):
        if len(p):
            d, j = tree.query(p[:, :2]); w = d < 4
            if w.any(): coh[a] = abs(np.mean(np.exp(1j * (p[w, 3] - np.angle(Nn[j[w]])))))
    cs = np.array([np.nanmean(coh[max(0, a - 5):a + 6]) for a in range(len(coh))]); out = []
    for a in range(len(PL) - 1):
        p, q = PL[a], PL[a + 1]
        if not len(p) or not len(q): continue
        d, _ = tree.query(p[:, :2]); w = d < 4
        if not w.any(): continue
        dq, k = cKDTree(q[:, :2]).query(p[w, :2]); dn, _ = tree.query(q[k, :2]); ok = dq < 6
        out += [(dd, n - dd, cs[a]) for dd, n in zip(d[w][ok], dn[ok])]
    return np.array(out)
def one(rel):
    try: return rel, samples(os.path.join(SW, rel))
    except Exception: return rel, None
edges = np.arange(0, 4.01, 0.5); cen = 0.5 * (edges[1:] + edges[:-1])
def well(a):
    if a is None or len(a) < 2000: return np.nan, None
    F = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        x = a[(a[:, 0] >= lo) & (a[:, 0] < hi), 1]; F.append(2 * x.mean() / x.var() if len(x) > 100 else np.nan)   # -U'/kT
    F = np.array(F); U = np.concatenate([[0], np.cumsum(-0.5 * (F[:-1] + F[1:]) * 0.5)]); U -= np.nanmin(U); return U[-1], F
if __name__ == "__main__":
    man = {r["dir"]: r for r in csv.DictReader(open(SW + "/manifest_shape_sweep.tsv"), delimiter="\t")}
    runs = sorted(os.path.relpath(os.path.dirname(f), SW) for f in glob.glob(SW + "/*/*/*/*/s*/DONE"))
    with Pool(5) as pool: res = pool.map(one, runs, chunksize=4)
    L = []; say = lambda x: (L.append(x), print(x, flush=True))
    bins = ((0, 0.5), (0.5, 0.7), (0.7, 0.85), (0.85, 1.01))
    say("well depth U(3.75) - U_min (kT_eff) and force at d = 2.75 (-U'/kT), by class, R_eq and layer coherence (100 t.u.)")
    for el in (False, True):
        for R in (1.0, 1.5, 2.0):
            A = [a for rel, a in res if a is not None and len(a) and (float(man[rel]["aspect"]) >= 1.25) == el and float(man[rel]["Req_la"]) == R]
            A = np.concatenate(A); line = f"  {'elongated' if el else 'circle-like':11s} R_eq {R}: "
            for lo, hi in bins:
                m = (A[:, 2] >= lo) & (A[:, 2] < hi); Wd, F = well(A[m]); line += f"coh {lo:.2f}-{hi:.2f}: well {Wd:4.2f}, F(2.75) {F[5] if F is not None else np.nan:+5.2f} ({m.mean():.0%})   "
            say(line)
    say("\nimposed-swirl ellipses and circles (second half): well, F(2.75), mean coherence")
    for rel in sorted(glob.glob(ROOT + "/swirl_ellipse*/*") + glob.glob(ROOT + "/swirl_threshold/passive_A0.01_s*")):
        if not os.path.exists(rel + "/DONE"): continue
        try: a = samples(rel); Wd, F = well(a); say(f"  {os.path.relpath(rel, ROOT):40s}: well {Wd:4.2f}  F(2.75) {F[5] if F is not None else np.nan:+5.2f}  coherence {np.nanmean(a[:, 2]):.2f}  (n {len(a)})")
        except Exception as ex: say(f"  {rel}: {ex}")
    open(os.path.join(H, "restoring_coh.txt"), "w").write("\n".join(L) + "\n")
