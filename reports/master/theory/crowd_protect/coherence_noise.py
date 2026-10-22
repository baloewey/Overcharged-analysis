#!/usr/bin/env python3
"""Is the low tilt coherence of crowded layers noise or structure? (user, 2026-10-22: 'the flow sits closer to the wall, thus smaller and more prone to MPCD noise?')
Shape-sweep circles (runs/shapes_server/circle, 135 runs), frames every 10 t.u. (as escape_events.py). Per frame: wall +1/2 (depth < 4) sorted by arc position;
tilt tl = orientation - outward normal (wrapped). Coherence C = |<exp(i tl)>| (as shape_runs.csv); majority fraction (sign of tl); nearest-neighbour sign
correlation <s_j s_j+1> - <s>^2 normalised (0 = independent flips, > 0 = domains); number of sign domains around the ring. Flow (wallflow.npz): signed
tangential flow <v_t> and <|v_t|> at 1-3 cells (coherent fraction) and the run circulation (shape_runs.csv). -> coherence_noise.txt"""
import csv, glob, json, os, sys
import numpy as np
from multiprocessing import Pool
from scipy.spatial import cKDTree
ROOT = "/home/baloewe/Physics/Overcharged/runs"; sys.path.insert(0, ROOT + "/cavity")
from edge_criteria import read_with_angle  # noqa
from wall_fourier import wall_polyline  # noqa
SW = ROOT + "/shapes_server"; H = os.path.dirname(os.path.abspath(__file__))
def one(rel):
    run = os.path.join(SW, rel); bc = json.load(open(run + "/input.json"))["BC"][0]; Q = bc["Q"]; W, s, P = wall_polyline(bc); tree = cKDTree(W)
    Wc = W[:, 0] + 1j * W[:, 1]; T = np.roll(Wc, -1) - np.roll(Wc, 1); T /= np.abs(T); Nn = -1j * T
    if np.mean((Nn * np.conj(Wc - Wc.mean())).real) < 0: Nn = -Nn
    t, fr = read_with_angle(run, np.array(Q[:2], float), 0.5); fr = fr[::5]; out = []
    for f in fr:
        if not len(f): continue
        p = f[f[:, 2] > 0]; d, j = tree.query(p[:, :2]); w = d < 4
        if w.sum() < 4: continue
        tl = np.angle(np.exp(1j * (p[w, 3] - np.angle(Nn[j[w]])))); o = np.argsort(np.asarray(s)[j[w]]); tl = tl[o]; sg = np.sign(tl)
        C = abs(np.mean(np.exp(1j * tl))); maj = max((sg > 0).mean(), (sg < 0).mean())
        var = 1 - np.mean(sg) ** 2; nn = (np.mean(sg * np.roll(sg, -1)) - np.mean(sg) ** 2) / var if var > 1e-9 else np.nan
        out.append((C, maj, nn, np.sum(sg != np.roll(sg, -1)), np.median(np.abs(tl)), w.sum()))
    return rel, np.array(out)
if __name__ == "__main__":
    runs = sorted(os.path.relpath(os.path.dirname(f), SW) for f in glob.glob(SW + "/circle/*/*/*/s*/DONE"))
    with Pool(8) as p: res = p.map(one, runs, chunksize=2)
    SR = {r["run"]: r for r in csv.DictReader(open(H + "/../shapes/shape_runs.csv"))}; LV = np.array([0.04, 0.06, 0.09, 0.12, 0.14])
    Wf = np.load(H + "/wallflow.npz"); wf = {r: a for r, a in zip(Wf["runs"], Wf["prof"])}; dc = 0.5 * (Wf["DB"][:-1] + Wf["DB"][1:]); band = (dc > 1) & (dc < 3)
    L = ["COHERENCE: NOISE OR STRUCTURE? Shape-sweep circles, per-frame values averaged over frames, then medians over runs per rho level",
         "rho | runs | C | majority fraction | NN sign correlation (0 = independent) | sign domains / ring (random expectation N/2) | N wall | median |tilt| | random-sign C floor cos|tilt| | |<v_t>|/<|v_t|> at 1-3 cells | |circulation|"]
    for i, r0 in enumerate(LV):
        rs = [(r, a) for r, a in res if len(a) and abs(float(SR[r]["rho"]) - r0) == min(abs(float(SR[r]["rho"]) - LV))]
        M = np.array([np.nanmean(a, 0) for r, a in rs]); coh_flow = [abs(wf[r][4][band].mean()) / wf[r][2][band].mean() for r, a in rs]
        circ = [float(SR[r]["circ"]) for r, a in rs]
        L.append(f"{r0:.2f} | {len(rs)} | {np.median(M[:, 0]):.2f} | {np.median(M[:, 1]):.2f} | {np.median(M[:, 2]):+.2f} | {np.median(M[:, 3]):.1f} | {np.median(M[:, 5]):.1f} | "
                 f"{np.median(M[:, 4]) * 180 / np.pi:.0f} | {np.cos(np.median(M[:, 4])):.2f} | {np.median(coh_flow):.2f} | {np.median(circ):.3f}")
    open(H + "/coherence_noise.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
