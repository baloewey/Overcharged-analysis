#!/usr/bin/env python3
"""Test 2 (2026-10-22): near-wall flow vs wall density, shape sweep (runs/shapes_server, 675 runs; flowfield.dat.gz every 100 t.u., second half).
Each lattice cell inside the wall: distance d to the wall polyline, outward normal n and ccw tangent t at the nearest wall point. Per run, in 0.5-cell
bins of d (0-12 cells): mean inward velocity <-v.n>, rms inward velocity, mean |v.t|, rms |v|, and the signed tangential flow <v.t>.
-> wallflow.npz (per-run profiles) and wallflow.txt (medians over runs per rho level and shape class)."""
import csv, glob, gzip, json, os, sys
import numpy as np
from multiprocessing import Pool
from scipy.spatial import cKDTree
from matplotlib.path import Path
ROOT = "/home/baloewe/Physics/Overcharged/runs"; sys.path.insert(0, ROOT + "/cavity")
from wall_fourier import wall_polyline  # noqa
SW = ROOT + "/shapes_server"; H = os.path.dirname(os.path.abspath(__file__)); DB = np.arange(0, 12.01, 0.5)
def one(rel):
    run = os.path.join(SW, rel)
    try:
        bc = json.load(open(run + "/input.json"))["BC"][0]; Q = bc["Q"]; W, s, P = wall_polyline(bc)
        Wc = W[:, 0] + 1j * W[:, 1]; T = np.roll(Wc, -1) - np.roll(Wc, 1); T /= np.abs(T)
        if np.sum(Wc.real * np.roll(Wc.imag, -1) - np.roll(Wc.real, -1) * Wc.imag) < 0: T = -T          # make T ccw
        Nn = -1j * T; tree = cKDTree(W); path = Path(W)
        rows = np.loadtxt((l for l in gzip.open(run + "/flowfield.dat.gz", "rt") if len(l.split()) == 7 and l.split()[0][0].isdigit()))
        ts = np.unique(rows[:, 0]); ts = ts[ts > ts.max() / 2]
        m0 = rows[:, 0] == ts[0]; x = rows[m0, 1] + 0.5 - Q[0]; y = rows[m0, 2] + 0.5 - Q[1]
        ins = path.contains_points(np.column_stack([x, y])); d, j = tree.query(np.column_stack([x, y])); b = np.digitize(d, DB) - 1
        sel = ins & (b >= 0) & (b < len(DB) - 1); nb = len(DB) - 1
        acc = np.zeros((5, nb)); cnt = np.zeros(nb)
        for tt in ts:
            m = rows[:, 0] == tt; v = rows[m, 4] + 1j * rows[m, 5]
            vin = -(v * np.conj(Nn[j])).real; vt = (v * np.conj(T[j])).real
            for q, arr in enumerate((vin, vin ** 2, np.abs(vt), np.abs(v) ** 2, vt)): acc[q] += np.bincount(b[sel], arr[sel], nb)
            cnt += np.bincount(b[sel], None, nb)
        mean = acc / np.maximum(cnt, 1)
        return rel, np.vstack([mean[0], np.sqrt(mean[1]), mean[2], np.sqrt(mean[3]), mean[4]]), None
    except Exception as ex:
        return rel, None, str(ex)
if __name__ == "__main__":
    runs = sorted(os.path.relpath(os.path.dirname(f), SW) for f in glob.glob(SW + "/*/*/*/*/s*/DONE"))
    with Pool(int(os.environ.get("NPROC", 8))) as p: res = p.map(one, runs, chunksize=2)
    ok = [(r, a) for r, a, e in res if a is not None]; print("errors", sum(1 for r in res if r[2]), [r[2] for r in res if r[2]][:3])
    np.savez_compressed(os.path.join(H, "wallflow.npz"), runs=np.array([r for r, a in ok]), prof=np.array([a for r, a in ok]), DB=DB)
