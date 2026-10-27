#!/usr/bin/env python3
"""Flow quantities for the 10-seed +k shape sweep, run ON THE SERVER (flow fields stay there). Per run (second half of the flow frames, every 100 t.u.):
circ = mean |circulation| (as shape_sweep_analysis.py), two = fraction of frames with two counter-rotating halves (as there), and the near-wall flow profile
(as crowd_protect/wallflow.py): per 0.5-cell bin of wall distance 0-12: mean inward v, rms inward v, mean |v_t|, rms |v|, signed v_t.
Usage (server): ~/Overcharged/venv_analysis/bin/python flow_server.py -> flow10.csv, wallflow10.npz"""
import csv, glob, gzip, json, os, sys, numpy as np
from multiprocessing import Pool
from scipy.spatial import cKDTree
from matplotlib.path import Path
def wall_polyline(bc, n=6000):                                           # copy of runs/cavity/wall_fourier.wall_polyline (no local imports on the server)
    th = np.linspace(0, 2 * np.pi, n, endpoint=False)
    if bc.get("shape", 0) == 1:
        a, b = bc.get("cassiniA", 0.0), bc["R"]; r2 = a * a * np.cos(2 * th) + np.sqrt(b ** 4 - a ** 4 * np.sin(2 * th) ** 2); x, y = np.sqrt(r2) * np.cos(th), np.sqrt(r2) * np.sin(th)
    else:
        A, B = bc["aInv"][0] * bc["R"], bc["aInv"][1] * bc["R"]; x, y = A * np.cos(th), B * np.sin(th)
    W = np.column_stack([x, y]); seg = np.hypot(*np.diff(np.vstack([W, W[:1]]), axis=0).T); return W, np.concatenate([[0], np.cumsum(seg)[:-1]]), seg.sum()
SW = os.path.expanduser("~/Overcharged/runs/shapes"); DB = np.arange(0, 12.01, 0.5)
def one(run):
    try:
        bc = json.load(open(run + "/input.json"))["BC"][0]; Q = bc["Q"]; W, s, P = wall_polyline(bc); Wc = W[:, 0] + 1j * W[:, 1]
        T = np.roll(Wc, -1) - np.roll(Wc, 1); T /= np.abs(T)
        if np.sum(Wc.real * np.roll(Wc.imag, -1) - np.roll(Wc.real, -1) * Wc.imag) < 0: T = -T
        Nn = -1j * T; tree = cKDTree(W); path = Path(W)
        rows = np.loadtxt((l for l in gzip.open(run + "/flowfield.dat.gz", "rt") if len(l.split()) == 7 and l.split()[0][0].isdigit()))
        ts = np.unique(rows[:, 0]); ts = ts[ts > ts.max() / 2]
        m0 = rows[:, 0] == ts[0]; x = rows[m0, 1] + 0.5 - Q[0]; y = rows[m0, 2] + 0.5 - Q[1]; ins = path.contains_points(np.column_stack([x, y]))
        d, j = tree.query(np.column_stack([x, y])); b = np.digitize(d, DB) - 1; sel = ins & (b >= 0) & (b < len(DB) - 1); nb = len(DB) - 1
        acc = np.zeros((5, nb)); cnt = np.zeros(nb); circ, two = [], []; r = np.hypot(x, y)
        for tt in ts:
            m = rows[:, 0] == tt; vx, vy = rows[m, 4], rows[m, 5]; v = vx + 1j * vy
            if len(vx) != len(x): continue
            vin = -(v * np.conj(Nn[j])).real; vt = (v * np.conj(T[j])).real
            for q, arr in enumerate((vin, vin ** 2, np.abs(vt), np.abs(v) ** 2, vt)): acc[q] += np.bincount(b[sel], arr[sel], nb)
            cnt += np.bincount(b[sel], None, nb)
            if tt in ts[::2]:
                circ.append(np.mean(((x * vy - y * vx) / np.maximum(r, 0.5))[ins]))
                def hr(se):
                    cx, cy = x[se].mean(), y[se].mean(); dx, dy = x[se] - cx, y[se] - cy
                    return np.mean((dx * vy[se] - dy * vx[se]) / np.maximum(np.hypot(dx, dy), 0.5))
                two.append(hr(ins & (x < 0)) * hr(ins & (x > 0)) < 0)
        mean = acc / np.maximum(cnt, 1)
        return os.path.relpath(run, SW), np.mean(np.abs(circ)), np.mean(two), np.vstack([mean[0], np.sqrt(mean[1]), mean[2], np.sqrt(mean[3]), mean[4]]), None
    except Exception as ex:
        return os.path.relpath(run, SW), np.nan, np.nan, None, str(ex)
if __name__ == "__main__":
    runs = sorted(os.path.dirname(f) for f in glob.glob(SW + "/*/*/*/*/s*/DONE"))
    with Pool(int(os.environ.get("NPROC", 64))) as p: res = p.map(one, runs, chunksize=2)
    with open("flow10.csv", "w", newline="") as fo:
        w = csv.writer(fo); w.writerow(["run", "circ", "two", "err"]); [w.writerow([r, c, t, e or ""]) for r, c, t, a, e in res]
    ok = [(r, a) for r, c, t, a, e in res if a is not None]
    np.savez_compressed("wallflow10.npz", runs=np.array([r for r, a in ok]), prof=np.array([a for r, a in ok]), DB=DB)
    print("runs", len(res), "errors", sum(1 for x in res if x[4]))
