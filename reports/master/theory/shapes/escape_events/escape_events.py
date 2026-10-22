#!/usr/bin/env python3
"""Per-event escape data, shape sweep (runs/shapes_server, 675 runs; frames every 10 t.u. as kinetics_carriers.py).
Every wall +1/2 (depth < 4) in every frame is followed 60 t.u. (late = still in the wall, depth < 4, at +40 t.u. and depth > 6 at +60 t.u.: a lagged outcome
free of the ~40 t.u. pre-escape run-up). Also followed 20 t.u. (nearest-neighbour chain, < 6 cells per 10 t.u.); escape = depth > 6 at the end.
Per event: run index, local tilt (defect orientation - outward wall normal at the nearest wall point, wrapped to (-pi, pi], deg), arc gaps to the two
nearest wall +1/2 neighbours (cells), depth, |curvature| x R_eq, outcome. Writes events.npz + runs.csv (run index -> family, aspect, R_eq, k, rho)."""
import csv, glob, json, os, sys
import numpy as np
from multiprocessing import Pool
from scipy.spatial import cKDTree
ROOT = "/home/baloewe/Physics/Overcharged/runs"; sys.path.insert(0, ROOT + "/cavity")
from edge_criteria import read_with_angle  # noqa
from wall_fourier import wall_polyline  # noqa
H = os.path.dirname(os.path.abspath(__file__)); SW = ROOT + "/shapes_server"
def one(rel):
    run = os.path.join(SW, rel)
    try:
        bc = json.load(open(run + "/input.json"))["BC"][0]; Q = bc["Q"]
        W, s, P = wall_polyline(bc); n = len(W); Wc = W[:, 0] + 1j * W[:, 1]
        T = np.roll(Wc, -1) - np.roll(Wc, 1); T /= np.abs(T); Nn = -1j * T; kap = np.gradient(np.unwrap(np.angle(Nn)), P / n)
        Req = np.sqrt(abs(0.5 * np.sum(W[:, 0] * np.roll(W[:, 1], -1) - np.roll(W[:, 0], -1) * W[:, 1])) / np.pi)
        tree = cKDTree(W); sarc = np.asarray(s) * P / np.max(s) if np.max(s) > 1.01 else np.asarray(s) * P
        t, fr = read_with_angle(run, np.array(Q[:2], float), 0.5); fr = fr[::5]
        PL, DEP, IDX = [], [], []
        for f in fr:
            p = f[f[:, 2] > 0] if len(f) else np.zeros((0, 4))
            d, i = tree.query(p[:, :2]) if len(p) else (np.zeros(0), np.zeros(0, int))
            PL.append(p); DEP.append(d); IDX.append(i)
        out = []
        for a in range(len(PL) - 6):
            w = DEP[a] < 4
            if w.sum() < 2: continue
            p0 = PL[a][w]; i0 = IDX[a][w]; sa = sarc[i0]
            o = np.argsort(sa); ss = sa[o]; g = np.diff(np.r_[ss, ss[0] + P]); gl = np.empty(len(ss)); gr = np.empty(len(ss))
            gr[o] = g; gl[o] = np.roll(g, 1)
            tl = np.degrees(np.angle(np.exp(1j * (p0[:, 3] - np.angle(Nn[i0])))))
            pos = p0[:, :2]; alive = np.ones(len(pos), bool); dfin = np.zeros(len(pos))
            for b in (1, 2):
                if not len(PL[a + b]): alive[:] = False; break
                dd, jj = cKDTree(PL[a + b][:, :2]).query(pos); alive &= dd < 6; pos = PL[a + b][jj, :2]; dfin = DEP[a + b][jj]
            esc = alive & (dfin > 6)
            for b in (3, 4, 5, 6):
                if not len(PL[a + b]): alive[:] = False; break
                dd, jj = cKDTree(PL[a + b][:, :2]).query(pos); alive &= dd < 6; pos = PL[a + b][jj, :2]; dfin = DEP[a + b][jj]
                if b == 4: alive &= dfin < 4
            late = alive & (dfin > 6)
            out.append(np.column_stack([np.full(len(pos), a), tl, gl, gr, DEP[a][w], np.abs(kap[i0]) * Req, esc, np.full(len(pos), w.sum()), late]))
        return rel, (np.vstack(out) if out else np.zeros((0, 9))), None
    except Exception as ex:
        return rel, None, str(ex)
if __name__ == "__main__":
    man = {r["dir"]: r for r in csv.DictReader(open(SW + "/manifest_shape_sweep.tsv"), delimiter="\t")}
    runs = sorted(os.path.relpath(os.path.dirname(f), SW) for f in glob.glob(SW + "/*/*/*/*/s*/DONE"))
    with Pool(8) as pool: res = pool.map(one, runs, chunksize=4)
    E, meta = [], []
    for j, (rel, ev, err) in enumerate(res):
        m = man.get(rel, {}); meta.append(dict(idx=j, run=rel, family=m.get("family"), aspect=m.get("aspect"), Req_la=m.get("Req_la"), err=err or ""))
        if ev is not None and len(ev): E.append(np.column_stack([np.full(len(ev), j), ev]))
    E = np.vstack(E)
    np.savez_compressed(os.path.join(H, "events.npz"), E=E, cols=np.array(["run", "frame", "tilt", "gap_l", "gap_r", "depth", "kapR", "esc", "nwall", "late"]))
    with open(os.path.join(H, "runs.csv"), "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=list(meta[0])); w.writeheader(); w.writerows(meta)
    print(f"{len(E)} events, {int(E[:, 7].sum())} escapes, errors {sum(1 for m in meta if m['err'])}", flush=True)
