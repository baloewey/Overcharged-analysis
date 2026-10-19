#!/usr/bin/env python3
"""Tilt sign vs circulation, +k (shapes_server, wall +1/2) and -k (shapes_mirror, wall -1/2), same cavities (user, 2026-10-19).
Per flow frame (every 100 t.u., second half; defect frame nearest in time): circulation Gamma (whole cavity) and per half (x < 0, x > 0, about each half's centroid);
common wall tilt tau_c: +1/2 arg<exp(i tau)>, -1/2 arg<exp(3 i tau)>/3 (tau = solver angle - outward normal; wall = < 4 cells), whole wall and per half.
Agreement = fraction of frames with sign(tau_c) = sign(Gamma) (frames with |coherence| < 0.2 skipped); signed tilt <tau_c sign(Gamma)>.
Writes tilt_sign.csv (one row per run and sign)."""
import csv, glob, gzip, json, os, sys
import numpy as np
from multiprocessing import Pool
from scipy.spatial import cKDTree
from matplotlib.path import Path
ROOT = "/home/baloewe/Physics/Overcharged/runs"; sys.path.insert(0, ROOT + "/cavity")
from edge_criteria import read_with_angle  # noqa
from wall_fourier import wall_polyline  # noqa
H = os.path.dirname(os.path.abspath(__file__))
def one(run):
    try:
        bc = json.load(open(run + "/input.json"))["BC"][0]; k = bc["customAnch"]; Q = bc["Q"]; sg = 1 if k > 0 else -1; m3 = 1 if k > 0 else 3
        W, s, P = wall_polyline(bc); Wc = W[:, 0] + 1j * W[:, 1]; T = np.roll(Wc, -1) - np.roll(Wc, 1); Nn = -1j * T / np.abs(T); tree = cKDTree(W)
        td, fr = read_with_angle(run, np.array(Q[:2], float), 0.5); td = np.array(td)
        rows = np.loadtxt((l for l in gzip.open(run + "/flowfield.dat.gz", "rt") if len(l.split()) == 7 and l.split()[0][0].isdigit()))
        ts = np.unique(rows[:, 0]); ts = ts[ts > ts.max() / 2]; xs = rows[:, 1] + 0.5 - Q[0]; ys = rows[:, 2] + 0.5 - Q[1]; path = Path(W)
        res = {"whole": [], "half": []}
        for tt in ts:
            m = rows[:, 0] == tt; x, y, vx, vy = xs[m], ys[m], rows[m, 4], rows[m, 5]; ins = path.contains_points(np.column_stack([x, y]))
            def gam(sel):
                cx, cy = x[sel].mean(), y[sel].mean(); dx, dy = x[sel] - cx, y[sel] - cy
                return np.mean((dx * vy[sel] - dy * vx[sel]) / np.maximum(np.hypot(dx, dy), 0.5))
            f = fr[int(np.argmin(np.abs(td - tt)))]
            if not len(f): continue
            d, i = tree.query(f[:, :2]); w = (d < 4) & (np.sign(f[:, 2]) == sg)
            if w.sum() < 3: continue
            z = np.exp(1j * m3 * (f[w, 3] - np.angle(Nn[i[w]]))); xw = f[w, 0]
            def tc(sel):
                if sel.sum() < 2: return None
                zz = z[sel].mean(); return (np.angle(zz) / m3, abs(zz)) if abs(zz) > 0.2 else None
            G = gam(ins); a = tc(np.ones(len(z), bool))
            if a: res["whole"].append((np.sign(a[0]) == np.sign(G), a[0] * np.sign(G)))
            for hs, hz in ((ins & (x < 0), xw < 0), (ins & (x > 0), xw > 0)):
                if hs.sum() < 10: continue
                b = tc(hz)
                if b: res["half"].append((np.sign(b[0]) == np.sign(gam(hs)), b[0] * np.sign(gam(hs))))
        o = dict(run=run.replace(ROOT + "/", ""), k=k)
        for key, v in res.items():
            v = np.array(v) if v else np.zeros((0, 2))
            o[f"agree_{key}"] = v[:, 0].mean() if len(v) else np.nan; o[f"signed_{key}_deg"] = np.degrees(v[:, 1].mean()) if len(v) else np.nan; o[f"n_{key}"] = len(v)
        return o
    except Exception as e:
        return dict(run=run, error=str(e))
if __name__ == "__main__":
    runs = sorted(os.path.dirname(f) for sw in ("shapes_server", "shapes_mirror") for f in glob.glob(f"{ROOT}/{sw}/*/*/*/*/s0[1-3]/DONE"))
    with Pool(16) as p: res = p.map(one, runs, chunksize=2)
    man = {}
    for sw, mf in (("shapes_server", "manifest_shape_sweep.tsv"), ("shapes_mirror", "manifest_shape_mirror.tsv")):
        for r in csv.DictReader(open(f"{ROOT}/{sw}/{mf}"), delimiter="\t"): man[f"{sw}/{r['dir']}"] = r
    cols = ["run", "family", "aspect", "Req_la", "rho", "k", "agree_whole", "signed_whole_deg", "n_whole", "agree_half", "signed_half_deg", "n_half"]
    with open(H + "/tilt_sign.csv", "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=cols, extrasaction="ignore"); w.writeheader()
        for o in res:
            m = man.get(o["run"], {}); o.update(family=m.get("family"), aspect=m.get("aspect"), Req_la=m.get("Req_la"), rho=m.get("rho")); w.writerow(o)
    print(len(res), "runs, errors", sum("error" in o for o in res))
