#!/usr/bin/env python3
"""Cavity-level predictors of the wall-defect escape rate (675 runs). Per run: escape rate per wall defect (kinetics_carriers.csv g = r_esc/nwall), aspect, R_eq, k,
coherence / circulation / P(two) (shape_runs.csv), bulk +1/2 and -1/2 numbers (kinetics_carriers.csv), and flow in the wall band from flowfield (second half, frames every 100 t.u.):
rms speed, rms normal and tangential components, at depth 2-8 cells, and the rms speed in the interior (depth > 8). Linear regression of log(escape rate) on standardised
predictors, with and without the aspect indicator; R^2 and the remaining aspect coefficient. -> runlevel.txt"""
import csv, glob, gzip, json, os, sys
import numpy as np
from multiprocessing import Pool
from scipy.spatial import cKDTree
ROOT = "/home/baloewe/Physics/Overcharged/runs"; sys.path.insert(0, ROOT + "/cavity")
from wall_fourier import wall_polyline  # noqa
SW = ROOT + "/shapes_server"; TH = "/home/baloewe/Physics/Overcharged/reports/master/theory"; H = os.path.dirname(os.path.abspath(__file__))
def flowstats(rel):
    try:
        run = os.path.join(SW, rel); bc = json.load(open(run + "/input.json"))["BC"][0]; Q = bc["Q"]; W, s, P = wall_polyline(bc); tree = cKDTree(W)
        Wc = W[:, 0] + 1j * W[:, 1]; Tn = np.roll(Wc, -1) - np.roll(Wc, 1); Nn = -1j * Tn / np.abs(Tn)
        if np.mean((Nn * np.conj(Wc)).real) < 0: Nn = -Nn
        rows = np.loadtxt((l for l in gzip.open(run + "/flowfield.dat.gz", "rt") if len(l.split()) == 7 and l.split()[0][0].isdigit()))
        t = rows[:, 0]; m = t > t.max() / 2; x = rows[m, 1] + 0.5 - Q[0]; y = rows[m, 2] + 0.5 - Q[1]; vx, vy = rows[m, 4], rows[m, 5]
        import matplotlib.path as mp
        ins = mp.Path(W).contains_points(np.column_stack([x, y])); d, j = tree.query(np.column_stack([x, y])); nn = Nn[j]
        vn = vx * nn.real + vy * nn.imag; vt = -vx * nn.imag + vy * nn.real; band = ins & (d > 2) & (d < 8); core = ins & (d > 8)
        return rel, (np.sqrt(np.mean(vx[band] ** 2 + vy[band] ** 2)), np.sqrt(np.mean(vn[band] ** 2)), np.sqrt(np.mean(vt[band] ** 2)), np.sqrt(np.mean(vx[core] ** 2 + vy[core] ** 2)))
    except Exception:
        return rel, None
if __name__ == "__main__":
    kc = {r["run"]: r for r in csv.DictReader(open(TH + "/shapes/kinetics_carriers.csv"))}; sr = {r["run"]: r for r in csv.DictReader(open(TH + "/shapes/shape_runs.csv"))}
    runs = [r for r in kc if r in sr and kc[r].get("r_esc") not in (None, "", "nan")]
    with Pool(8) as pool: fs = dict(pool.map(flowstats, runs, chunksize=4))
    rows = []
    for r in runs:
        a, b, f = kc[r], sr[r], fs.get(r)
        if f is None: continue
        try:
            g = float(a["r_esc"]) / float(a["nwall"])
            if g <= 0: continue
            rows.append(dict(lg=np.log(g), elong=float(float(a["aspect"]) >= 1.25), aspect=float(a["aspect"]), Req=float(a["Req_la"]), k=float(a["k"]), coh=float(b["coh"]), circ=abs(float(b["circ"])),
                             two=float(b["two"]), nbp=float(a["nbp"]), nbm=float(a["nbm"]), vband=f[0], vn=f[1], vt=f[2], vcore=f[3]))
        except (ValueError, KeyError): continue
    keys = ["Req", "k", "coh", "circ", "two", "nbp", "nbm", "vband", "vn", "vt", "vcore"]; y = np.array([r["lg"] for r in rows]); L = []
    say = lambda x: (L.append(x), print(x, flush=True))
    def reg(cols):
        X = np.column_stack([np.ones(len(y))] + [(np.array([r[c] for r in rows]) - np.mean([r[c] for r in rows])) / (np.std([r[c] for r in rows]) or 1) if c != "elong" else np.array([r[c] for r in rows]) for c in cols])
        b, *_ = np.linalg.lstsq(X, y, rcond=None); res = y - X @ b; R2 = 1 - res.var() / y.var(); return b, R2
    b, R2 = reg(["elong"]); say(f"{len(y)} runs; log(escape rate) ~ elong: coefficient {b[1]:.2f} (x{np.exp(b[1]):.1f}), R^2 {R2:.2f}")
    for c in keys:
        b, R2 = reg(["elong", c]); say(f"  + {c:6s}: elong coefficient {b[1]:+.2f} (x{np.exp(b[1]):.1f}), {c} coefficient {b[2]:+.2f} per s.d., R^2 {R2:.2f}")
    b, R2 = reg(keys); say(f"  all predictors WITHOUT elong: R^2 {R2:.2f}; coefficients " + ", ".join(f"{c} {x:+.2f}" for c, x in zip(keys, b[1:])))
    b, R2 = reg(["elong"] + keys); say(f"  all predictors WITH elong: elong coefficient {b[1]:+.2f} (x{np.exp(b[1]):.1f}), R^2 {R2:.2f}")
    open(os.path.join(H, "runlevel.txt"), "w").write("\n".join(L) + "\n")
