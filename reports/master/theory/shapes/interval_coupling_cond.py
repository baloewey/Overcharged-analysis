#!/usr/bin/env python3
"""Does the wall layer follow the instantaneous flow? Elongated shapes of the server sweep (aspect >= 1.25, ellipses and ovals, all R_eq, k, seeds).
Per flow interval (100 t.u.; second half of each run): normalised circulation c = |<v_theta>| / v_rms, two-vortex indicator, wall +1/2 count n_w (mean over the
interval's defect frames), tilt coherence of the wall layer. Within each run the series are z-scored; pooled over runs:
  (1) n_w in one- vs two-vortex intervals (z-score and raw, relative to 2k - 2);
  (2) lagged correlation corr(c(i), n_w(i + lag)) and corr(coherence(i), n_w(i + lag)), lag = -5..+5 intervals (+: wall count after the flow state);
  (3) the same for the circle-like shapes (aspect <= 1.22) as a control. Writes interval_coupling.txt."""
import csv, glob, gzip, json, os, sys
import numpy as np
from multiprocessing import Pool
ROOT = "/home/baloewe/Physics/Overcharged/runs"; sys.path.insert(0, ROOT + "/cavity")
from edge_criteria import read_with_angle  # noqa
from wall_fourier import wall_polyline  # noqa
from scipy.spatial import cKDTree
import matplotlib.path as mpath
H = os.path.dirname(os.path.abspath(__file__)); SW = ROOT + "/shapes_server"
def one(rel):
    run = os.path.join(SW, rel)
    try:
        bc = json.load(open(run + "/input.json"))["BC"][0]; k = bc["customAnch"]; Q = bc["Q"]
        W, s, P = wall_polyline(bc); Wc = W[:, 0] + 1j * W[:, 1]; tree = cKDTree(W)
        Tn = np.roll(Wc, -1) - np.roll(Wc, 1); Nn = -1j * Tn / np.abs(Tn); path = mpath.Path(W)
        rows = np.loadtxt((l for l in gzip.open(run + "/flowfield.dat.gz", "rt") if len(l.split()) == 7 and l.split()[0][0].isdigit()))
        ts = np.unique(rows[:, 0]); ts = ts[ts > ts.max() / 2]
        x0 = rows[:, 1] + 0.5 - Q[0]; y0 = rows[:, 2] + 0.5 - Q[1]
        t, fr = read_with_angle(run, np.array(Q[:2], float), 0.45); t = np.array(t)
        C, TW, NW, CO = [], [], [], []
        for tt in ts:
            m = rows[:, 0] == tt; x, y, vx, vy = x0[m], y0[m], rows[m, 4], rows[m, 5]; ins = path.contains_points(np.column_stack([x, y])); r = np.hypot(x, y)
            vt = (x * vy - y * vx) / np.maximum(r, 0.5); vrms = np.sqrt(np.mean((vx ** 2 + vy ** 2)[ins]))
            def hr(sel):
                cx, cy = x[sel].mean(), y[sel].mean(); dx, dy = x[sel] - cx, y[sel] - cy
                return np.mean((dx * vy[sel] - dy * vx[sel]) / np.maximum(np.hypot(dx, dy), 0.5))
            sel = (t > tt - 100) & (t <= tt); idx = np.flatnonzero(sel)[::5]
            nws, cos = [], []
            for j in idx:
                f = fr[j]; p = f[f[:, 2] > 0]
                if not len(p): nws.append(0); continue
                d, i = tree.query(p[:, :2]); w = d < 4; nws.append(w.sum())
                if w.any(): cos.append(np.mean(np.exp(1j * (p[w, 3] - np.angle(Nn[i[w]])))))
            if not nws: continue
            C.append(abs(np.mean(vt[ins])) / max(vrms, 1e-9)); TW.append(int(hr(ins & (x < 0)) * hr(ins & (x > 0)) < 0)); NW.append(np.mean(nws))
            CO.append(abs(np.mean(cos)) if cos else np.nan)
        return dict(rel=rel, k=k, C=np.array(C), TW=np.array(TW), NW=np.array(NW), CO=np.array(CO))
    except Exception as ex:
        return dict(rel=rel, error=str(ex))
if __name__ == "__main__":
    rows = list(csv.DictReader(open(os.path.join(H, "shape_runs.csv"))))
    elong = [r["run"] for r in rows if r["family"] in ("ellipse", "oval") and float(r["aspect"]) >= 1.25]
    ctrl = [r["run"] for r in rows if r["family"] in ("ellipse", "oval", "circle") and float(r["aspect"]) <= 1.22]
    import pickle
    with Pool(6) as pool: RE = pool.map(one, elong, chunksize=2); RC = pool.map(one, ctrl, chunksize=2)
    pickle.dump(dict(RE=RE, RC=RC), open(os.path.join(H, "interval_coupling_cond.pkl"), "wb"))
    L = []
    def say(s): print(s, flush=True); L.append(s)
    def z(a): a = np.asarray(a, float); return (a - np.nanmean(a)) / (np.nanstd(a) if np.nanstd(a) > 0 else 1)
    for lab, RS in (("elongated (aspect >= 1.25)", RE), ("circle-like control (aspect <= 1.22)", RC)):
        RS = [r for r in RS if "error" not in r and len(r["NW"]) > 10]
        say(f"\n== {lab}: {len(RS)} runs, {sum(len(r['NW']) for r in RS)} intervals")
        one_z = np.concatenate([z(r["NW"])[r["TW"] == 0] for r in RS]); two_z = np.concatenate([z(r["NW"])[r["TW"] == 1] for r in RS])
        one_v = np.concatenate([(2 * r["k"] - 2) - r["NW"][r["TW"] == 0] for r in RS]); two_v = np.concatenate([(2 * r["k"] - 2) - r["NW"][r["TW"] == 1] for r in RS])
        say(f"  wall count z-score: one-vortex intervals {np.mean(one_z):+.3f} (n {len(one_z)}), two-vortex {np.mean(two_z) if len(two_z) else np.nan:+.3f} (n {len(two_z)}); "
            f"vacancies: one {np.mean(one_v):.2f}, two {np.mean(two_v) if len(two_v) else np.nan:.2f}")
        for nm, key in (("circulation c", "C"), ("tilt coherence", "CO")):
            line = f"  corr({nm}(i), wall count(i+lag)): "
            for lag in range(-5, 6):
                xs, ys = [], []
                for r in RS:
                    a, b = z(r[key]), z(r["NW"]); n = len(a)
                    for i in range(n):
                        if 0 <= i + lag < n and np.isfinite(a[i]) and np.isfinite(b[i + lag]): xs.append(a[i]); ys.append(b[i + lag])
                line += f"{lag:+d}:{np.corrcoef(xs, ys)[0, 1]:+.2f} "
            say(line + "(lags in 100 t.u.)")
    # conditioning on circle-like coherence: wall count (vacancies) in intervals binned by tilt coherence (and the next interval)
    for lab, RS in (("elongated", RE), ("circle-like", RC)):
        RS = [r for r in RS if "error" not in r and len(r["NW"]) > 10]
        co = np.concatenate([r["CO"] for r in RS]); vac = np.concatenate([(2 * r["k"] - 2) - r["NW"] for r in RS])
        vac_next = np.concatenate([np.r_[(2 * r["k"] - 2) - r["NW"][1:], np.nan] for r in RS])
        for lo, hi in ((0, 0.5), (0.5, 0.7), (0.7, 0.85), (0.85, 0.95), (0.95, 1.01)):
            m = (co >= lo) & (co < hi)
            say(f"  [{lab}] tilt coherence {lo:.2f}-{hi:.2f}: {m.sum():6d} intervals ({m.mean():.1%}); vacancies {np.nanmean(vac[m]) if m.any() else np.nan:.2f}, next interval {np.nanmean(vac_next[m]) if m.any() else np.nan:.2f}")
    # same binning per R_eq (size confound), and within-run: vacancies minus the run mean (removes shape/size/k entirely)
    Req = {r["run"]: r["Req_la"] for r in rows}
    for lab, RS in (("elongated", RE), ("circle-like", RC)):
        RS = [r for r in RS if "error" not in r and len(r["NW"]) > 10]
        for R in sorted(set(Req[r["rel"]] for r in RS)):
            S = [r for r in RS if Req[r["rel"]] == R]
            co = np.concatenate([r["CO"] for r in S]); vac = np.concatenate([(2 * r["k"] - 2) - r["NW"] for r in S])
            say(f"  [{lab} R_eq {R}] " + "  ".join(f"coh {lo:.2f}-{hi:.2f}: vac {np.nanmean(vac[(co >= lo) & (co < hi)]):.2f} ({((co >= lo) & (co < hi)).mean():.0%})"
                                                for lo, hi in ((0, 0.5), (0.5, 0.7), (0.7, 0.85), (0.85, 1.01)) if ((co >= lo) & (co < hi)).any()))
        co = np.concatenate([r["CO"] for r in RS]); dv = np.concatenate([-(r["NW"] - np.nanmean(r["NW"])) for r in RS])
        say(f"  [{lab} within-run] " + "  ".join(f"coh {lo:.2f}-{hi:.2f}: vac - run mean {np.nanmean(dv[(co >= lo) & (co < hi)]):+.2f}"
                                               for lo, hi in ((0, 0.5), (0.5, 0.7), (0.7, 0.85), (0.85, 1.01)) if ((co >= lo) & (co < hi)).any()))
    open(os.path.join(H, "interval_coupling_conditioned.txt"), "w").write("\n".join(L) + "\n")
