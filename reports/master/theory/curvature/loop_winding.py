#!/usr/bin/env python3
"""Direct winding test of the topological form of the 2k-2 selection (user, 2026-11-23; RESEARCH 'TOPOLOGICAL FORM'). Claim: on a closed loop where
the director keeps an angle alpha to the loop tangent that never reaches +-90 deg, the director winding along the loop is the loop's turning number +1,
so the enclosed charge is exactly +1. Self-contained (server copy, ~/Overcharged/analysis_curv): circles (cavity/main_act0.07) and the shape sweep
(shapes/circle, ellipse, oval; eMode 2 / Cassini), director every 100 t.u., second half of each run.
Per director frame and loop (offset curve at depth y = 3 and 5 cells inside the wall; nearest-cell director):
  alpha_j = director angle - tangent angle (mod pi, in (-90, 90]); w = winding of the director along the loop (unwrap 2 theta) / 2;
  Q_in = net charge of the defects inside the loop (defects.dat at the same time); aligned = max |alpha| < 75 deg along the loop.
Per run: frames, P(w = Q_in) (estimator / topology check), P(aligned), P(Q_in = 1 | aligned), P(Q_in = 1 | not aligned), <Q_in>.
-> loop_winding.csv (per run), loop_winding.txt (by shape class)"""
import csv, glob, gzip, json, os, sys, numpy as np
from multiprocessing import Pool
from matplotlib.path import Path
H = os.path.dirname(os.path.abspath(__file__)); RUNS = os.path.expanduser("~/Overcharged/runs"); YS = (3.0, 5.0); AMAX = np.radians(75)


def wall(bc, n=4000):
    th = np.linspace(0, 2 * np.pi, n, endpoint=False)
    if bc.get("shape", 0) == 1:
        a, b = bc.get("cassiniA", 0.0), bc["R"]; r2 = a * a * np.cos(2 * th) + np.sqrt(b ** 4 - a ** 4 * np.sin(2 * th) ** 2)
        x, y = np.sqrt(r2) * np.cos(th), np.sqrt(r2) * np.sin(th)
    else:
        A, B = bc.get("aInv", [1, 1])[0] * bc["R"], bc.get("aInv", [1, 1])[1] * bc["R"]; x, y = A * np.cos(th), B * np.sin(th)
    W = np.column_stack([x, y]); t = np.roll(W, -1, 0) - np.roll(W, 1, 0); t /= np.hypot(*t.T)[:, None]   # counterclockwise tangent
    nin = np.column_stack([-t[:, 1], t[:, 0]]); nin *= np.sign(np.sum(nin * (-W), 1))[:, None]
    return W, t, nin


def frames_dir(fn):
    with gzip.open(fn, "rb") as f: raw = f.read()
    lines = raw.split(b"\n"); i = next(j for j, l in enumerate(lines) if len(l.split()) == 8 and l.split()[0][:1].isdigit())
    a = np.array(b" ".join(lines[i:]).split(), dtype=np.float32)
    a = a[:len(a) // 8 * 8].reshape(-1, 8); return a


def frames_def(fn):
    out, cur, tcur, nexp = {}, [], None, 0
    with gzip.open(fn, "rt") as f:
        for line in f:
            p = line.split()
            if len(p) == 2:
                try: t, n = float(p[0]), int(p[1])
                except ValueError: continue
                if tcur is not None and len(cur) == nexp: out[round(tcur, 1)] = np.array(cur).reshape(-1, 4)
                cur, tcur, nexp = [], t, n
            elif len(p) == 4 and tcur is not None: cur.append([float(v) for v in p])
    if tcur is not None and len(cur) == nexp: out[round(tcur, 1)] = np.array(cur).reshape(-1, 4)
    return out


def one(run):
    try:
        J = json.load(open(run + "/input.json")); bc = J["BC"][0]; Q = np.array(bc["Q"][:2], float); k = bc["customAnch"]; nx_, ny_ = J["domain"][:2]
        W, T, N = wall(bc); A = frames_dir(run + "/directorfield.dat.gz"); F = frames_def(run + "/defects.dat.gz")
        ts = np.unique(A[:, 0]); ts = ts[ts >= ts[0] + 0.5 * (ts[-1] - ts[0])]
        loops = {}
        for y in YS:
            P = W + y * N; step = max(1, int(len(P) / (2 * np.sum(np.hypot(*np.diff(P, axis=0).T)))))   # ~2 points per cell along the loop
            idx = np.arange(0, len(P), step); loops[y] = (P[idx], T[idx], Path(P))
        st = {y: dict(n=0, eq=0, al=0, q1_al=0, q1_nal=0, nal=0, qsum=0.0) for y in YS}
        for tt in ts:
            d = A[A[:, 0] == tt]; G = np.full((nx_, ny_, 3), np.nan, np.float32)
            G[d[:, 1].astype(int), d[:, 2].astype(int)] = d[:, [4, 5, 7]]
            f = F.get(round(float(tt), 1))
            if f is None: continue
            for y in YS:
                P, Tg, path = loops[y]; ci = np.clip(np.floor(P + Q).astype(int), 0, [nx_ - 1, ny_ - 1])
                g = G[ci[:, 0], ci[:, 1]]; ok = np.isfinite(g[:, 0]) & (g[:, 2] > 0.05)
                if ok.mean() < 0.9: continue
                th = np.arctan2(g[ok, 1], g[ok, 0]); ps = np.arctan2(Tg[ok, 1], Tg[ok, 0])
                al = (th - ps + np.pi / 2) % np.pi - np.pi / 2
                w = np.sum(np.angle(np.exp(1j * np.diff(np.r_[2 * th, 2 * th[0]])))) / (4 * np.pi)
                qin = f[path.contains_points(f[:, :2] - Q), 2].sum() if len(f) else 0.0          # defects.dat positions are absolute
                if os.environ.get("DEBUG") and y == 3.0: print(f"t {tt:.0f} w {w:+.2f} Qin {qin:+.1f} max|alpha| {np.degrees(np.max(np.abs(al))):.0f} ok {ok.mean():.2f}")
                s = st[y]; s["n"] += 1; s["eq"] += abs(w - qin) < 0.25; s["qsum"] += qin; aligned = np.max(np.abs(al)) < AMAX
                if aligned: s["al"] += 1; s["q1_al"] += abs(qin - 1) < 0.25
                else: s["nal"] += 1; s["q1_nal"] += abs(qin - 1) < 0.25
        kind = run.split("/runs/")[1].split("/")[0] if "/shapes/" not in run else "shapes/" + run.split("/shapes/")[1].split("/")[0]
        geo = run.split("/runs/")[1]
        row = dict(run=geo, kind=kind, k=k, R=bc["R"])
        for y in YS:
            s = st[y]; n = max(s["n"], 1)
            row.update({f"n{y:g}": s["n"], f"Peq{y:g}": s["eq"] / n, f"Pal{y:g}": s["al"] / n, f"Q1al{y:g}": s["q1_al"] / max(s["al"], 1) if s["al"] else np.nan,
                        f"Q1nal{y:g}": s["q1_nal"] / max(s["nal"], 1) if s["nal"] else np.nan, f"Qin{y:g}": s["qsum"] / n})
        return row
    except Exception as ex:
        return dict(run=run, err=str(ex)[:120])


if __name__ == "__main__":
    runs = sorted(os.path.dirname(p) for p in glob.glob(RUNS + "/shapes/*/*/*/k*/s*/directorfield.dat.gz") + glob.glob(RUNS + "/shapes/circle/circle/*/k*/s*/directorfield.dat.gz"))
    runs += sorted(os.path.dirname(p) for p in glob.glob(RUNS + "/cavity/main_act0.07/rho*/k*/s[01][0-9]/directorfield.dat.gz"))
    runs = sorted(set(runs))
    with Pool(int(os.environ.get("NPROC", 24))) as p: res = p.map(one, runs, chunksize=4)
    err = [r for r in res if "err" in r]; res = [r for r in res if "err" not in r]
    with open(os.path.join(H, "loop_winding.csv"), "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=list(res[0])); w.writeheader(); w.writerows(res)
    print("runs", len(res), "errors", len(err), err[:3])
