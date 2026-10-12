#!/usr/bin/env python3
"""
Wall Fourier, static + dynamic, for the oval (Cassini) and ellipse runs — the non-circular version of
runs/cavity/wall_fourier_all.py (same definitions, positions by arclength along the real wall).

+1/2 within --band cells of the wall, at the arclength s of the nearest wall point (ccw); P = perimeter.
rho_m(n) = sum_j exp(2 pi i m s_j / P) per defect frame (every 2 time units, first 20% dropped).
Static:  S(m) = <|rho_m|^2>/<N>, m* = argmax over m >= 2; f_static = |<rho_m*>|^2 / <|rho_m*|^2>;
         m_lat = 2k - 2 (wall count for a +1 bulk). Also S, f_static and F1000 at m = m_lat (the ordered period
         is not always the peak in the weakly ordered ovals).
Dynamic at m*: F(t) = |<rho(n+t) rho*(n)>| / <|rho|^2> (co-moving lattice memory), F1000 = F(1000 time units),
         tauF = first 1/e crossing; tau_fl = 1/e time of rho - <rho>.
Sliding: per flow window (flow frames every 100 time units) Omega = arg(sum rho(n+1) rho*(n)) / (m* dt),
         v_slide = Omega P / (2 pi) (arclength speed), signed with the sign of the near-wall tangential flow v_w
         of that window (> 0: lattice slides with the wall flow); v_ratio = v_slide / <|v_w|>.
Flow (from flowfield.dat.gz, cells fully inside): v_w = mean tangential velocity (ccw) of cells 1 < d < 5 from the wall;
         circ = mean (x vy - y vx)/r over all inside cells (single-vortex measure, as for circles); vrms.
         Note: for circles wall_fourier_all.py normalises by the interior |<v_theta>|, here by the near-wall flow.
Per-defect drift: wall +1/2 matched frame to frame (nearest neighbour within 1.5 cells), ds/dt along the wall, signed with
         sign(v_w) of the flow frame; drift = mean over 100-time-unit blocks, drift_se = standard error over blocks
         (robust where the lattice is weak and the m* phase is noise). circ_agree = fraction of flow frames with sign v_w = sign circ.
Output: wall_fourier_dyn.csv (one row per run). Usage: wall_fourier_dyn.py [--workers 12] [--band 4]
"""
import argparse, csv, glob, gzip, json, os, sys
from multiprocessing import Pool
import numpy as np
from matplotlib.path import Path
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "cavity"))
from edge_criteria import read_with_angle            # noqa: E402
from wall_fourier import wall_polyline, project      # noqa: E402

COLS = ["set", "run", "shape", "k", "aspect", "P", "Nwall", "mstar", "mlat", "Sstar", "fstatic", "F1000", "tauF", "taufl",
        "S_lat", "f_lat", "F1000_lat", "v_slide", "v_abs", "v_ratio", "vw_abs", "vw_signed_frac", "circ_abs", "circ_agree", "vrms",
        "drift", "drift_se", "drift_ratio"]
BAND = 4.0


def acorr(x, L):
    n = len(x); f = np.fft.fft(x, 2 * n)
    return np.fft.ifft(f * np.conj(f))[:L] / (n - np.arange(L))


def read_flow(run, Q):
    t, rows = [], []
    with gzip.open(os.path.join(run, "flowfield.dat.gz"), "rt") as f:
        for line in f:
            p = line.split()
            if len(p) == 7 and p[0][0] in "0123456789":
                rows.append([float(x) for x in p])
    A = np.array(rows)
    ts = np.unique(A[:, 0]); nc = len(A) // len(ts)
    A = A[: len(ts) * nc].reshape(len(ts), nc, 7)
    x = A[0, :, 1] + 0.5 - Q[0]; y = A[0, :, 2] + 0.5 - Q[1]
    return ts, x, y, A[:, :, 4], A[:, :, 5]


def aspect(bc):
    if bc.get("shape", 0) == 1:
        W, _, _ = wall_polyline(bc); return np.ptp(W[:, 0]) / np.ptp(W[:, 1])
    return bc["aInv"][0] / bc["aInv"][1]


def one(run):
    try:
        bc = json.load(open(os.path.join(run, "input.json")))["BC"][0]
        k = int(bc["customAnch"]); Q = np.array(bc["Q"][:2], float)
        W, s, P = wall_polyline(bc)
        t, fr = read_with_angle(run, Q, 0.2); dt = t[1] - t[0]
        mmax = max(60, 2 * k + 10); m = np.arange(1, mmax + 1)
        rho = np.zeros((len(fr), mmax), complex); N = np.zeros(len(fr))
        prev, dsv, dst = None, [], []
        for n, f in enumerate(fr):
            p = f[f[:, 2] > 0, :2]; cur = None
            if len(p):
                d, sa = project(W, s, p); q = d < BAND; sel = sa[q]; N[n] = len(sel); cur = (p[q], sel)
                if len(sel): rho[n] = np.exp(2j * np.pi * np.outer(sel, m) / P).sum(0)
            if prev is not None and cur is not None and len(prev[1]) and len(cur[1]):
                dd, jj = cKDTree(cur[0]).query(prev[0], distance_upper_bound=1.5); ok = np.isfinite(dd)
                dsv.append((cur[1][jj[ok]] - prev[1][ok] + P / 2) % P - P / 2); dst.append(np.full(ok.sum(), t[n]))
            prev = cur
        pw = (np.abs(rho) ** 2).mean(0); S = pw / max(N.mean(), 1e-9); mean = rho.mean(0)
        fst = np.abs(mean) ** 2 / np.maximum(pw, 1e-12)
        ms = 2 + int(np.argmax(S[1:])); i = ms - 1; x = rho[:, i]
        mlat = 2 * k - 2; il = mlat - 1
        rel = os.path.relpath(run, HERE)
        row = dict(set=rel.split("/")[0], run=rel, shape="oval" if bc.get("shape", 0) == 1 else "ellipse", k=k,
                   aspect=aspect(bc), P=P, Nwall=N.mean(), mstar=ms, mlat=mlat, Sstar=S[i], fstatic=fst[i])
        L = min(len(x) // 2, int(2000 / dt) + 1); l1000 = int(round(1000 / dt))
        F = np.abs(acorr(x, L)) / max(pw[i], 1e-12)
        row["F1000"] = F[l1000] if l1000 < L else np.nan
        b = np.where(F < np.exp(-1))[0]; row["tauF"] = b[0] * dt if len(b) else np.inf
        y = x - mean[i]; Cy = acorr(y, L).real / max(np.mean(np.abs(y) ** 2), 1e-12)
        b = np.where(Cy < np.exp(-1))[0]; row["taufl"] = b[0] * dt if len(b) else np.inf
        if il >= 1:
            Fl = np.abs(acorr(rho[:, il], L)) / max(pw[il], 1e-12)
            row.update(S_lat=S[il], f_lat=fst[il], F1000_lat=Fl[l1000] if l1000 < L else np.nan)
        # flow
        tf, cx, cy, VX, VY = read_flow(run, Q)
        inside = Path(W).contains_points(np.column_stack([cx, cy]))
        dw, sw = project(W, s, np.column_stack([cx, cy]))
        inside &= dw > 0.75                                  # drop cells cut by the wall
        band = inside & (dw > 1) & (dw < 5)
        jw = np.searchsorted(s, sw[band]); jw %= len(W)
        tang = W[(jw + 1) % len(W)] - W[(jw - 1) % len(W)]; tang /= np.hypot(*tang.T)[:, None]
        vw = (VX[:, band] * tang[:, 0] + VY[:, band] * tang[:, 1]).mean(1)
        r = np.hypot(cx[inside], cy[inside]); r[r < 1e-9] = 1e-9
        circ = ((cx[inside] * VY[:, inside] - cy[inside] * VX[:, inside]) / r).mean(1)
        vrms = np.sqrt((VX[:, inside] ** 2 + VY[:, inside] ** 2).mean(1))
        n0 = int(0.2 * len(tf))
        row.update(vw_abs=np.abs(vw[n0:]).mean(), circ_abs=np.abs(circ[n0:]).mean(), vrms=vrms[n0:].mean(),
                   circ_agree=np.mean(np.sign(vw[n0:]) == np.sign(circ[n0:])),
                   vw_signed_frac=abs(np.mean(np.sign(vw[n0:]))))   # 1 = wall flow never reverses
        prod = x[1:] * np.conj(x[:-1]); tt = t[1:]; vs = []
        for j in range(n0 + 1, len(tf)):
            w = (tt > tf[j - 1]) & (tt <= tf[j])
            if w.sum() < 5 or vw[j] == 0: continue
            vs.append(np.sign(vw[j]) * np.angle(prod[w].sum()) / (ms * dt) * P / (2 * np.pi))
        if dsv:
            ds = np.concatenate(dsv) / dt; tt = np.concatenate(dst)
            sg = np.sign(vw[np.clip(np.searchsorted(tf, tt), 0, len(tf) - 1)])
            blk = np.floor(tt / 100); bm = np.array([(ds * sg)[blk == u].mean() for u in np.unique(blk)])
            row.update(drift=bm.mean(), drift_se=bm.std() / np.sqrt(len(bm)), drift_ratio=bm.mean() / row["vw_abs"])
        if vs:
            row.update(v_slide=np.mean(vs), v_abs=np.mean(np.abs(vs)), v_ratio=np.mean(vs) / row["vw_abs"])
        return row
    except Exception:
        import traceback
        return dict(set="err", run=run, k=traceback.format_exc()[-400:])


def main():
    global BAND
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=12); ap.add_argument("--band", type=float, default=4.0)
    ap.add_argument("--only", default=None, help="substring filter on run paths (testing)")
    a = ap.parse_args(); BAND = a.band
    runs = sorted(os.path.dirname(p) for p in glob.glob(os.path.join(HERE, "*", "**", "DONE"), recursive=True))
    runs = [r for r in runs if "passive" not in r and "server_validation" not in r]
    if a.only: runs = [r for r in runs if a.only in r]
    print(len(runs), "runs", flush=True)
    out = os.path.join(HERE, "wall_fourier_dyn.csv" if not a.only else "wall_fourier_dyn_test.csv")
    with Pool(a.workers) as pool, open(out, "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=COLS, extrasaction="ignore"); w.writeheader()
        for i, row in enumerate(pool.imap_unordered(one, runs)):
            if row["set"] == "err":
                print("FAILED", row, flush=True); continue
            w.writerow(row); fo.flush(); print(i, row["run"], flush=True)
    print("written", out)


if __name__ == "__main__":
    main()
