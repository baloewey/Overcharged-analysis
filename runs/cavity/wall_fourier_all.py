#!/usr/bin/env python3
"""
Wall Fourier analysis for every circle run (main sweep, rho <= 0.15, +k and mirror -k), static and dynamic.

Wall-sign defects (+1/2 for +k, -1/2 for -k) within --band cells of the wall, at arclength s = R phi (ccw).
rho_m(n) = sum_j exp(i m phi_j) per defect frame (every 2 time units, first 20% dropped).
Static (as wall_fourier.py):
  S(m) = <|rho_m|^2>/<N>, m* = argmax over m >= 2; f_static = |<rho_m*>|^2 / <|rho_m*|^2>;
  m_lat = number of wall defects for a +1 bulk: 2k - 2 (+k), 2|k| + 2 (-k); ratio m*/m_lat
Dynamic at m*:
  F(t) = |<rho(n+t) rho*(n)>| / <|rho|^2>   lattice memory in the co-moving frame (sliding does not decorrelate it);
         F_1000 = F at t = 1000 time units (pinned or rigidly sliding lattice: ~1; melting lattice: -> 0);
         tau_F = first 1/e crossing (inf if none within 2000)
  tau_fl = 1/e time of the fluctuating part rho - <rho> (wall_fourier.py's tau)
  v_slide = lattice sliding speed along the wall: per flow window (100 time units) Omega = arg(sum rho(n+1) rho*(n)) / (m* dt),
            times R, signed with sign<v_theta> of that flow frame (> 0: the wall lattice slides with the circulation);
            v_ratio = v_slide / |<v_theta>|
Output analysis_fourier/wall_fourier_all.csv. Usage: wall_fourier_all.py [--workers 12] [--band 4]
"""
import argparse, csv, json, os, sys
from multiprocessing import Pool
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
OUT = os.path.join(HERE, "analysis_fourier")
COLS = ["set", "run", "k", "R", "rho", "Nwall", "mstar", "mlat", "mratio", "Sstar", "fstatic", "F1000", "tauF", "taufl",
        "v_slide", "v_abs", "v_ratio", "vth_abs"]
BAND = 4.0


def acorr(x, L):
    n = len(x); f = np.fft.fft(x, 2 * n)
    c = np.fft.ifft(f * np.conj(f))[:L] / (n - np.arange(L))
    return c                                  # c[l] = <x(n+l) x*(n)>


def one(run):
    from edge_criteria import read_with_angle
    try:
        full = os.path.join(HERE, "server_main", run)
        bc = json.load(open(os.path.join(full, "input.json")))["BC"][0]
        k = int(bc["customAnch"]); R = bc["R"]; sw = 1 if k > 0 else -1
        t, fr = read_with_angle(full, np.array(bc["Q"][:2], float), 0.2)
        dt = t[1] - t[0]; mmax = max(60, 2 * abs(k) + 10); m = np.arange(1, mmax + 1)
        rho = np.zeros((len(fr), mmax), complex); N = np.zeros(len(fr))
        for n, f in enumerate(fr):
            p = f[np.sign(f[:, 2]) == sw]
            if not len(p): continue
            r = np.hypot(p[:, 0], p[:, 1]); q = R - r < BAND
            ph = np.arctan2(p[q, 1], p[q, 0]); N[n] = q.sum()
            if q.any(): rho[n] = np.exp(1j * np.outer(ph, m)).sum(0)
        pw = (np.abs(rho) ** 2).mean(0); S = pw / max(N.mean(), 1e-9); mean = rho.mean(0)
        ms = 2 + int(np.argmax(S[1:])); i = ms - 1; x = rho[:, i]
        mlat = 2 * abs(k) - 2 if k > 0 else 2 * abs(k) + 2
        row = dict(set="pos" if k > 0 else "neg", run=run, k=k, R=R, rho=abs(k) / (np.pi * R), Nwall=N.mean(), mstar=ms, mlat=mlat,
                   mratio=ms / mlat if mlat > 0 else np.nan, Sstar=S[i], fstatic=abs(mean[i]) ** 2 / max(pw[i], 1e-12))
        L = min(len(x) // 2, int(2000 / dt) + 1)
        F = np.abs(acorr(x, L)) / max(pw[i], 1e-12)
        l1000 = int(round(1000 / dt)); row["F1000"] = F[l1000] if l1000 < L else np.nan
        b = np.where(F < np.exp(-1))[0]; row["tauF"] = b[0] * dt if len(b) else np.inf
        y = x - mean[i]; Cy = acorr(y, L).real / max(np.mean(np.abs(y) ** 2), 1e-12)
        b = np.where(Cy < np.exp(-1))[0]; row["taufl"] = b[0] * dt if len(b) else np.inf
        fl = np.load(os.path.join(HERE, "analysis_flow", "flow_red", run.replace("/", "__") + ".npz"))
        tf, vth = fl["t"], fl["vth"]; n0 = int(0.2 * len(tf))
        prod = x[1:] * np.conj(x[:-1]); tt = t[1:]
        vs, va = [], []
        for j in range(n0 + 1, len(tf)):
            w = (tt > tf[j - 1]) & (tt <= tf[j])
            if w.sum() < 5 or vth[j] == 0: continue
            om = np.angle(prod[w].sum()) / (ms * dt)
            vs.append(np.sign(vth[j]) * om * R); va.append(abs(vth[j]))
        row["v_slide"] = np.mean(vs) if vs else np.nan; row["vth_abs"] = np.mean(va) if va else np.nan
        row["v_abs"] = np.mean(np.abs(vs)) if vs else np.nan
        row["v_ratio"] = row["v_slide"] / row["vth_abs"] if vs and row["vth_abs"] > 0 else np.nan
        return row
    except Exception as ex:
        import traceback
        return dict(set="err", run=run, k=traceback.format_exc()[-300:])


def main():
    global BAND
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=12); ap.add_argument("--band", type=float, default=4.0)
    a = ap.parse_args(); BAND = a.band
    os.makedirs(OUT, exist_ok=True)
    runs = [l.strip() for l in open(os.path.join(HERE, "analysis_flow", "runlist.txt")) if l.strip()]
    with Pool(a.workers) as pool, open(os.path.join(OUT, "wall_fourier_all.csv"), "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=COLS, extrasaction="ignore"); w.writeheader()
        for i, row in enumerate(pool.imap_unordered(one, runs, chunksize=2)):
            if row["set"] == "err":
                print("FAILED", row, flush=True); continue
            w.writerow(row); fo.flush()
            if i % 100 == 0: print(i, flush=True)
    print("written")


if __name__ == "__main__":
    main()
