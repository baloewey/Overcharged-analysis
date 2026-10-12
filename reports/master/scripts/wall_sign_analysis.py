#!/usr/bin/env python3
"""
Sign-aware wall-layer analysis for any cavity (master document; made for the negative-winding ellipses, 2026-10-05).

Wall sign = sign of the winding k; wall defects = defects of that sign within BAND = 4 cells of the wall polyline (nearest point).
Frames after the first 20%. Per run:
  N_wall          wall defects per frame; m_lat = 2k - 2 (k > 0) or 2|k| + 2 (k < 0), the wall count of a +1 bulk; occupancy N_wall/m_lat
  S_lat, f_lat    instantaneous order S(m_lat) = <|rho_m|^2>/<N> and static fraction |<rho_m>|^2/<|rho_m|^2> at m = m_lat,
                  rho_m = sum_j exp(2 pi i m s_j / P) over wall defects (s = arclength of the nearest wall point)
  tips/flanks     wall defects per 100 cells of wall with |x| > A/sqrt(2) vs the rest (A = largest |x| of the wall)
  travel          for pairs of frames 50 and 500 time units apart with the track in the band throughout: median along-wall |ds|
                  and fraction > half the site spacing P/(2 m_lat); tracks: build_tracks(rlink 8, gap 4)
  survival        Kaplan-Meier S(1000) and t_1/2 of wall visits of wall-sign tracks (wall s < 3, bulk s > 6, hysteresis;
                  visits with observed start; escape = event, track end or run end = censored); escape rate per wall defect
  charge          Q(c) from lobe_profiles.profile (levels = scaled copies of the wall for circles/ellipses, Cassini levels for ovals);
                  bulk charge = maximum of Q over c/b < 0.97 for k < 0 (as the circle mirror runs), Q at c/b = 0.85 for both signs
Usage (from runs/): wall_sign_analysis.py RUN [RUN ...] [--labels ...] [--workers 8] [--out FILE.tsv]
"""
import argparse, json, os, sys
from multiprocessing import Pool
import numpy as np
sys.path.insert(0, "cavity"); sys.path.insert(0, "cassini")
from edge_criteria import read_with_angle       # noqa: E402
from wall_fourier import wall_polyline, project  # noqa: E402
from defect_tracks import build_tracks          # noqa: E402
from lobe_profiles import profile               # noqa: E402

BAND, SW, SB = 4.0, 3.0, 6.0


def km(d, e):
    o = np.argsort(d); d, e = d[o], e[o]; ts, S, s, nr = [0.0], [1.0], 1.0, len(d)
    for tt in np.unique(d):
        m = d == tt; ne = e[m].sum()
        if ne: s *= 1 - ne / nr; ts.append(tt); S.append(s)
        nr -= m.sum()
    return np.array(ts), np.array(S)


def one(run):
    bc = json.load(open(os.path.join(run, "input.json")))["BC"][0]
    k = int(bc["customAnch"]); sg = 1 if k > 0 else -1; mlat = 2 * k - 2 if k > 0 else 2 * abs(k) + 2
    W, s, P = wall_polyline(bc); A = np.abs(W[:, 0]).max()
    t, fr = read_with_angle(run, np.array(bc["Q"][:2], float), 0.2); dt = t[1] - t[0]; nfr = len(fr)
    tipw = np.abs(W[:, 0]) > A / np.sqrt(2)
    seg = np.hypot(*np.diff(np.vstack([W, W[:1]]), axis=0).T); Lt, Lf = seg[tipw].sum(), seg[~tipw].sum()
    rho = np.zeros(nfr, complex); N = np.zeros(nfr); wt = wf = 0.0
    for n, f in enumerate(fr):
        p = f[np.sign(f[:, 2]) == sg, :2]
        if not len(p): continue
        d, sa = project(W, s, p); q = d < BAND; N[n] = q.sum()
        if q.any():
            rho[n] = np.exp(2j * np.pi * mlat * sa[q] / P).sum()
            tip = tipw[np.searchsorted(s, sa[q]) % len(W)]; wt += tip.sum(); wf += (~tip).sum()
    pw = np.mean(np.abs(rho) ** 2)
    out = dict(run=run, k=k, P=P, mlat=mlat, N_wall=N.mean(), occ=N.mean() / mlat, S_lat=pw / max(N.mean(), 1e-9),
               f_lat=abs(rho.mean()) ** 2 / max(pw, 1e-12), tips=wt / nfr / Lt * 100, flanks=wf / nfr / Lf * 100)
    out["tip_ratio"] = out["tips"] / max(out["flanks"], 1e-9)
    # tracks: travel and survival of wall-sign defects
    tracks, _ = build_tracks(fr, dt, 8.0, 4); half = P / (2 * mlat)
    for lag_t in (50, 500):
        lag = int(round(lag_t / dt)); ds = []
        for T in tracks.values():
            if np.sign(T["q"]) != sg or len(T["n"]) < lag + 1: continue
            d, sa = project(W, s, T["xy"]); inw = d < BAND
            su = np.unwrap(sa / P * 2 * np.pi) * P / (2 * np.pi); idx = {v: i for i, v in enumerate(T["n"])}
            for i, ni in enumerate(T["n"]):
                j = idx.get(ni + lag)
                if j is not None and inw[i:j + 1].all(): ds.append(abs(su[j] - su[i]))
        ds = np.array(ds)
        out[f"pairs{lag_t}"] = len(ds); out[f"med{lag_t}"] = np.median(ds) if len(ds) else np.nan
        out[f"far{lag_t}"] = np.mean(ds > half) if len(ds) else np.nan
    dur, ev, wocc, nesc = [], [], 0, 0
    for T in tracks.values():
        if np.sign(T["q"]) != sg: continue
        d, _ = project(W, s, T["xy"]); st = np.full(len(d), -1); cur = -1
        for i, di in enumerate(d):
            if di < SW: cur = 1
            elif di > SB: cur = 0
            st[i] = cur
        wocc += np.sum(st == 1); i = 0; n = T["n"]
        while i < len(st):
            if st[i] < 0: i += 1; continue
            j = i
            while j + 1 < len(st) and st[j + 1] == st[i]: j += 1
            nxt = st[j + 1] if j + 1 < len(st) else None
            if st[i] == 1:
                if nxt == 0: nesc += 1
                if i > 0 and st[i - 1] >= 0: dur.append((n[j] - n[i] + 1) * dt); ev.append(1 if nxt == 0 else 0)
            i = j + 1
    ts, S = km(np.array(dur), np.array(ev)) if dur else (np.array([0.0]), np.array([1.0]))
    out["visits"] = len(dur); out["S1000"] = S[np.searchsorted(ts, 1000, side="right") - 1]
    out["t_half"] = ts[np.argmax(S <= 0.5)] if (S <= 0.5).any() else np.inf
    out["esc_rate"] = nesc / max(wocc / nfr, 1e-9) / (nfr * dt)
    # charge profile
    p = profile(run, 200, 0.2); x = p["edges"][1:]; Q = p["Q"]
    out["Q085"] = float(np.interp(0.85, x, Q)); out["Qpeak"] = float(Q[x < 0.97].max()) if k < 0 else np.nan
    return out


COLS = ["label", "k", "P", "mlat", "N_wall", "occ", "S_lat", "f_lat", "tips", "flanks", "tip_ratio", "pairs50", "med50", "far50",
        "pairs500", "med500", "far500", "visits", "S1000", "t_half", "esc_rate", "Q085", "Qpeak"]


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("runs", nargs="+"); ap.add_argument("--labels", nargs="+")
    ap.add_argument("--workers", type=int, default=8); ap.add_argument("--out", default=None)
    a = ap.parse_args(); labels = a.labels or a.runs
    with Pool(a.workers) as pool: res = pool.map(one, a.runs)
    lines = ["\t".join(COLS)]
    for lab, r in zip(labels, res):
        r["label"] = lab
        lines.append("\t".join(f"{r[c]:.4g}" if isinstance(r[c], float) else str(r[c]) for c in COLS))
    txt = "\n".join(lines); print(txt)
    if a.out: open(a.out, "w").write(txt + "\n")


if __name__ == "__main__":
    main()
