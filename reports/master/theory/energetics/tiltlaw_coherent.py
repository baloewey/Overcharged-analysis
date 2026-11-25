#!/usr/bin/env python3
"""Tilt law across lambda, three readings of the wall-defect tilt (user, 2026-11-25; runs/lambda_circle, second half of each run, same wall selection as tiltlaw.py):
 (1) <|tau|> over all frames (as Fig. 3b);
 (2) <|tau|> over COHERENT frames only (per-frame C = |<exp(i tau_j)>_j| > 0.9), and over coherent frames that also hold the full 2k-2 ring;
 (3) |<tau>|: per frame the absolute value of the SIGNED mean tilt, |(1/N) sum_j tau_j|, then averaged over frames (a common tilt survives, a +/- mixture cancels).
Pinning check: along-wall speed of the wall +1/2 (nearest-neighbour matching between consecutive defect frames, matches within 1.5 cells, wall defects only),
 signed with the frame's majority tilt (+ = sliding in the sense set by the tilt), and the fraction of wall defects with |tau| < 15 deg (sitting at a homeotropic site).
-> tiltlaw_coherent.txt, tiltlaw_coherent.npy"""
import os, sys, json, glob, numpy as np
from multiprocessing import Pool
from scipy.spatial import cKDTree
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
import energetics as En  # noqa
ROOT = "/home/baloewe/Physics/Overcharged"


def one(run):
    J = json.load(open(run + "/input.json")); bc = J["BC"][0]; Q = bc["Q"]; k = bc["customAnch"]; lam = J["species"][0].get("tumble", 2.0); N = 2 * k - 2
    W, s, P = En.wall_polyline(bc); tree = cKDTree(W); Wc = W[:, 0] + 1j * W[:, 1]; Tn = np.roll(Wc, -1) - np.roll(Wc, 1); Tn /= np.abs(Tn); Nn = -1j * Tn
    if np.mean((Nn * np.conj(Wc)).real) < 0: Nn = -Nn
    t, fr = En.read_with_angle(run, np.array(Q[:2], float), 0.5)
    absall, abscoh, absfull, signed, C, full, small, v_signed, v_abs = [], [], [], [], [], [], [], [], []
    prev = None
    for i, f in enumerate(fr):
        if not len(f): prev = None; continue
        d, j = tree.query(f[:, :2]); w = (f[:, 2] > 0) & (d < 4)
        if not w.any(): prev = None; continue
        tau = np.angle(np.exp(1j * (f[w, 3] - np.angle(Nn[j[w]])))); c = abs(np.mean(np.exp(1j * tau))); isfull = w.sum() == N
        absall += list(np.abs(tau)); C.append(c); full.append(isfull); signed.append(abs(np.mean(tau))); small += list(np.abs(tau) < np.radians(15))
        if c > 0.9:
            abscoh += list(np.abs(tau))
            if isfull: absfull += list(np.abs(tau))
        z = f[w, 0] + 1j * f[w, 1]; sg = np.sign(np.mean(tau)) or 1.0
        if prev is not None and t[i] - prev[0] > 0:
            zp = prev[1]; dd, jj = cKDTree(np.c_[zp.real, zp.imag]).query(np.c_[z.real, z.imag]); m = dd < 1.5
            if m.any():
                dphi = np.angle(z[m] / zp[jj[m]]); v = np.abs(z[m]) * dphi / (t[i] - prev[0])   # along-wall speed (+ = counter-clockwise)
                v_signed += list(sg * v); v_abs += list(np.abs(v))
        prev = (t[i], z)
    deg = lambda x: np.degrees(np.mean(x)) if len(x) else np.nan
    return dict(k=k, lam=lam, run=os.path.basename(run), all=deg(absall), coh=deg(abscoh), cohfull=deg(absfull), signed=deg(signed),
                fcoh=np.mean(np.array(C) > 0.9), fcohfull=np.mean((np.array(C) > 0.9) & np.array(full)), C=np.mean(C), small=np.mean(small),
                vs=np.mean(v_signed) if v_signed else np.nan, va=np.mean(v_abs) if v_abs else np.nan, dt=np.median(np.diff(t)))


if __name__ == "__main__":
    runs = sorted(r for r in glob.glob(ROOT + "/runs/lambda_circle/k*_lam*_s*") if os.path.exists(r + "/DONE"))
    with Pool(int(os.environ.get("NPROC", 12))) as p: res = p.map(one, runs, chunksize=1)
    G = {}
    for r in res: G.setdefault((r["k"], r["lam"]), []).append(r)
    keys = ["all", "coh", "cohfull", "signed", "fcoh", "fcohfull", "C", "small", "vs", "va"]
    out = {}; L = ["TILT LAW, three readings (deg; mean +- s.e. over seeds). coh = frames with C > 0.9; coh+full = also exactly 2k-2 wall +1/2; |<tau>| = per-frame |signed mean|.",
                   "Pinning: along-wall defect speed (cells/t.u.), signed with the frame's majority tilt, and |v|; fraction of wall defects with |tau| < 15 deg.",
                   "k lam | runs | arccos | <|tau|> all | <|tau|> coh | <|tau|> coh+full | |<tau>| | frac coh | frac coh+full | C | v signed | |v| | frac |tau|<15"]
    for key in sorted(G):
        g = G[key]; A = np.array([[r[x] for x in keys] for r in g], float); m = np.nanmean(A, 0); e = np.nanstd(A, 0) / np.sqrt(max(len(g) - 1, 1))
        ac = np.degrees(np.arccos(1 / key[1])) if key[1] >= 1 else np.nan; out[key] = dict(zip(keys, m)); out[key].update({x + "_se": v for x, v in zip(keys, e)}); out[key]["arccos"] = ac
        L.append(f"{key[0]} {key[1]:3.1f} | {len(g)} | {ac:5.1f} | {m[0]:5.1f}±{e[0]:.1f} | {m[1]:5.1f}±{e[1]:.1f} | {m[2]:5.1f}±{e[2]:.1f} | {m[3]:5.1f}±{e[3]:.1f} | "
                 f"{m[4]:.2f} | {m[5]:.2f} | {m[6]:.2f} | {m[8]:+.4f}±{e[8]:.4f} | {m[9]:.4f} | {m[7]:.3f}")
    L.append(f"defect frame interval (median): {np.median([r['dt'] for r in res]):g} t.u.")
    np.save(os.path.join(H, "tiltlaw_coherent.npy"), out, allow_pickle=True)
    open(os.path.join(H, "tiltlaw_coherent.txt"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
