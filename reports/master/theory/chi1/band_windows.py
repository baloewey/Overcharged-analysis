#!/usr/bin/env python3
"""Band angle read (i) in the fixed window 2-3 cells and (ii) JUST OUTSIDE each layer's near field, window [edge, edge + 1) with edge = d + a/2pi from the
frame's own median wall-defect depth d and spacing a = 2 pi (R - d)/(2k - 2) (user, 2026-11-24: the field the layer sees beyond its own structure).
Both static-subtracted (static_band.py: 2 alpha_corr = 2 alpha_meas - 2 alpha_static + tau_static). Coherent frames (C > 0.9) with exactly 2k-2 wall +1/2.
Per run: tilt, d, a, edge, and for each window 2 alpha_meas, 2 alpha_corr. Sets: MPCD chi 0.5 (R 31.8) and chi 1 (R 21.6).
Band thickness l_b per k: depth where the ring-averaged alpha(y) profile (curvature_profiles.npz / chi1_band.npz) falls to half its maximum over y 1-4.
-> band_windows.npy (per k means), band_windows.txt"""
import glob, gzip, json, os, sys, warnings, numpy as np
from multiprocessing import Pool
warnings.filterwarnings("ignore")
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
import static_band as SB  # noqa  (TS, read_with_angle)
ROOT = SB.ROOT


def vals(X, Y, R, th2, S, th2s, sg, lo, hi, ts_):
    rr = np.hypot(X, Y); ph = np.arctan2(Y, X); u = R - rr; b = (u >= lo) & (u < hi) & np.isfinite(th2)
    if b.sum() < 10: return np.nan, np.nan
    am = np.angle(np.sum(S[b] * np.exp(1j * sg * (th2[b] - 2 * ph[b] - np.pi)))) / 2
    ast = np.angle(np.sum(np.exp(1j * sg * (th2s[b] - 2 * ph[b] - np.pi)))) / 2
    return 2 * am, 2 * am - 2 * ast + ts_


def run_one(args):
    run, stride = args
    try:
        bc = json.load(open(run + "/input.json"))["BC"][0]; R = bc["R"]; k = bc["customAnch"]; Q = np.array(bc["Q"][:2], float); N = 2 * k - 2
        with gzip.open(run + "/directorfield.dat.gz", "rt") as fh: D = np.loadtxt((l for l in fh if len(l.split()) == 8 and l.split()[0][0].isdigit()))
        td = np.unique(D[:, 0]); td = td[td >= td[0] + 0.5 * (td[-1] - td[0])][::stride]; t, fr = SB.read_with_angle(run, Q, 0.0); t = np.array(t)
        n = int(D[:, 1].max()) + 1; m = int(D[:, 2].max()) + 1; X, Y = np.meshgrid(np.arange(n) + 0.5 - Q[0], np.arange(m) + 0.5 - Q[1], indexing="ij")
        rows = []
        for tt in td:
            i = int(np.argmin(np.abs(t - tt)))
            if abs(t[i] - tt) > 1e-6: continue
            f = fr[i]; zs = f[:, 0] + 1j * f[:, 1]; qs = f[:, 2]; r = np.abs(zs); w = (qs > 0) & (R - r < 4)
            if w.sum() != N or abs(qs.sum() - k) > 1e-6: continue
            ez = np.mean(np.exp(1j * (f[w, 3] - np.angle(zs[w]))))
            if abs(ez) < 0.9: continue
            sg = np.sign(np.angle(ez)) or 1.0; d = np.median(R - r[w]); a = 2 * np.pi * (R - d) / N; edge = d + a / (2 * np.pi)
            dd = D[D[:, 0] == tt]; th2 = np.full((n, m), np.nan); S = np.zeros((n, m)); ix, iy = dd[:, 1].astype(int), dd[:, 2].astype(int)
            th2[ix, iy] = 2 * np.arctan2(dd[:, 5], dd[:, 4]); S[ix, iy] = dd[:, 7]; th2[np.hypot(X, Y) >= R] = np.nan
            th2s, _, _ = SB.TS.fields(X, Y, list(zs), list(qs), R)
            nt = SB.TS.natural_tilts(list(zs), list(qs), R, len(zs)); ts_ = sg * np.angle(np.mean(np.exp(1j * nt[np.flatnonzero(w)])))
            fm, fc = vals(X, Y, R, th2, S, th2s, sg, 2.0, 3.0, ts_); em, ec = vals(X, Y, R, th2, S, th2s, sg, edge, edge + 1.0, ts_)
            rows.append((abs(np.angle(ez)), d, a, edge, fm, fc, em, ec))
        if len(rows) < 3: return None
        A = np.array(rows); deg = lambda c: np.degrees(np.nanmedian(A[:, c]))
        return dict(run=run, k=k, R=R, n=len(rows), tau=deg(0), d=np.median(A[:, 1]), a=np.median(A[:, 2]), edge=np.median(A[:, 3]),
                    fm=deg(4), fc=deg(5), em=deg(6), ec=deg(7))
    except Exception as ex:
        return dict(run=run, err=repr(ex)[:150])


def band_thickness(prof_y, prof_a):
    """depth where alpha(y) falls to half of its maximum over y 1-4 (linear interpolation beyond the maximum)"""
    m = (prof_y >= 1) & (prof_y <= 4); amax = np.nanmax(prof_a[m]); j0 = int(np.flatnonzero(m)[np.nanargmax(prof_a[m])])
    for j in range(j0, len(prof_y) - 1):
        if prof_a[j + 1] < amax / 2: return prof_y[j] + (prof_y[j + 1] - prof_y[j]) * (prof_a[j] - amax / 2) / max(prof_a[j] - prof_a[j + 1], 1e-9)
    return np.nan


if __name__ == "__main__":
    c05 = [(os.path.dirname(p), 1) for p in sorted(glob.glob(ROOT + "/runs/cavity/server_main/main_act0.07/rho*/k*/s[01][0-9]/directorfield.dat.gz"))
           if abs(json.load(open(os.path.dirname(p) + "/input.json"))["BC"][0]["R"] - 31.83) < 0.5]
    c1 = [(os.path.dirname(p), 5) for p in sorted(glob.glob(ROOT + "/runs/chi1_matched/R1.43_k*_s*/directorfield.dat.gz"))]
    with Pool(int(os.environ.get("NPROC", 12))) as p: r05 = p.map(run_one, c05, chunksize=1); r1 = p.map(run_one, c1, chunksize=1)
    Zp = np.load(os.path.join(H, "..", "curvature", "curvature_profiles.npz")); Zb = np.load(os.path.join(H, "chi1_band.npz"))
    prof = {"MPCD, chi = 0.5": (Zp["YC"], Zp["k"], Zp["alpha"], np.abs(Zp["R"] - 31.83) < 0.5), "MPCD, chi = 1": (Zb["YC"], Zb["k"], Zb["alpha"], np.ones(len(Zb["k"]), bool))}
    L = ["BAND WINDOWS (static-subtracted; coherent frames, exactly 2k-2 wall defects). Means over seeds of per-run medians.",
         "set | k | runs | tilt | d | a | edge d + a/2pi | l_b (band half-depth) | fixed [2,3): 2a_corr, gap | edge [edge, edge+1): 2a_meas, 2a_corr, gap"]
    out = {}
    for lab, rr in (("MPCD, chi = 0.5", r05), ("MPCD, chi = 1", r1)):
        ok = [r for r in rr if r and "err" not in r]; Y_, K_, A_, S_ = prof[lab]
        for e in [r for r in rr if r and "err" in r][:2]: L.append(f"  ERROR {e['run'][-40:]}: {e['err']}")
        for k in sorted({r["k"] for r in ok}):
            g = [r for r in ok if r["k"] == k]; f = lambda key: np.nanmean([r[key] for r in g])
            sel = S_ & (K_ == k); lb = band_thickness(Y_, np.nanmean(A_[sel], 0)) if sel.any() else np.nan
            row = dict(k=k, tau=f("tau"), d=f("d"), a=f("a"), edge=f("edge"), lb=lb, fc=f("fc"), em=f("em"), ec=f("ec"), n=len(g)); out.setdefault(lab, []).append(row)
            L.append(f"{lab:17s} | {k:2d} | {len(g):2d} | {row['tau']:5.1f} | {row['d']:4.2f} | {row['a']:5.1f} | {row['edge']:5.2f} | {lb:5.2f} | {row['fc']:5.1f}, {row['tau'] - row['fc']:+5.1f} | "
                     f"{row['em']:5.1f}, {row['ec']:5.1f}, {row['tau'] - row['ec']:+5.1f}")
    np.save(os.path.join(H, "band_windows.npy"), out, allow_pickle=True)
    open(os.path.join(H, "band_windows.txt"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
