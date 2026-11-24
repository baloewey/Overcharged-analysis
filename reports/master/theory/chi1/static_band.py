#!/usr/bin/env python3
"""Static-subtracted band angle (user, 2026-11-24). The ring average of the director just inside the layer contains the wall defects' own near field
(decays over a/2pi; biased in sparse walls). Per coherent frame (C > 0.9, exactly 2k-2 wall +1/2, total charge k), the exact one-constant static texture
of the SAME tracked defects (ejection/tilted_structures.fields: anchoring theta = k phi, images) is ring-averaged on the same cells as the measured field
(depth window [2, 3) cells, as fig_band_vs_tilt). The static average differs from the zero-mode value tau_s/2 (tau_s = natural tilt of the same defects)
only by the near field, so the corrected band is
    alpha_corr = alpha_meas - (alpha_static - tau_s / 2)     ->  2 alpha_corr = 2 alpha_meas - 2 alpha_static + tau_s
Gap tilt - 2 alpha_corr = (tau - tau_s) - 2 (alpha_meas - alpha_static): the flow-induced rotation of the field relative to the static texture.
Also reported: the plain ring average 2 alpha_meas (fig_band_vs_tilt) and the near-field bias 2 alpha_static - tau_s.
Sets: MPCD chi 0.5 (main sweep, R 31.83, local director fields) and chi 1 (chi1_matched R1.43, R 21.6). Continuum: continuum_static_band() reads
runs/continuum_chi/snap_*_snaps.npz (Q snapshots + defect rows). -> static_band.npz (per run), static_band.txt"""
import glob, gzip, json, os, sys, warnings, numpy as np
from multiprocessing import Pool
warnings.filterwarnings("ignore")
ROOT = "/home/baloewe/Physics/Overcharged"; T_ = ROOT + "/reports/master/theory"; sys.path.insert(0, ROOT + "/runs/cavity")
from edge_criteria import read_with_angle  # noqa
import importlib.util
spec = importlib.util.spec_from_file_location("ts", T_ + "/ejection/tilted_structures.py"); TS = importlib.util.module_from_spec(spec)
exec(compile(open(T_ + "/ejection/tilted_structures.py").read().split('\nif __name__ == "__main__":')[0], "ts", "exec"), TS.__dict__)
H = os.path.dirname(os.path.abspath(__file__)); Y0, Y1 = 2.0, 3.0


def frame_values(X, Y, R, k, th2, S, zs, qs, tails):
    """th2: measured 2 theta on the grid (nan outside); zs, qs: all defects; tails: tail angle of each defect (nan if unknown). Returns tau, 2a_meas, 2a_stat, tau_s or None."""
    zs = np.asarray(zs); qs = np.asarray(qs); r = np.abs(zs); N = 2 * k - 2; w = (qs > 0) & (R - r < 4)
    if w.sum() != N or abs(qs.sum() - k) > 1e-6: return None
    ez = np.mean(np.exp(1j * (tails[w] - np.angle(zs[w]))))
    if abs(ez) < 0.9: return None
    tau = np.angle(ez); sg = np.sign(tau) or 1.0; rr = np.hypot(X, Y); ph = np.arctan2(Y, X); u = R - rr; band = (u >= Y0) & (u < Y1) & np.isfinite(th2)
    th2s, _, _ = TS.fields(X, Y, list(zs), list(qs), R)
    am = np.angle(np.sum(S[band] * np.exp(1j * sg * (th2[band] - 2 * ph[band] - np.pi)))) / 2
    ast = np.angle(np.sum(np.exp(1j * sg * (th2s[band] - 2 * ph[band] - np.pi)))) / 2
    nt = TS.natural_tilts(list(zs), list(qs), R, len(zs)); ts_ = sg * np.angle(np.mean(np.exp(1j * nt[np.flatnonzero(w)])))
    return abs(tau), 2 * am, 2 * ast, ts_


def mpcd_run(args):
    run, stride = args
    try:
        bc = json.load(open(run + "/input.json"))["BC"][0]; R = bc["R"]; k = bc["customAnch"]; Q = np.array(bc["Q"][:2], float)
        with gzip.open(run + "/directorfield.dat.gz", "rt") as fh: D = np.loadtxt((l for l in fh if len(l.split()) == 8 and l.split()[0][0].isdigit()))
        td = np.unique(D[:, 0]); td = td[td >= td[0] + 0.5 * (td[-1] - td[0])][::stride]; t, fr = read_with_angle(run, Q, 0.0); t = np.array(t)
        n = int(D[:, 1].max()) + 1; m = int(D[:, 2].max()) + 1; X, Y = np.meshgrid(np.arange(n) + 0.5 - Q[0], np.arange(m) + 0.5 - Q[1], indexing="ij")
        rows = []
        for tt in td:
            i = int(np.argmin(np.abs(t - tt)))
            if abs(t[i] - tt) > 1e-6: continue
            f = fr[i]; dd = D[D[:, 0] == tt]; th2 = np.full((n, m), np.nan); S = np.zeros((n, m)); ix, iy = dd[:, 1].astype(int), dd[:, 2].astype(int)
            th2[ix, iy] = 2 * np.arctan2(dd[:, 5], dd[:, 4]); S[ix, iy] = dd[:, 7]; th2[np.hypot(X, Y) >= R] = np.nan
            v = frame_values(X, Y, R, k, th2, S, f[:, 0] + 1j * f[:, 1], f[:, 2], f[:, 3])
            if v: rows.append(v)
        if len(rows) < 3: return None
        a = np.degrees(np.array(rows)); return dict(run=run, k=k, R=R, n=len(rows), tau=np.median(a[:, 0]), am=np.median(a[:, 1]), ast=np.median(a[:, 2]), tss=np.median(a[:, 3]))
    except Exception as ex:
        return dict(run=run, err=repr(ex))


def cont_run(fn):
    """continuum snapshot file: director from Q, defects from the saved rows, tails from div Q at each core (as ellipse_active_K.observe)."""
    Z = np.load(fn); x, y = Z["x"], Z["y"]; X, Y = np.meshgrid(x, y, indexing="ij"); R = float(Z["R"]); h = x[1] - x[0]
    J = json.load(open(fn.replace("_snaps.npz", ".json"))); k = J["args"]["k"]; off = np.r_[0, np.cumsum(Z["n"])]; rows = []
    dX = lambda f: (np.roll(f, -1, 0) - np.roll(f, 1, 0)) / (2 * h); dY = lambda f: (np.roll(f, -1, 1) - np.roll(f, 1, 1)) / (2 * h)
    for i, tt in enumerate(Z["t"]):
        if tt < 500: continue
        q1, q2 = Z["q1"][i].astype(float), Z["q2"][i].astype(float); S = 2 * np.hypot(q1, q2); th2 = np.arctan2(q2, q1); th2[np.hypot(X, Y) >= R] = np.nan
        d = Z["rows"][off[i]:off[i + 1]]; vx, vy = dX(q1) + dY(q2), dX(q2) - dY(q1); tails = []
        for px, py, q in d:
            mk = np.hypot(X - px, Y - py) < 1.6; tails.append(np.arctan2(vy[mk].mean(), vx[mk].mean()) if q > 0 else np.nan)
        v = frame_values(X, Y, R, k, th2, S, d[:, 0] + 1j * d[:, 1], d[:, 2], np.array(tails))
        if v: rows.append(v)
    if len(rows) < 2: return None
    a = np.degrees(np.array(rows)); return dict(run=fn, k=k, R=R, n=len(rows), tau=np.median(a[:, 0]), am=np.median(a[:, 1]), ast=np.median(a[:, 2]), tss=np.median(a[:, 3]))


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "mpcd"
    if which == "mpcd":
        c05 = [(os.path.dirname(p), 1) for p in sorted(glob.glob(ROOT + "/runs/cavity/server_main/main_act0.07/rho*/k*/s[01][0-9]/directorfield.dat.gz"))
               if abs(json.load(open(os.path.dirname(p) + "/input.json"))["BC"][0]["R"] - 31.83) < 0.5]
        c1 = [(os.path.dirname(p), 5) for p in sorted(glob.glob(ROOT + "/runs/chi1_matched/R1.43_k*_s*/directorfield.dat.gz"))]
        with Pool(int(os.environ.get("NPROC", 12))) as p: r05 = p.map(mpcd_run, c05, chunksize=1); r1 = p.map(mpcd_run, c1, chunksize=1)
        res = {"MPCD, chi = 0.5": r05, "MPCD, chi = 1": r1}
    else:
        fs = sorted(glob.glob(ROOT + "/runs/continuum_chi/snap_*_snaps.npz"))
        with Pool(8) as p: rc = p.map(cont_run, fs)
        res = {}
        for f, r in zip(fs, rc):
            if not r: continue
            b = os.path.basename(f); lab = f"Continuum {'laminar' if '_z1_' in b else 'turbulent'}, chi = {'0.5' if 'chi0.5' in b else '1'}"; res.setdefault(lab, []).append(r)
    L = [f"STATIC-SUBTRACTED BAND ({which}); depth window [{Y0:g}, {Y1:g}) cells; medians over coherent frames per run, then mean over seeds",
         "set | k | runs | tilt | 2 a_meas (plain ring avg) | near-field bias 2 a_stat - tau_s | 2 a_corr (static-subtracted) | gap tilt - 2 a_corr"]
    out = {}
    for lab, rr in res.items():
        ok = [r for r in rr if r and "err" not in r]; err = [r for r in rr if r and "err" in r]
        for k in sorted({r["k"] for r in ok}):
            g = [r for r in ok if r["k"] == k]; f = lambda key: np.mean([r[key] for r in g])
            corr = f("am") - f("ast") + f("tss"); out.setdefault(lab, []).append((k, f("tau"), f("am"), corr, len(g)))
            L.append(f"{lab:30s} | {k:2d} | {len(g):2d} | {f('tau'):5.1f} | {f('am'):5.1f} | {f('ast') - f('tss'):+5.1f} | {corr:5.1f} | {f('tau') - corr:+5.1f}")
        for e in err[:2]: L.append(f"  ERROR {e['run'][-40:]}: {e['err'][:100]}")
    np.save(os.path.join(H, f"static_band_{which}.npy"), out, allow_pickle=True)
    open(os.path.join(H, f"static_band_{which}.txt"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
