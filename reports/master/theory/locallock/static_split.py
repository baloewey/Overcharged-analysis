#!/usr/bin/env python3
"""Test A (user, 2026-10-22): how much of the split tilt = 2 alpha_1 + excess is static geometry? For each director frame (second half) with exactly
2k-2 wall +1/2 and total charge k, the exact static one-constant field of the SAME tracked defects (ejection/tilted_structures.fields; anchoring theta
= k phi, images) is analysed with the same estimator as the measured director: alpha_1 = ring average (all cells at u in [d+0.5, d+1.5), d = median
wall-defect depth of the frame) of the angle from the local tangent; tilt of the static field = mean natural tilt of its wall defects
(tilted_structures.natural_tilts). Everything signed by the measured tilt. Per run: measured tilt, 2 alpha_1, excess; static tilt, 2 alpha_1, excess.
If the static split reproduces the crowding trend of the measured one, the split is geometry (lock + images + depth) and the physics is in tau itself.
Coherent frames only (C > 0.9). Also 2 alpha at u = d + a/2pi (the range of the layer's periodic distortion), to test whether a layer locks to the field at the edge of its own pattern. Runs as local_lock.py. -> static_split.txt, static_split.npz"""
import os, sys, glob, gzip, json, warnings, numpy as np
from multiprocessing import Pool
warnings.filterwarnings("ignore")
ROOT = "/home/baloewe/Physics/Overcharged"; T_ = ROOT + "/reports/master/theory"; sys.path.insert(0, ROOT + "/runs/cavity")
from edge_criteria import read_with_angle  # noqa
import importlib.util
spec = importlib.util.spec_from_file_location("ts", T_ + "/ejection/tilted_structures.py"); TS = importlib.util.module_from_spec(spec)
exec(compile(open(T_ + "/ejection/tilted_structures.py").read().split('\nif __name__ == "__main__":')[0], "ts", "exec"), TS.__dict__)
def pick(run, name):
    for f in (run + "/" + name + ".gz", run + "/" + name):
        if os.path.exists(f): return f
def one(run):
    try:
        bc = json.load(open(run + "/input.json"))["BC"][0]; R = bc["R"]; k = bc["customAnch"]; Q = np.array(bc["Q"][:2], float); N = int(2 * k - 2)
        if N < 4: return None
        fn = pick(run, "directorfield.dat")
        with (gzip.open(fn, "rt") if fn.endswith(".gz") else open(fn)) as fh: D = np.loadtxt((l for l in fh if len(l.split()) == 8 and l.split()[0][0].isdigit()))
        td = np.unique(D[:, 0]); td = td[td >= td[0] + 0.5 * (td[-1] - td[0])]; t, fr = read_with_angle(run, Q, 0.0); t = np.array(t)
        n = int(D[:, 1].max()) + 1; m = int(D[:, 2].max()) + 1; X, Y = np.meshgrid(np.arange(n) + 0.5 - Q[0], np.arange(m) + 0.5 - Q[1], indexing="ij")
        r = np.hypot(X, Y); ph = np.arctan2(Y, X); u = R - r; ins = r < R
        rows = []
        for tt in td:
            i = int(np.argmin(np.abs(t - tt)))
            if abs(t[i] - tt) > 1e-6: continue
            f = fr[i]
            if abs(f[:, 2].sum() - k) > 1e-6: continue
            rr = np.hypot(f[:, 0], f[:, 1]); w = (f[:, 2] > 0) & (R - rr < 4)
            if w.sum() != N: continue
            ez = np.mean(np.exp(1j * (f[w, 3] - np.arctan2(f[w, 1], f[w, 0]))))
            if abs(ez) < 0.9: continue
            tau = np.angle(ez); sg = np.sign(tau) or 1.0; d = np.median(R - rr[w]); band = ins & (u >= d + 0.5) & (u < d + 1.5)
            dd = D[D[:, 0] == tt]; th2 = np.full((n, m), np.nan); S = np.zeros((n, m)); ix, iy = dd[:, 1].astype(int), dd[:, 2].astype(int)
            th2[ix, iy] = 2 * np.arctan2(dd[:, 5], dd[:, 4]); S[ix, iy] = dd[:, 7]
            zs = list(f[:, 0] + 1j * f[:, 1]); qs = list(f[:, 2]); th2s, _, _ = TS.fields(X, Y, zs, qs, R)
            sel = band & np.isfinite(th2)
            am = np.angle(np.sum(S[sel] * np.exp(1j * sg * (th2[sel] - 2 * ph[sel] - np.pi)))) / 2
            ast = np.angle(np.sum(np.exp(1j * sg * (th2s[band] - 2 * ph[band] - np.pi)))) / 2
            nt = TS.natural_tilts(zs, qs, R, len(zs)); wi = np.flatnonzero(w); ts_ = sg * np.angle(np.mean(np.exp(1j * nt[wi])))
            sp = 2 * np.pi * (R - d) / N; uq = d + sp / (2 * np.pi); bq = ins & (u >= uq - 0.5) & (u < uq + 0.5) & np.isfinite(th2)
            aq = np.angle(np.sum(S[bq] * np.exp(1j * sg * (th2[bq] - 2 * ph[bq] - np.pi)))) / 2
            rows.append((abs(tau), 2 * am, ts_, 2 * ast, 2 * aq, d))
        if len(rows) < 3: return None
        a = np.degrees(np.array(rows)[:, :5]); dep = np.median(np.array(rows)[:, 5])
        return dict(run=os.path.relpath(run, ROOT + "/runs"), k=k, a=2 * np.pi * (R - dep) / N, d=dep, n=len(rows),
                    tau=np.median(a[:, 0]), f_m=np.median(a[:, 1]), tau_s=np.median(a[:, 2]), f_s=np.median(a[:, 3]), f_q=np.median(a[:, 4]))
    except Exception as ex:
        return dict(run=run, err=repr(ex))
if __name__ == "__main__":
    SETS = {"active": sorted(os.path.dirname(p) for p in glob.glob(ROOT + "/runs/cavity/server_main/main_act0.07/rho*/k*/s0[1-3]/directorfield.dat.gz")) + [ROOT + "/runs/switchon2/dir_k06_alignX/active", ROOT + "/runs/switchon2/dir_k06_rand_s3/active"],
            "passive": sorted(r for r in glob.glob(ROOT + "/runs/swirl_crowd/*k*_s*") + glob.glob(ROOT + "/runs/swirl_shear/k*_A*") if os.path.exists(r + "/DONE"))}
    L, OUT = ["Measured split: tilt = 2 alpha_1 + excess; static split: static tilt (natural tilts of the same defects) = 2 alpha_1,static + excess_static (deg; medians per spacing bin)"], {}
    for nm, runs in SETS.items():
        with Pool(8) as p: res = p.map(one, runs, chunksize=1)
        ok = [x for x in res if x and "err" not in x]; err = [x for x in res if x and "err" in x]
        A = {q: np.array([x[q] for x in ok]) for q in ("k", "a", "d", "tau", "f_m", "tau_s", "f_s", "f_q")}; OUT[nm] = A
        L.append(f"== {nm}: {len(ok)} runs ({len(err)} errors)")
        L.append("  spacing | n | tilt | 2a1 meas | excess meas | static tilt | 2a1 static | excess static | flow part of the field (meas - static) | 2 alpha at u = d + a/2pi (pattern range) | tilt - that")
        for lo, hi in ((0, 8), (8, 12), (12, 20), (20, 80)):
            m = (A["a"] >= lo) & (A["a"] < hi)
            if m.sum(): L.append(f"  [{lo:2d},{hi:2d}) | {m.sum():3d} | {np.median(A['tau'][m]):5.1f} | {np.median(A['f_m'][m]):5.1f} | {np.median(A['tau'][m]-A['f_m'][m]):+5.1f} | {np.median(A['tau_s'][m]):5.1f} | {np.median(A['f_s'][m]):5.1f} | {np.median(A['tau_s'][m]-A['f_s'][m]):+5.1f} | {np.median(A['f_m'][m]-A['f_s'][m]):+5.1f} | {np.median(A['f_q'][m]):5.1f} | {np.median(A['tau'][m]-A['f_q'][m]):+5.1f}")
        L.append(f"  corr(tilt, 2 alpha_1) {np.corrcoef(A['tau'], A['f_m'])[0,1]:+.2f}; corr(tilt, 2 alpha(d + a/2pi)) {np.corrcoef(A['tau'], A['f_q'])[0,1]:+.2f}; rms(tilt - 2 alpha(d + a/2pi)) {np.std(A['tau'] - A['f_q']):.1f}")
        for x in err[:3]: L.append(f"  ERROR {x['run']}: {x['err']}")
    np.savez("static_split.npz", **{f"{nm}|{q}": v for nm, A in OUT.items() for q, v in A.items()})
    open("static_split.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
