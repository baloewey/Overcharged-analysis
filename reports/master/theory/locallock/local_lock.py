#!/usr/bin/env python3
"""Does the layer lock to the near-wall field, run by run? And what sets the depth of the wall defects? (user, 2026-10-22)
TEST 1 (lock to the local field). Per run: ring-averaged director angle from the local tangent, alpha(u) = arg<S e^{2i(theta - phi - pi/2)}>/2 signed by the
circulation (tiltcrowd/band_profile.py estimator, 1-cell bins of distance u from the wall, director frames of the second half); alpha_1 = alpha in the bin
just inside the layer (u in [d+0.5, d+1.5), d = median wall-defect depth), alpha_pk = max of alpha over u = 1-5 cells. Tilt tau of the coherent layer,
spacing a, depth d from tiltcrowd/defect_flow. Prediction 'lock to the local field': tau = 2 alpha_1 (or 2 alpha_pk); the excess tau - 2 alpha is what the
layer adds beyond the local field.
TEST 2 (depth). The one-constant equilibrium radius of a ring of N = 2k-2 wall +1/2 in a disk of radius R with a central +1 (master depth law, as used in
lockvalley/tiltcrowd: r_s = R ((N + 3)/(3N + 3))^(1/2N) - 0.61) gives a predicted depth R - r_s; compared with the measured median depth.
Runs: active main sweep with director output (server_main/main_act0.07, seeds 1-3), passive imposed vortex (swirl_crowd, swirl_shear), active switch-on reruns.
-> local_lock.txt, local_lock.npz"""
import os, sys, glob, gzip, json, warnings, numpy as np
from multiprocessing import Pool
warnings.filterwarnings("ignore")
ROOT = "/home/baloewe/Physics/Overcharged"; sys.path.insert(0, ROOT + "/reports/master/theory/tiltcrowd"); sys.path.insert(0, ROOT + "/runs/cavity")
import defect_flow as DF  # noqa
EDG = np.arange(0, 12.01, 1.0); UC = 0.5 * (EDG[1:] + EDG[:-1])
def frames(fn, ncol):
    op = gzip.open(fn, "rt") if fn.endswith(".gz") else open(fn)
    with op as fh: return np.loadtxt((l for l in fh if len(l.split()) == ncol and l.split()[0][0].isdigit()))
def pick(run, name):
    for f in (run + "/" + name + ".gz", run + "/" + name):
        if os.path.exists(f): return f
def one(run):
    try:
        bc = json.load(open(run + "/input.json"))["BC"][0]; R = bc["R"]; k = bc["customAnch"]; Q = np.array(bc["Q"][:2], float); N = int(2 * k - 2)
        if N < 4: return None
        a = DF.one(run)
        if not np.isfinite(a.get("tilt", np.nan)): return None
        D = frames(pick(run, "directorfield.dat"), 8); Fl = frames(pick(run, "flowfield.dat"), 7)
        td = np.unique(D[:, 0]); td = td[td >= td[0] + 0.5 * (td[-1] - td[0])]; Z = np.zeros(len(EDG) - 1, complex); n = 0
        for tt in td:
            d = D[D[:, 0] == tt]; f = Fl[np.abs(Fl[:, 0] - tt) < 1e-6]
            if not len(f): continue
            xf = f[:, 1] + 0.5 - Q[0]; yf = f[:, 2] + 0.5 - Q[1]; rf = np.hypot(xf, yf)
            circ = np.sign(((xf * f[:, 5] - yf * f[:, 4]) / np.maximum(rf, 1e-9))[rf < R - 1].mean())
            if circ == 0: continue
            x = d[:, 1] + 0.5 - Q[0]; y = d[:, 2] + 0.5 - Q[1]; r = np.hypot(x, y); ph = np.arctan2(y, x); u = R - r; m = (u >= 0) & (u < EDG[-1]) & (d[:, 7] > 0)
            z = d[m, 7] * np.exp(2j * circ * (np.arctan2(d[m, 5], d[m, 4]) - ph[m] - np.pi / 2)); b = np.digitize(u[m], EDG) - 1
            for i in range(len(EDG) - 1):
                q = b == i
                if q.any(): Z[i] += z[q].mean()
            n += 1
        if n < 2: return None
        al = np.degrees(np.angle(Z / n) / 2); dep = a["depth"]
        i1 = np.argmin(np.abs(UC - (dep + 1.0))); m5 = (UC >= 1) & (UC <= 5)
        rs = R * ((N + 3) / (3 * N + 3)) ** (1 / (2 * N)) - 0.61
        return dict(run=os.path.relpath(run, ROOT + "/runs"), k=k, R=R, a=a["a"], d=dep, d_law=R - rs, tau=abs(a["tilt"]), a1=al[i1], apk=np.max(al[m5]), prof=al)
    except Exception as ex:
        return dict(run=run, err=repr(ex))
if __name__ == "__main__":
    SETS = {"active": sorted(os.path.dirname(p) for p in glob.glob(ROOT + "/runs/cavity/server_main/main_act0.07/rho*/k*/s0[1-3]/directorfield.dat.gz")) + [ROOT + "/runs/switchon2/dir_k06_alignX/active", ROOT + "/runs/switchon2/dir_k06_rand_s3/active"],
            "passive": sorted(r for r in glob.glob(ROOT + "/runs/swirl_crowd/*k*_s*") + glob.glob(ROOT + "/runs/swirl_shear/k*_A*") if os.path.exists(r + "/DONE"))}
    L, OUT = [], {}
    for nm, runs in SETS.items():
        with Pool(8) as p: res = p.map(one, runs, chunksize=1)
        ok = [x for x in res if x and "err" not in x]; err = [x for x in res if x and "err" in x]
        A = {q: np.array([x[q] for x in ok]) for q in ("k", "a", "d", "d_law", "tau", "a1", "apk")}; OUT[nm] = A
        L.append(f"== {nm}: {len(ok)} runs ({len(err)} errors, {len(runs) - len(ok) - len(err)} without a coherent layer or director frames)")
        for q in ("a1", "apk"):
            x2 = 2 * A[q]; sl, ic = np.polyfit(x2, A["tau"], 1); res_ = A["tau"] - x2
            L.append(f"  tau vs 2*{q}: corr {np.corrcoef(x2, A['tau'])[0,1]:+.2f}, fit tau = {ic:+.1f} + {sl:.2f} (2 {q}); excess tau - 2 {q}: median {np.median(res_):+.1f}, corr with log spacing {np.corrcoef(res_, np.log(A['a']))[0,1]:+.2f}")
        L.append("  by spacing: n | tilt | 2 alpha_1 | excess (tau - 2 alpha_1) | 2 alpha_pk | depth measured | depth law")
        for lo, hi in ((0, 8), (8, 12), (12, 20), (20, 80)):
            m = (A["a"] >= lo) & (A["a"] < hi)
            if m.sum(): L.append(f"  [{lo:2d},{hi:2d}) {m.sum():3d} | {np.median(A['tau'][m]):5.1f} | {np.median(2*A['a1'][m]):5.1f} | {np.median(A['tau'][m]-2*A['a1'][m]):+5.1f} | {np.median(2*A['apk'][m]):5.1f} | {np.median(A['d'][m]):.2f} | {np.median(A['d_law'][m]):.2f}")
        sl, ic = np.polyfit(A["d_law"], A["d"], 1)
        L.append(f"  DEPTH: measured vs law: corr {np.corrcoef(A['d_law'], A['d'])[0,1]:+.2f}, fit d = {ic:+.2f} + {sl:.2f} d_law, rms {np.std(A['d'] - (ic + sl*A['d_law'])):.2f}; median measured - law {np.median(A['d'] - A['d_law']):+.2f}; corr(depth, log spacing) {np.corrcoef(A['d'], np.log(A['a']))[0,1]:+.2f}")
        for x in err[:3]: L.append(f"  ERROR {x['run']}: {x['err']}")
    np.savez("local_lock.npz", **{f"{nm}|{q}": v for nm, A in OUT.items() for q, v in A.items()})
    open("local_lock.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
