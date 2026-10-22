#!/usr/bin/env python3
"""1D boundary-layer model of the near-wall director with the MEASURED shear profile (user, 2026-10-22; candidate B).
Director angle alpha(u) from the local tangent (u = distance from the wall, signed so that alpha > 0 is the flow-aligned sense). For a steady azimuthal
flow v(r) the lab-frame director advected by the flow rotates relative to the local tangent only through the pure shear e(r) = (r/2) d(v/r)/dr (rigid
rotation gives no torque), so in the one-constant Ericksen-Leslie limit (curvature terms dropped)
    0 = K alpha'' + chi gamma1 e(u) (lambda cos 2 alpha - 1),     K = 44, gamma1 = 60, chi = 0.5, lambda = 2 (alpha_L = 30 deg where e > 0).
The wall-defect layer at depth d is a boundary that carries the lock: alpha(d) = tau/2; interior boundary alpha'(U) = 0 at U = 0.85 R.
 M0 (no torque at the layer, alpha'(d) = 0): predicted tilt tau0 = 2 alpha(d).
 M1 (measured tilt imposed, alpha(d) = tau/2): predicted field one cell inside the layer, 2 alpha(d+1), hence the lead tau - 2 alpha(d+1), and the torque
    per unit wall length the layer must receive to hold that tilt, G = K alpha'(d) (kBT/cell; > 0 = the layer is pulled back by the outer field).
Compared with the measured coherent-frame tilt, near-wall field 2 alpha_1 and lead (locallock/static_split.one), and with the band-ansatz torque
chi gamma1 e_w d / 2 (e_w = shear at the layer). Shear profile: azimuthally averaged v_phi of the flow frames of the second half, signed by the
circulation, 1-cell bins, e(u) from a smoothed profile. Runs: as static_split.py. -> bl1d.txt, bl1d.npz"""
import os, sys, glob, gzip, json, warnings, numpy as np
from multiprocessing import Pool
warnings.filterwarnings("ignore")
ROOT = "/home/baloewe/Physics/Overcharged"; sys.path.insert(0, ROOT + "/reports/master/theory/locallock")
import static_split as SS  # noqa
K, G1, CHI, LAM = 44.0, 60.0, 0.5, 2.0; DU = 0.25
def flow_profile(run, R, Q):
    fn = SS.pick(run, "flowfield.dat")
    with (gzip.open(fn, "rt") if fn.endswith(".gz") else open(fn)) as fh: F = np.loadtxt((l for l in fh if len(l.split()) == 7 and l.split()[0][0].isdigit()))
    tf = np.unique(F[:, 0]); tf = tf[tf >= tf[0] + 0.5 * (tf[-1] - tf[0])]; edges = np.arange(0, R + 0.01, 1.0); rc = 0.5 * (edges[1:] + edges[:-1]); P = []
    for tt in tf:
        f = F[F[:, 0] == tt]; x = f[:, 1] + 0.5 - Q[0]; y = f[:, 2] + 0.5 - Q[1]; r = np.hypot(x, y); m = r < R
        vt = (x * f[:, 5] - y * f[:, 4]) / np.maximum(r, 1e-9); h, _ = np.histogram(r[m], edges, weights=vt[m]); c, _ = np.histogram(r[m], edges)
        p = h / np.maximum(c, 1); sg = np.sign(np.sum(p * rc)) or 1.0; P.append(sg * p)
    vm = np.mean(P, 0); vs = np.convolve(vm, np.ones(3) / 3, mode="same"); w = vs / np.maximum(rc, 1e-9)
    e = 0.5 * rc * np.gradient(w, rc)                                  # pure shear e(r) = (r/2) d(v/r)/dr
    return R - rc, -e                                                   # in u = R - r; sign so that wall-sign shear (v grows away from the wall) is > 0
def solve(ug, eg, mode, a0, dt=None):
    a = np.full(len(ug), np.radians(25.0)); a[0] = a0 if mode == "M1" else a[0]; D = K / G1 / CHI if False else K / G1
    dt = 0.2 * DU ** 2 / D
    for _ in range(200000):
        lap = np.zeros_like(a); lap[1:-1] = (a[2:] - 2 * a[1:-1] + a[:-2]) / DU ** 2
        src = CHI * eg * (LAM * np.cos(2 * a) - 1)
        if mode == "M0": lap[0] = 2 * (a[1] - a[0]) / DU ** 2           # alpha'(d) = 0
        lap[-1] = 2 * (a[-2] - a[-1]) / DU ** 2                          # alpha'(U) = 0
        da = dt * (D * lap + src)
        if mode == "M1": da[0] = 0.0
        a += da
        if np.max(np.abs(da)) < 1e-9: break
    return a
def one(run):
    try:
        r = SS.one(run)
        if not r or "err" in r: return None
        bc = json.load(open(run + "/input.json"))["BC"][0]; R = bc["R"]; Q = np.array(bc["Q"][:2], float)
        u, e = flow_profile(run, R, Q); d = r["d"]; U = 0.85 * R
        ug = np.arange(d, U + 1e-9, DU); eg = np.interp(ug, u[::-1], e[::-1])
        a0 = solve(ug, eg, "M0", None); a1 = solve(ug, eg, "M1", np.radians(r["tau"] / 2))
        i1 = int(round(1.0 / DU)); G = K * (a1[1] - a1[0]) / DU * -1      # K alpha'(d): > 0 when the field falls inward (pulls the layer back)
        ew = float(np.interp(d, u[::-1], e[::-1]))
        return dict(run=r["run"], k=r["k"], a=r["a"], d=d, tau=r["tau"], f_m=r["f_m"], tau0=np.degrees(2 * a0[0]), f1=np.degrees(2 * a1[i1]),
                    G=G, ew=ew, Gans=CHI * G1 * ew * d / 2, ug=ug, a1=np.degrees(a1), eg=eg)
    except Exception as ex:
        return dict(run=run, err=repr(ex))
if __name__ == "__main__":
    SETS = {"active": sorted(os.path.dirname(p) for p in glob.glob(ROOT + "/runs/cavity/server_main/main_act0.07/rho*/k*/s0[1-3]/directorfield.dat.gz")) + [ROOT + "/runs/switchon2/dir_k06_alignX/active", ROOT + "/runs/switchon2/dir_k06_rand_s3/active"],
            "passive": sorted(r for r in glob.glob(ROOT + "/runs/swirl_crowd/*k*_s*") + glob.glob(ROOT + "/runs/swirl_shear/k*_A*") if os.path.exists(r + "/DONE"))}
    L, OUT = ["1D BOUNDARY-LAYER MODEL with the measured shear profile (deg; medians per spacing bin). M0: no torque at the layer -> predicted tilt. M1: measured tilt imposed -> predicted near-wall field 2 alpha(d+1), predicted lead, required torque G = K alpha'(d) (kBT/cell) vs band ansatz chi gamma1 e_w d/2."], {}
    for nm, runs in SETS.items():
        with Pool(8) as p: res = p.map(one, runs, chunksize=1)
        ok = [x for x in res if x and "err" not in x]; err = [x for x in res if x and "err" in x]
        A = {q: np.array([x[q] for x in ok]) for q in ("k", "a", "d", "tau", "f_m", "tau0", "f1", "G", "ew", "Gans")}; OUT[nm] = A
        L.append(f"== {nm}: {len(ok)} runs ({len(err)} errors)")
        L.append("  spacing | n | tilt meas | M0 tilt | 2a1 meas | M1 2a(d+1) | lead meas | M1 lead | wall shear e_w | required torque G | band ansatz")
        for lo, hi in ((0, 8), (8, 12), (12, 20), (20, 80)):
            m = (A["a"] >= lo) & (A["a"] < hi)
            if m.sum(): L.append(f"  [{lo:2d},{hi:2d}) | {m.sum():3d} | {np.median(A['tau'][m]):5.1f} | {np.median(A['tau0'][m]):5.1f} | {np.median(A['f_m'][m]):5.1f} | {np.median(A['f1'][m]):5.1f} | {np.median(A['tau'][m]-A['f_m'][m]):+5.1f} | {np.median(A['tau'][m]-A['f1'][m]):+5.1f} | {np.median(A['ew'][m]):.3f} | {np.median(A['G'][m]):5.2f} | {np.median(A['Gans'][m]):5.2f}")
        L.append(f"  corr(tilt, M0 tilt) {np.corrcoef(A['tau'], A['tau0'])[0,1]:+.2f}; corr(lead meas, M1 lead) {np.corrcoef(A['tau']-A['f_m'], A['tau']-A['f1'])[0,1]:+.2f}; corr(required torque, band ansatz) {np.corrcoef(A['G'], A['Gans'])[0,1]:+.2f}; corr(required torque, depth) {np.corrcoef(A['G'], A['d'])[0,1]:+.2f}")
        for x in err[:3]: L.append(f"  ERROR {x['run']}: {x['err']}")
    np.savez("bl1d.npz", **{f"{nm}|{q}": v for nm, A in OUT.items() for q, v in A.items()})
    open("bl1d.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
