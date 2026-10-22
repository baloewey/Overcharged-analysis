#!/usr/bin/env python3
"""Does the defect core rotate further than its surroundings? (user, 2026-10-22). Both models, same estimator.
For a wall +1/2 at z0: core orientation psi_c = direction of div Q averaged within 1.6 cells (the tracker/tilt convention); ring orientation
psi_r(r1, r2) = arg < exp(i (2 theta - phi')) > over cells at distance r1-r2 from z0 (phi' = polar angle about z0; for a +1/2 with director
theta = phi'/2 + theta0 this is the radial-arm direction 2 theta0), cells inside the cavity and > 1.5 from other defects. The same two numbers on the
exact static one-constant field of the same defect positions remove the wall/image distortion:
core rotation = (psi_c - psi_r)_measured - (psi_c - psi_r)_static, signed by the tilt (+ = core rotated further in the flow sense).
MPCD: passive imposed-vortex runs (swirl_crowd, swirl_shear) and active switch-on reruns (switchon2 dir_k06_*), director frames of the second half.
Continuum: flat wall, flatwall_band.run machinery re-done here (s 0.03/0.06, a_s 8/12/20), static field flat_pinning.setup.
-> core_rotation.txt"""
import os, sys, glob, gzip, json, warnings, numpy as np
from multiprocessing import Pool
warnings.filterwarnings("ignore")
ROOT = "/home/baloewe/Physics/Overcharged"; T_ = ROOT + "/reports/master/theory"; sys.path.insert(0, ROOT + "/runs/cavity"); sys.path.insert(0, T_ + "/flatpin"); sys.path.insert(0, T_ + "/shear_tilt")
from edge_criteria import read_with_angle  # noqa
import importlib.util
spec = importlib.util.spec_from_file_location("ts", T_ + "/ejection/tilted_structures.py"); TS = importlib.util.module_from_spec(spec)
exec(compile(open(T_ + "/ejection/tilted_structures.py").read().split('\nif __name__ == "__main__":')[0], "ts", "exec"), TS.__dict__)
RINGS = ((2.0, 3.0), (3.0, 4.5))
def divq_dir(q1, q2, X, Y, x0, y0, h):
    dvx = np.gradient(q1, h, axis=0) + np.gradient(q2, h, axis=1); dvy = np.gradient(q2, h, axis=0) - np.gradient(q1, h, axis=1)
    m = np.hypot(X - x0, Y - y0) < 1.6; return np.arctan2(dvy[m].mean(), dvx[m].mean())
def ring_dir(th2, X, Y, x0, y0, valid, r1, r2):
    rr = np.hypot(X - x0, Y - y0); m = valid & (rr >= r1) & (rr < r2); ph = np.arctan2(Y - y0, X - x0)
    return np.angle(np.mean(np.exp(1j * (th2[m] - ph[m])))) if m.sum() > 4 else np.nan
def per_defect(th2m, Sm, th2s, X, Y, x0, y0, valid, h):
    q1m, q2m = Sm * np.cos(th2m), Sm * np.sin(th2m); q1s, q2s = np.cos(th2s), np.sin(th2s)
    cm = divq_dir(q1m, q2m, X, Y, x0, y0, h); cs = divq_dir(q1s, q2s, X, Y, x0, y0, h); out = []
    for r1, r2 in RINGS:
        rm = ring_dir(th2m, X, Y, x0, y0, valid, r1, r2); rs = ring_dir(th2s, X, Y, x0, y0, valid, r1, r2)
        out.append(np.angle(np.exp(1j * ((cm - rm) - (cs - rs)))))
    return out
def mpcd(run):
    try:
        bc = json.load(open(run + "/input.json"))["BC"][0]; R = bc["R"]; k = bc["customAnch"]; Q = np.array(bc["Q"][:2], float); N = int(2 * k - 2)
        fn = run + "/directorfield.dat.gz" if os.path.exists(run + "/directorfield.dat.gz") else run + "/directorfield.dat"
        with (gzip.open(fn, "rt") if fn.endswith(".gz") else open(fn)) as fh: D = np.loadtxt((l for l in fh if len(l.split()) == 8 and l.split()[0][0].isdigit()))
        td = np.unique(D[:, 0]); td = td[td >= td[0] + 0.5 * (td[-1] - td[0])]; t, fr = read_with_angle(run, Q, 0.0); t = np.array(t)
        n = int(D[:, 1].max()) + 1; m = int(D[:, 2].max()) + 1; X, Y = np.meshgrid(np.arange(n) + 0.5 - Q[0], np.arange(m) + 0.5 - Q[1], indexing="ij"); ins = np.hypot(X, Y) < R
        vals = [[], []]; tilts = []
        for tt in td:
            i = int(np.argmin(np.abs(t - tt)))
            if abs(t[i] - tt) > 1e-6: continue
            f = fr[i]
            if abs(f[:, 2].sum() - k) > 1e-6: continue
            rr = np.hypot(f[:, 0], f[:, 1]); w = (f[:, 2] > 0) & (R - rr < 4)
            if w.sum() != N: continue
            tau = np.angle(np.mean(np.exp(1j * (f[w, 3] - np.arctan2(f[w, 1], f[w, 0]))))); sg = np.sign(tau) or 1.0; tilts.append(abs(tau))
            d = D[D[:, 0] == tt]; th2 = np.zeros((n, m)); S = np.zeros((n, m)); ix, iy = d[:, 1].astype(int), d[:, 2].astype(int)
            th2[ix, iy] = 2 * np.arctan2(d[:, 5], d[:, 4]); S[ix, iy] = d[:, 7]
            th2s, _, _ = TS.fields(X, Y, list(f[:, 0] + 1j * f[:, 1]), list(f[:, 2]), R)
            for x0, y0 in zip(f[w, 0], f[w, 1]):
                valid = ins.copy()
                for x1, y1 in zip(f[:, 0], f[:, 1]):
                    if (x1, y1) != (x0, y0): valid &= np.hypot(X - x1, Y - y1) > 1.5
                for j, v in enumerate(per_defect(th2, S, th2s, X, Y, x0, y0, valid, 1.0)): vals[j].append(sg * v)
        if not tilts: return None
        return dict(run=os.path.relpath(run, ROOT + "/runs"), k=k, tilt=np.degrees(np.median(tilts)), n=len(vals[0]), r23=np.degrees(np.nanmedian(vals[0])), r345=np.degrees(np.nanmedian(vals[1])))
    except Exception as ex:
        return dict(run=run, err=repr(ex))
def continuum(job):
    import flatwall_band as FB, flatwall as F, flat_pinning as FP, scipy.ndimage as nd
    a_s, s, lam = job; S0, h, c = F.S0, F.h, F.c; H = 40.0; dt = 0.04; nx, ny = int(round(a_s / h)), int(round(H / h)) + 2
    x = (np.arange(nx) + 0.5) * h; y = (np.arange(ny) - 1) * h; X, Y = np.meshgrid(x, y, indexing="ij")
    thw = np.pi / 2 + np.pi * x / a_s; qw1, qw2 = 0.5 * S0 * np.cos(2 * thw), 0.5 * S0 * np.sin(2 * thw)
    w = np.exp(-np.maximum(Y, 0) / 3.0); th = w * (np.pi / 2 + np.pi * X / a_s); q1, q2 = 0.5 * S0 * np.cos(2 * th), 0.5 * S0 * np.sin(2 * th)
    def bc(q1, q2): q1[:, 0] = qw1; q2[:, 0] = qw2; q1[:, -1] = q1[:, -2]; q2[:, -1] = q2[:, -2]; return q1, q2
    lap = lambda f: (np.roll(f, 1, 0) + np.roll(f, -1, 0) - 2 * f) / h ** 2 + (np.roll(f, 1, 1) + np.roll(f, -1, 1) - 2 * f) / h ** 2
    gx = lambda f: (np.roll(f, -1, 0) - np.roll(f, 1, 0)) / (2 * h); vx = s * np.maximum(Y, 0); Exy = Oxy = 0.5 * s; lam_c = lam * S0
    def step(q1, q2, flow):
        tr = 2 * (q1 ** 2 + q2 ** 2); H1 = lap(q1) + q1 - c * tr * q1; H2 = lap(q2) + q2 - c * tr * q2
        if flow: d1 = -vx * gx(q1) + 2 * Oxy * q2 + H1; d2 = -vx * gx(q2) + lam_c * Exy - 2 * Oxy * q1 + H2
        else: d1, d2 = H1, H2
        return bc(q1 + dt * d1, q2 + dt * d2)
    q1, q2 = bc(q1, q2)
    for _ in range(int(300 / dt)): q1, q2 = step(q1, q2, False)
    for _ in range(int(3000 / dt)): q1, q2 = step(q1, q2, True)
    m = F.measure(q1, q2, X, Y, a_s); S = 2 * np.hypot(q1, q2)
    lab, n = nd.label((S < 0.5 * S0) & (Y > 0.3), structure=np.ones((3, 3))); i = np.argmax(nd.sum(S0 - S, lab, range(1, n + 1))) + 1; mk = lab == i
    ang = np.angle(np.sum(np.exp(2j * np.pi * X[mk] / a_s) * (S0 - S[mk]))); cx = (ang % (2 * np.pi)) * a_s / (2 * np.pi); cy = float(np.sum(Y[mk] * (S0 - S[mk])) / np.sum(S0 - S[mk]))
    # periodic copies: shift x so the defect sits mid-period, tile 3 periods for rings
    Xs, Ys, ths = FP.setup(a_s, cy, cx, hh=h, ymax=H)[:3]
    th2m = np.angle(q1 + 1j * q2); Sm = S / S0
    th2s = np.full(X.shape, np.nan)
    for j in range(1, ny):
        yy = Y[0, j]
        if yy < h / 2 or yy > Ys[0, -1]: continue
        jj = min(max(np.searchsorted(Ys[0], yy), 1), Ys.shape[1] - 1); fr_ = (yy - Ys[0, jj - 1]) / h
        th2s[:, j] = np.angle((1 - fr_) * np.exp(2j * ths[:, jj - 1]) + fr_ * np.exp(2j * ths[:, jj]))
    T3 = lambda A: np.concatenate([A, A, A], axis=0); XX = np.concatenate([X - a_s, X, X + a_s], axis=0); YY = T3(Y)
    valid = (YY > 0) & np.isfinite(T3(th2s)); r = per_defect(T3(th2m), T3(Sm), np.nan_to_num(T3(th2s)), XX, YY, cx, cy, valid, h)
    sg = np.sign(m["tilt"]) or 1.0
    return dict(job=job, tilt=m["tilt"], depth=cy, r23=np.degrees(sg * r[0]), r345=np.degrees(sg * r[1]))
if __name__ == "__main__":
    runs = sorted(r for r in glob.glob(ROOT + "/runs/swirl_crowd/*k*_s*") + glob.glob(ROOT + "/runs/swirl_shear/k*_A*") if os.path.exists(r + "/DONE")) + \
           [ROOT + "/runs/switchon2/dir_k06_alignX/active", ROOT + "/runs/switchon2/dir_k06_rand_s3/active"]
    with Pool(6) as p:
        M = p.map(mpcd, runs); C = p.map(continuum, [(a_s, s, 2.0) for s in (0.03, 0.06) for a_s in (8.0, 12.0, 20.0)])
    L = ["CORE ROTATION (deg; core orientation minus ring orientation, measured minus static; + = core rotated further in the flow sense). Ring 2-3 cells | ring 3-4.5 cells",
         "MPCD (passive imposed vortex unless marked active): run | k | tilt | wall-defect samples | ring 2-3 | ring 3-4.5"]
    for x in M:
        if not x: continue
        if "err" in x: L.append(f"ERROR {x['run']}: {x['err']}"); continue
        L.append(f"{x['run']:34s} | {x['k']:2d} | {x['tilt']:5.1f} | {x['n']:4d} | {x['r23']:+6.1f} | {x['r345']:+6.1f}")
    L.append("CONTINUUM flat wall: a_s | s | tilt | depth | ring 2-3 | ring 3-4.5")
    for x in C: L.append(f"{x['job'][0]:4.0f} | {x['job'][1]:.2f} | {x['tilt']:5.1f} | {x['depth']:.2f} | {x['r23']:+6.1f} | {x['r345']:+6.1f}")
    open("core_rotation.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
