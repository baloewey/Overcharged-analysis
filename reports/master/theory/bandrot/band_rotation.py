#!/usr/bin/env python3
"""Flow-induced rotation of the director near the wall layer, MPCD passive imposed-vortex runs (user, 2026-10-22).
For each director frame (second half; frames whose tracked defects sum to the wall charge k): the exact static one-constant field for the tracked
defect positions (ejection/tilted_structures.fields: all defects + images, wall anchoring theta = k phi as in the runs) is subtracted from the
measured cell director: dtheta = theta_meas - theta_static (mod pi, in (-90, 90] deg), signed by the circulation. With no flow and K1 = K3 this
is zero up to noise; with flow it is the rotation the flow adds on top of the static texture of the same defects.
Profiles against distance from the wall u = R - r (0.5-cell bins), cells farther than 1.5 cells from every core, for three azimuthal classes:
'at defects' (within 0.15 spacing of a wall +1/2), 'mid-gap' (within 0.15 spacing of the midpoint between neighbours), and all.
Per run also: wall shear s (tiltcrowd/wallshear), median wall-defect depth d, tilt; the band estimate chi gamma1 s d^2/(4K) (vorticity acting alone
on a band of width d held at the wall), chi 0.5, gamma1 60, K 44. -> band_rotation.txt, band_rotation.npz"""
import os, sys, glob, gzip, json, warnings, numpy as np
from multiprocessing import Pool
warnings.filterwarnings("ignore")
ROOT = "/home/baloewe/Physics/Overcharged"; T_ = ROOT + "/reports/master/theory"
sys.path.insert(0, T_ + "/tiltcrowd"); sys.path.insert(0, ROOT + "/runs/cavity")
import bulk_spiral as BS, wallshear as WS  # noqa
from edge_criteria import read_with_angle  # noqa
import importlib.util
spec = importlib.util.spec_from_file_location("ts", T_ + "/ejection/tilted_structures.py"); TS = importlib.util.module_from_spec(spec)
exec(compile(open(T_ + "/ejection/tilted_structures.py").read().split('\nif __name__ == "__main__":')[0], "ts", "exec"), TS.__dict__)
UB = np.arange(0, 10.01, 0.5); CHI, G1, K = 0.5, 60.0, 44.0
def one(run):
    try:
        bc = json.load(open(run + "/input.json"))["BC"][0]; R = bc["R"]; k = bc["customAnch"]; Q = np.array(bc["Q"][:2], float); N = int(2 * k - 2)
        fn = run + "/directorfield.dat.gz" if os.path.exists(run + "/directorfield.dat.gz") else run + "/directorfield.dat"
        op = gzip.open(fn, "rt") if fn.endswith(".gz") else open(fn)
        with op as fh: D = np.loadtxt((l for l in fh if len(l.split()) == 8 and l.split()[0][0].isdigit()))
        td = np.unique(D[:, 0]); td = td[td >= td[0] + 0.5 * (td[-1] - td[0])]
        t, fr = read_with_angle(run, Q, 0.0); t = np.array(t); n = int(D[:, 1].max()) + 1; m = int(D[:, 2].max()) + 1
        X, Y = np.meshgrid(np.arange(n) + 0.5 - Q[0], np.arange(m) + 0.5 - Q[1], indexing="ij"); r = np.hypot(X, Y); u = R - r; phi = np.arctan2(Y, X); ins = r < R
        acc = {c: [np.zeros(len(UB) - 1), np.zeros(len(UB) - 1)] for c in ("at", "mid", "all")}; tilts, depths, nfr = [], [], 0
        for tt in td:
            i = int(np.argmin(np.abs(t - tt)))
            if abs(t[i] - tt) > 1e-6: continue
            f = fr[i]
            if abs(f[:, 2].sum() - k) > 1e-6: continue
            rr = np.hypot(f[:, 0], f[:, 1]); w = (f[:, 2] > 0) & (R - rr < 4)
            if w.sum() != N: continue
            zs = list(f[:, 0] + 1j * f[:, 1]); qs = list(f[:, 2]); e2, _, _ = TS.fields(X, Y, zs, qs, R)
            d = D[D[:, 0] == tt]; th = np.full((n, m), np.nan); th[d[:, 1].astype(int), d[:, 2].astype(int)] = np.arctan2(d[:, 5], d[:, 4])
            dth = np.angle(np.exp(1j * (2 * th - e2))) / 2
            tau = np.angle(np.mean(np.exp(1j * (f[w, 3] - np.arctan2(f[w, 1], f[w, 0]))))); sg = np.sign(tau) or 1.0; tilts.append(abs(tau)); depths.append(np.median(R - rr[w]))
            far = np.ones((n, m), bool)
            for z in zs: far &= np.abs(X + 1j * Y - z) > 1.5
            pw = np.sort(np.arctan2(f[w, 1], f[w, 0])); gaps = np.diff(np.r_[pw, pw[0] + 2 * np.pi]); mids = pw + gaps / 2; sp = 2 * np.pi / N
            dat = np.min(np.abs(np.angle(np.exp(1j * (phi[..., None] - pw[None, None, :])))), axis=-1); dmid = np.min(np.abs(np.angle(np.exp(1j * (phi[..., None] - mids[None, None, :])))), axis=-1)
            for cls, mk in (("at", dat < 0.15 * sp), ("mid", dmid < 0.15 * sp), ("all", np.ones_like(far))):
                sel = ins & far & mk & np.isfinite(dth)
                h, _ = np.histogram(u[sel], UB, weights=sg * dth[sel]); c_, _ = np.histogram(u[sel], UB); acc[cls][0] += h; acc[cls][1] += c_
            nfr += 1
        if nfr < 3: return None
        prof = {cls: np.degrees(v[0] / np.maximum(v[1], 1)) for cls, v in acc.items()}
        ws = WS.one(run); s = ws["e_band"] + ws["w_band"]; dm = float(np.median(depths))
        return dict(run=os.path.relpath(run, ROOT + "/runs"), k=k, a=2 * np.pi * (R - dm) / N, s=s, d=dm, tilt=np.degrees(np.median(tilts)), nfr=nfr,
                    est=np.degrees(CHI * G1 * s * dm ** 2 / (4 * K)), **{f"p_{c}": v for c, v in prof.items()})
    except Exception as ex:
        return dict(run=run, err=repr(ex))
if __name__ == "__main__":
    runs = sorted(r for r in glob.glob(ROOT + "/runs/swirl_crowd/*k*_s*") + glob.glob(ROOT + "/runs/swirl_shear/k*_A*") if os.path.exists(r + "/DONE"))
    with Pool(6) as p: res = p.map(one, runs)
    ok = [x for x in res if x and "err" not in x]; bad = [x for x in res if x and "err" in x]
    uc = 0.5 * (UB[1:] + UB[:-1]); L = [f"{len(ok)} runs ({len(bad)} errors, {sum(1 for x in res if x is None)} with < 3 usable frames). dtheta in deg, signed by circulation (positive = rotated further in the flow sense)."]
    L.append("run | k | spacing | depth d | wall shear s | tilt | frames | band estimate | dtheta (all) at u = 0.75, 1.25, 1.75, 2.25 (band) | 3.25, 4.75, 7.25 (beyond) | mid-gap u 0.75-1.75 | at-defect u 2.25-3.25")
    for x in ok:
        pa, pm, pt = x["p_all"], x["p_mid"], x["p_at"]; ix = lambda v: int(np.argmin(np.abs(uc - v)))
        L.append(f"{x['run']:24s} | {x['k']:2d} | {x['a']:5.1f} | {x['d']:.2f} | {x['s']:.3f} | {x['tilt']:5.1f} | {x['nfr']:2d} | {x['est']:+5.1f} | " +
                 " ".join(f"{pa[ix(v)]:+5.1f}" for v in (0.75, 1.25, 1.75, 2.25)) + " | " + " ".join(f"{pa[ix(v)]:+5.1f}" for v in (3.25, 4.75, 7.25)) +
                 f" | {np.mean(pm[ix(0.75):ix(1.75)+1]):+5.1f} | {np.mean(pt[ix(2.25):ix(3.25)+1]):+5.1f}")
    for x in bad: L.append(f"ERROR {x['run']}: {x['err']}")
    np.savez("band_rotation.npz", uc=uc, **{f"{x['run']}|{c}": x[f"p_{c}"] for x in ok for c in ("all", "mid", "at")})
    open("band_rotation.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
