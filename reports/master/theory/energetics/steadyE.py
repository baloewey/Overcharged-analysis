#!/usr/bin/env python3
"""Is the depleted wall layer of elongated cavities ELASTICALLY preferred? (user, 2026-10-12). Frank energy (binding.Ellipse, anisotropic K1 42 / K3 23.1, cores relaxed,
passive) of MEASURED steady-state configurations (second half, one frame per 100 t.u.; frames with total charge = k only) vs two single-defect moves:
  REFILL: each bulk +1/2 (depth > 6) moved to the middle of the largest arclength gap between neighbouring wall +1/2, at depth 1.6 -> dE_refill (min and mean over movers);
  EJECT: each wall +1/2 moved 7 cells inward along its normal -> dE_eject (min and mean).
dE_refill < 0: the wall would like to be refilled (depletion held up by kinetics/activity); > 0: the depleted layer is an elastic minimum for that geometry.
Usage: steadyE.py TAG RUNDIR [...] (env NPROC=4) -> steadyE_<TAG>.npy, steadyE_<TAG>.txt"""
import os, sys, numpy as np
from multiprocessing import Pool
from scipy.spatial import cKDTree
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
import energetics as En  # noqa
bd = En.bd
sys.path.insert(0, os.path.join(os.path.dirname(H), 'cassini_binding'))
from cassini_binding import CassiniAB  # noqa
from energetics_general import EllipseAB  # noqa
def one(run):
    try:
        J, bc, k, t, fr, ev = En.events(run, max_events=0); eps = float(bc.get("aInv", [1])[0] or 1); lam = J["species"][0].get("tumble", 2.0); bd.k = k
        if bc.get("shape", 0) == 1: E = CassiniAB(float(bc["cassiniA"]), float(bc["R"])); eps = round(E.A / E.B, 2)
        else: B_ = float(bc["R"]); E = EllipseAB(B_ * float(bc.get("aInv", [1])[0] or 1), B_)
        inside = (lambda z: E.contains(z)) if hasattr(E, "contains") else (lambda z: (z.real / E.A) ** 2 + (z.imag / E.B) ** 2 < 1)
        sg = 1 if k > 0 else -1; W, s, P = En.wall_polyline(bc); tree = cKDTree(E.W.real[:, None] * [1, 0] + E.W.imag[:, None] * [0, 1])
        t = np.array(t); h = len(fr) // 2; step = max(1, int(round(100 / np.median(np.diff(t))))); out = []
        def en(z, q): return E.energy(list(z), list(q), [None] * len(z))[0]
        for n in range(h, len(fr), step):
            f = fr[n]
            if not len(f): continue
            z = f[:, 0] + 1j * f[:, 1]; q = np.sign(f[:, 2]) * 0.5
            if abs(q.sum() - k) > 1e-6: continue
            d, j = tree.query(f[:, :2]); wall = np.flatnonzero((np.sign(q) == sg) & (d < 4)); bulk = np.flatnonzero((np.sign(q) == sg) & (d > 6))
            if len(wall) < 2: continue
            E0 = en(z, q)
            sw = np.sort([E.s[j[i]] for i in wall]); gaps = np.diff(np.r_[sw, sw[0] + E.P]); g = int(np.argmax(gaps)); smid = (sw[g] + gaps[g] / 2) % E.P
            im = int(np.searchsorted(E.s, smid) % len(E.s)); target = E.W[im] - 1.6 * E.Nn[im]
            dR = []
            for i in bulk:
                z2 = z.copy(); z2[i] = target; dR.append(en(z2, q) - E0)
            dEj = []
            for i in wall:
                z2 = z.copy(); z2[i] = z[i] - 7.0 * E.Nn[j[i]]
                if inside(z2[i]) and np.min(np.abs(E.W - z2[i])) > 2.0: dEj.append(en(z2, q) - E0)
            out.append((lam, eps, len(wall), len(bulk), gaps.max(), min(dR) if dR else np.nan, np.mean(dR) if dR else np.nan, min(dEj) if dEj else np.nan, np.mean(dEj) if dEj else np.nan))
        return run, out
    except Exception as ex:
        return run, [("ERR", str(ex))]
if __name__ == "__main__":
    tag = sys.argv[1]; runs = [os.path.abspath(r) for r in sys.argv[2:]]
    with Pool(int(os.environ.get("NPROC", 4))) as pool: res = pool.map(one, runs, chunksize=1)
    rows = [r + (run,) for run, o in res for r in o if r[0] != "ERR"]; errs = [(run, o) for run, o in res if o and o[0][0] == "ERR"]
    np.save(os.path.join(H, f"steadyE_{tag}.npy"), np.array(rows, dtype=object), allow_pickle=True)
    a = np.array([r[:9] for r in rows], float)
    L = ["GROUPED: lambda eps | frames | <n_wall> <n_bulk+> | largest gap (cells) | dE_refill min / mean | dE_eject min / mean  (kBT; median over frames)"]
    for key in sorted(set(map(tuple, a[:, :2]))):
        S = a[(a[:, 0] == key[0]) & (a[:, 1] == key[1])]; md = lambda c: np.nanmedian(S[:, c])
        L.append(f"  λ {key[0]:3.1f} ε {key[1]:3.1f} | {len(S)} | {S[:, 2].mean():.2f} {S[:, 3].mean():.2f} | {md(4):.1f} | {md(5):+.1f} / {md(6):+.1f} | {md(7):+.1f} / {md(8):+.1f}  (P(dE_refill_min<0) {np.mean(S[:, 5] < 0):.2f})")
    L += [f"ERR {r}: {o[0][1]}" for r, o in errs]
    print("\n".join(L)); open(os.path.join(H, f"steadyE_{tag}.txt"), "w").write("\n".join(L) + "\n")
