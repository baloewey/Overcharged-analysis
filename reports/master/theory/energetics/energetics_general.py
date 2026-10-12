#!/usr/bin/env python3
"""Generalisation of energetics.py / energetics_ts.py to any ellipse or circle (area, k) of the shape sweep (energy only; the sweep's flow is every 100 t.u.).
Ellipse geometry from the run input: circle (aInv 1,1) A = B = R; ellipse aInv = [eps, 1, 1], R = B, A = eps * B. The binding.py Ellipse is rebuilt for that A, B and the module's
winding k is set to the run's k before each energy. Lags (frames of 2 t.u.): -20, -10, -5, -2, 0, +2; controls at -20, -10, 0. Up to 10 depleting escapes per run.
Usage: energetics_general.py TAG RUNDIR [...] -> etsg_<TAG>.npy (rows: run, kind, lag, F0, E2, E6, cos tau, aspect, R_eq, region 0 tip / 1 mid / 2 flank of the escape site or control)"""
import json, os, sys
import numpy as np
from multiprocessing import Pool
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
import energetics as En  # noqa
from tipcount import reg  # noqa  (0 tip, 1 mid, 2 flank)
bd = En.bd
LAGS = (-20, -10, -5, -2, 0, 2)
class EllipseAB(bd.Ellipse):
    def __init__(self, A, B):
        self.A, self.B = A, B; h = bd.h
        import scipy.sparse as sp
        x = np.arange(-A - 2, A + 2 + h / 2, h); y = np.arange(-B - 2, B + 2 + h / 2, h); X, Y = np.meshgrid(x, y, indexing="ij"); self.X, self.Y = X, Y; self.Z = X + 1j * Y
        nx, ny = len(x), len(y); Dx1 = sp.diags([-1, 1], [-1, 1], shape=(nx, nx)) / (2 * h); Dy1 = sp.diags([-1, 1], [-1, 1], shape=(ny, ny)) / (2 * h)
        self.Dx = sp.kron(Dx1, sp.identity(ny), format="csr"); self.Dy = sp.kron(sp.identity(nx), Dy1, format="csr")
        rho = np.sqrt((X / A) ** 2 + (Y / B) ** 2); self.inside = rho < 1; self.t, self.W, self.s, self.P, self.Nn = bd.wall_param(A, B)
        out = ~self.inside & (rho < 1 + 1.6 / min(A, B)); self.band = out; idx = np.flatnonzero(out.ravel()); zb = self.Z.ravel()[idx]
        self.band_idx, self.band_near = idx, np.array([np.argmin(np.abs(self.W - zz)) for zz in zb])
def one(run):
    try:
        J, bc, k, t, fr, ev = En.events(run, max_events=10); bd.k = k
        ainv = bc.get("aInv", [1, 1, 1]); B = float(bc["R"]); A = B * float(ainv[0]) if ainv[0] else B
        E = EllipseAB(A, B); sg = 1 if k > 0 else -1; out = []
        from scipy.spatial import cKDTree
        W, s, P = En.wall_polyline(bc); tree = cKDTree(W); Amax = np.abs(W[:, 0]).max(); rng = np.random.default_rng(13)
        req = np.sqrt(A * B) / 22.21; asp = A / B
        def pcos(frame, p):
            zs = frame[:, 0] + 1j * frame[:, 1]; i = int(np.argmin(np.abs(zs - (p[0] + 1j * p[1])))); nrm = E.Nn[np.argmin(np.abs(E.W - zs[i]))]
            return np.cos(np.angle(np.exp(1j * (frame[i, 3] - np.angle(nrm)))))
        for e in ev:
            n = e["n"]; chain = dict(e["chain"]); c = chain[n]
            for m in range(n + 1, n + 3):
                h = fr[m]; hp = h[np.sign(h[:, 2]) == sg, :2] if len(h) else np.zeros((0, 2))
                if not len(hp): break
                d2 = np.hypot(*(hp - c).T); kk = np.argmin(d2)
                if d2[kk] > 2.5: break
                c = hp[kk]; chain[m] = c
            for lag in LAGS:
                m = n + lag
                if m not in chain: continue
                f = fr[m]; L = En.landscape(E, f, chain[m], k)
                if L is None: continue
                out.append((os.path.relpath(run, En.ROOT), "esc", lag, L["F0"], L["prof"][3], L["E6"], pcos(f, chain[m]), asp, req, reg(chain[n], A, B)))
                if lag in (-20, -10, 0):
                    d, j = tree.query(f[:, :2]); res = [jj for jj in np.flatnonzero((np.sign(f[:, 2]) == sg) & (d < 4)) if np.hypot(*(f[jj, :2] - chain[m])) > 3 and abs(abs(f[jj, 0]) - abs(chain[m][0])) / Amax < 0.15]
                    if res:
                        pc = f[rng.choice(res), :2]; Lc = En.landscape(E, f, pc, k)
                        if Lc is not None: out.append((os.path.relpath(run, En.ROOT), "ctl", lag, Lc["F0"], Lc["prof"][3], Lc["E6"], pcos(f, pc), asp, req, reg(pc, A, B)))
        return out
    except Exception as ex:
        return [("ERR " + run + " " + str(ex),)]
if __name__ == "__main__":
    tag = sys.argv[1]; runs = [os.path.abspath(r) for r in sys.argv[2:]]
    with Pool(8) as pool: res = pool.map(one, runs, chunksize=1)
    rows = sum(res, []); errs = [r for r in rows if len(r) == 1]; rows = [r for r in rows if len(r) > 1]
    np.save(os.path.join(H, f"etsg_{tag}.npy"), np.array(rows, dtype=object), allow_pickle=True); print("done", tag, len(rows), "errors", len(errs)); [print(e[0]) for e in errs[:5]]
