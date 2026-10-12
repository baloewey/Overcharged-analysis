#!/usr/bin/env python3
"""Self-limiting bulk charge? (user, 2026-10-12): more bulk charge pushes wall +1/2 closer to the wall -> weaker escape (less jet transmitted, stronger elastic hold) ->
the exit rate falls as the bulk charges up (and entries rise) -> steady state. Within-run test. Second half, frames every 2 t.u.
Per frame: n_wall (+1/2 < 4 cells), bulk charge Q_b = sum of charges deeper than 4 cells, mean depth of wall +1/2. Exit and entry events as in turnover.py, tagged with the frame's
Q_b. Everything binned by dQ = Q_b - <Q_b>_run (fluctuation within the run). Hazards: exits per wall defect per 1000 t.u., entries per vacancy (2k-2 - n_wall) per 1000 t.u.
Usage: hazard.py TAG RUNDIR [...] (env NPROC) -> hazard_<TAG>.npz, hazard_<TAG>.txt"""
import os, sys, numpy as np
from multiprocessing import Pool
from scipy.spatial import cKDTree
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
import energetics as En  # noqa
def one(run):
    try:
        J, bc, k, t, fr, ev = En.events(run, max_events=0); eps = float(bc.get("aInv", [1])[0] or 1); lam = J["species"][0].get("tumble", 2.0)
        N0 = 2 * abs(k) - 2; sg = 1 if k > 0 else -1; W, s, P = En.wall_polyline(bc); tree = cKDTree(W); t = np.array(t); h = len(fr) // 2
        D = [tree.query(f[:, :2])[0] if len(f) else np.zeros(0) for f in fr]
        def nxt(m, c, lim=2.5):
            g = fr[m]; gp = g[np.sign(g[:, 2]) == sg, :2] if len(g) else np.zeros((0, 2))
            if not len(gp): return None
            d = np.hypot(*(gp - c).T); kk = np.argmin(d); return gp[kk] if d[kk] < lim else None
        def crosses(n, c, out_):
            for m in range(n + 1, min(n + 11, len(fr))):
                c = nxt(m, c)
                if c is None: return False
                d = tree.query(c)[0]
                if (d > 6 if out_ else d < 4):
                    for mm in range(m + 1, min(m + 11, len(fr))):
                        c = nxt(mm, c)
                        if c is None or ((tree.query(c)[0] < 4) if out_ else (tree.query(c)[0] > 4)): return False
                    return True
            return False
        fq = []; evs = []
        for n in range(h, len(fr) - 22):
            f = fr[n]
            if not len(f): continue
            w = (np.sign(f[:, 2]) == sg) & (D[n] < 4); nw = w.sum(); qb = 0.5 * np.sign(f[D[n] >= 4, 2]).sum(); dep = D[n][w].mean() if nw else np.nan
            fq.append((n, nw, qb, dep))
            for i in np.flatnonzero(w):
                if crosses(n, f[i, :2], True):
                    pv = nxt(n - 1, f[i, :2])
                    if pv is not None and tree.query(pv)[0] < 4 and not crosses(n - 1, pv, True): evs.append((n, 1))
            for i in np.flatnonzero((np.sign(f[:, 2]) == sg) & (D[n] > 6)):
                if crosses(n, f[i, :2], False):
                    pv = nxt(n - 1, f[i, :2])
                    if pv is not None and tree.query(pv)[0] > 4 and not crosses(n - 1, pv, False): evs.append((n, -1))
        return dict(run=run, eps=eps, lam=lam, k=k, N0=N0, dt=float(np.median(np.diff(t))), fq=np.array(fq, float), ev=np.array(evs, int).reshape(-1, 2))
    except Exception as ex:
        return dict(run=run, err=str(ex))
if __name__ == "__main__":
    tag = sys.argv[1]; runs = [os.path.abspath(r) for r in sys.argv[2:]]
    with Pool(int(os.environ.get("NPROC", 6))) as pool: res = pool.map(one, runs, chunksize=1)
    ok = [r for r in res if "err" not in r]; np.save(os.path.join(H, f"hazard_{tag}.npy"), np.array(ok, dtype=object), allow_pickle=True)
    L = [f"{len(ok)} runs ({len(res) - len(ok)} errors). Within-run fluctuation of bulk charge dQ = Q_b - run mean (bins of 0.5)."]
    for grp, sel in (("circle-like (eps < 1.3)", lambda r: r["eps"] < 1.3), ("elongated (eps >= 1.3)", lambda r: r["eps"] >= 1.3)):
        S = [r for r in ok if sel(r) and len(r["fq"])]
        if not S: continue
        B = np.arange(-2.25, 2.3, 0.5); expo_w = np.zeros(len(B) - 1); expo_v = np.zeros(len(B) - 1); nex = np.zeros(len(B) - 1); nen = np.zeros(len(B) - 1); dsum = np.zeros(len(B) - 1); dn = np.zeros(len(B) - 1)
        for r in S:
            fq = r["fq"]; dq = fq[:, 2] - fq[:, 2].mean(); b = np.digitize(dq, B) - 1; idx = {int(n): j for j, n in enumerate(fq[:, 0])}
            for j in range(len(fq)):
                if 0 <= b[j] < len(B) - 1:
                    expo_w[b[j]] += fq[j, 1] * r["dt"]; expo_v[b[j]] += max(r["N0"] - fq[j, 1], 0) * r["dt"]
                    if np.isfinite(fq[j, 3]): dsum[b[j]] += fq[j, 3] - np.nanmean(fq[:, 3]); dn[b[j]] += 1
            for n, typ in r["ev"]:
                j = idx.get(int(n))
                if j is None or not (0 <= b[j] < len(B) - 1): continue
                if typ == 1: nex[b[j]] += 1
                else: nen[b[j]] += 1
        L.append(f"\n{grp}: {len(S)} runs, {int(nex.sum())} exits, {int(nen.sum())} entries")
        L.append("  dQ_bulk | frames | wall depth - run mean | exit hazard per wall defect /1000 t.u. (n) | entry hazard per vacancy /1000 t.u. (n)")
        for i in range(len(B) - 1):
            if dn[i] < 50: continue
            L.append(f"  {0.5 * (B[i] + B[i + 1]):+5.2f} | {int(dn[i]):6d} | {dsum[i] / dn[i]:+.3f} | {1000 * nex[i] / max(expo_w[i], 1e-9):6.2f} ± {1000 * np.sqrt(nex[i]) / max(expo_w[i], 1e-9):.2f} ({int(nex[i])}) | {1000 * nen[i] / max(expo_v[i], 1e-9):6.2f} ± {1000 * np.sqrt(nen[i]) / max(expo_v[i], 1e-9):.2f} ({int(nen[i])})")
    print("\n".join(L)); open(os.path.join(H, f"hazard_{tag}.txt"), "w").write("\n".join(L) + "\n")
