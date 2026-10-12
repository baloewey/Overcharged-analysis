#!/usr/bin/env python3
"""All-channel turnover of the wall layer (follow-up to refill.py: depleting escapes are a minority of losses). Second half of each run, frames every 2 t.u.
EXIT: a wall +1/2 (< 4 cells) whose chain reaches depth > 6 within 10 frames and stays > 4 for 10 more frames (any residence time).
ENTRY: a bulk +1/2 (> 6) whose chain reaches < 4 within 10 frames and stays < 4 for 10 more frames.
ANNIHILATION at the wall: a wall +1/2 with no successor within 2.5 cells in the next frame while a -1/2 within 3 cells also vanishes.
Per run: rates per wall defect (exit, annihilation) and per vacancy (entry), the depleting subset (energetics.events), and the predicted balance n*.
Usage: turnover.py TAG RUNDIR [...] -> turnover_<TAG>.txt"""
import os, sys, numpy as np
from multiprocessing import Pool
from scipy.spatial import cKDTree
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
import energetics as En  # noqa
def one(run):
    J, bc, k, t, fr, ev = En.events(run, max_events=10 ** 6); lam = J["species"][0].get("tumble", 2.0); eps = float(bc.get("aInv", [1])[0] or 1)
    N0 = 2 * abs(k) - 2; sg = 1 if k > 0 else -1; W, s, P = En.wall_polyline(bc); tree = cKDTree(W); t = np.array(t)
    D = [tree.query(f[:, :2])[0] if len(f) else np.zeros(0) for f in fr]; h = len(fr) // 2; T = t[-1] - t[h]
    nw = np.array([((np.sign(f[:, 2]) == sg) & (D[n] < 4)).sum() if len(f) else 0 for n, f in enumerate(fr)], float)
    def nxt(m, c, sgn=sg, lim=2.5):
        g = fr[m]; gp = g[np.sign(g[:, 2]) == sgn, :2] if len(g) else np.zeros((0, 2))
        if not len(gp): return None
        d = np.hypot(*(gp - c).T); kk = np.argmin(d); return gp[kk] if d[kk] < lim else None
    def crosses(n, c, target_out):
        for m in range(n + 1, n + 11):
            c = nxt(m, c)
            if c is None: return False
            d = tree.query(c)[0]
            if (d > 6 if target_out else d < 4):
                for mm in range(m + 1, m + 11):
                    c = nxt(mm, c)
                    if c is None: return False
                    if (tree.query(c)[0] < 4) if target_out else (tree.query(c)[0] > 4): return False
                return True
        return False
    nex = nen = nan = 0
    for n in range(h, len(fr) - 22, 1):
        f = fr[n]
        if not len(f): continue
        wi = np.flatnonzero((np.sign(f[:, 2]) == sg) & (D[n] < 4)); bi = np.flatnonzero((np.sign(f[:, 2]) == sg) & (D[n] > 6))
        # count each crossing once: require the defect to be on its side in frame n but crossing starts at n (previous frame also on that side)
        for i in wi:
            if crosses(n, f[i, :2], True):
                pv = nxt(n - 1, f[i, :2]) if n > 0 else None
                if pv is not None and tree.query(pv)[0] < 4 and not crosses(n - 1, pv, True): nex += 1
            elif nxt(n + 1, f[i, :2]) is None:
                mi = np.flatnonzero((np.sign(f[:, 2]) == -sg) & (np.hypot(*(f[:, :2] - f[i, :2]).T) < 3))
                if len(mi) and all(nxt(n + 1, f[j, :2], -sg) is None for j in mi): nan += 1
        for i in bi:
            if crosses(n, f[i, :2], False):
                pv = nxt(n - 1, f[i, :2]) if n > 0 else None
                if pv is not None and tree.query(pv)[0] > 4 and not crosses(n - 1, pv, False): nen += 1
    nd = sum(1 for e in ev if t[e["n"]] >= t[h]); nm = nw[h:].mean(); vac = max(N0 - nm, 1e-9)
    return (lam, eps, nm, nex, nen, nan, nd, 1000 * nex / (nm * T), 1000 * nan / (nm * T), 1000 * nen / (vac * T), 1000 * nd / (nm * T), N0)
if __name__ == "__main__":
    tag = sys.argv[1]; runs = [os.path.abspath(r) for r in sys.argv[2:]]
    with Pool(min(6, len(runs))) as pool: res = pool.map(one, runs, chunksize=1)
    L = ["GROUPED: lambda eps | runs | <n_wall> | exits entries annih. depleting | k_exit k_annih (per wall defect /1000 t.u.) | k_entry (per vacancy) | depleting share of exits | n* = N0 k_en/(k_ex + k_an + k_en)"]
    for key in sorted(set((r[0], r[1]) for r in res)):
        S = np.array([r[2:] for r in res if (r[0], r[1]) == key], float); m = S.mean(0); tot = S.sum(0); N0 = m[-1]
        L.append(f"  λ {key[0]:3.1f} ε {key[1]:3.1f} | {len(S)} | {m[0]:5.2f} | {int(tot[1])} {int(tot[2])} {int(tot[3])} {int(tot[4])} | {m[5]:.2f} {m[6]:.2f} | {m[7]:.2f} | {tot[4] / max(tot[1], 1):.2f} | {N0 * m[7] / (m[5] + m[6] + m[7]):.2f}")
    print("\n".join(L)); open(os.path.join(H, f"turnover_{tag}.txt"), "w").write("\n".join(L) + "\n")
