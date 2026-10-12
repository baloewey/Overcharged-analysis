#!/usr/bin/env python3
"""Loss vs refill of the wall layer (user, 2026-10-12: with higher lambda escapes are rarer, but restoring events may get rarer too -> same stationary state).
Second half of each run. LOSS events = depleting escapes (energetics.events). REFILL events = a bulk +1/2 (depth > 6) that enters the band (< 4 cells) and then STAYS there
for >= 40 t.u. (20 frames, chain followed with gap bridging as in events), while the smoothed wall count is <= 2k-2 (fills a vacancy) and rises afterwards.
Rates: k_loss = losses / (<n_wall> T) per wall defect; k_refill = refills / (<N0 - n_wall> T) per vacancy (both per 1000 t.u.). Balance: n* = N0 k_refill / (k_loss + k_refill).
Usage: refill.py TAG RUNDIR [...] -> refill_<TAG>.txt"""
import os, sys, numpy as np
from multiprocessing import Pool
from scipy.spatial import cKDTree
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
import energetics as En  # noqa
def one(run):
    J, bc, k, t, fr, ev = En.events(run, max_events=10 ** 6); lam = J["species"][0].get("tumble", 2.0); eps = float(bc.get("aInv", [1])[0] or 1)
    N0 = 2 * abs(k) - 2; sg = 1 if k > 0 else -1; W, s, P = En.wall_polyline(bc); tree = cKDTree(W); t = np.array(t)
    D = [tree.query(f[:, :2])[0] if len(f) else np.zeros(0) for f in fr]
    nw = np.array([((np.sign(f[:, 2]) == sg) & (D[n] < 4)).sum() if len(f) else 0 for n, f in enumerate(fr)], float); nws = np.convolve(nw, np.ones(11) / 11, mode="same")
    h = len(fr) // 2; T = t[-1] - t[h]; nref = 0
    def nxt(m, c, lim=2.5):
        g = fr[m]; gp = g[np.sign(g[:, 2]) == sg, :2] if len(g) else np.zeros((0, 2))
        if not len(gp): return None
        d = np.hypot(*(gp - c).T); kk = np.argmin(d); return gp[kk] if d[kk] < lim else None
    for n in range(max(h, 11), len(fr) - 26):
        f = fr[n]
        if not len(f): continue
        for i in np.flatnonzero((np.sign(f[:, 2]) == sg) & (D[n] > 6)):
            # does it reach the band within 10 frames, then stay 20 frames?
            c = f[i, :2]; ok = False; m = n
            for m in range(n + 1, n + 11):
                c = nxt(m, c)
                if c is None: break
                if tree.query(c)[0] < 4: ok = True; break
            if not ok: continue
            stay = True; c2 = c
            for mm in range(m + 1, m + 21):
                if mm >= len(fr): stay = False; break
                c2 = nxt(mm, c2)
                if c2 is None or tree.query(c2)[0] > 4: stay = False; break
            if stay and nws[max(m - 5, 0)] <= N0 + 0.25 and m + 25 < len(nws) and nws[m + 5:m + 26].min() > nws[max(m - 5, 0)] + 0.25: nref += 1
    nloss = sum(1 for e in ev if t[e["n"]] >= t[h]); nmean = nw[h:].mean(); vac = max(N0 - nmean, 1e-9)
    return (os.path.basename(os.path.dirname(run)) + "/" + os.path.basename(run), lam, eps, nmean, nloss, nref, 1000 * nloss / (nmean * T), 1000 * nref / (vac * T), N0)
if __name__ == "__main__":
    tag = sys.argv[1]; runs = [os.path.abspath(r) for r in sys.argv[2:]]
    with Pool(min(6, len(runs))) as pool: res = pool.map(one, runs, chunksize=1)
    L = ["GROUPED: lambda eps | runs | <n_wall> | losses refills | k_loss (per wall defect /1000 t.u.) | k_refill (per vacancy /1000 t.u.) | n* = N0 k_r/(k_l+k_r)"]
    for key in sorted(set((r[1], r[2]) for r in res)):
        S = [r for r in res if (r[1], r[2]) == key]; n = np.mean([r[3] for r in S]); nl = sum(r[4] for r in S); nr = sum(r[5] for r in S)
        kl = np.mean([r[6] for r in S]); kr = np.mean([r[7] for r in S]); N0 = S[0][8]
        L.append(f"  λ {key[0]:3.1f} ε {key[1]:3.1f} | {len(S)} | {n:5.2f} | {nl:3d} {nr:3d} | {kl:.3f} | {kr:.3f} | {N0 * kr / max(kl + kr, 1e-9):.2f}")
    print("\n".join(L)); open(os.path.join(H, f"refill_{tag}.txt"), "w").write("\n".join(L) + "\n")
