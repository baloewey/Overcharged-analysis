#!/usr/bin/env python3
"""Free-defect mobility in the bulk from +-1/2 pair annihilation (runs/rotvisc/quench_L96_s*: passive 96^2 periodic box quenched from random orientations,
defects every 1 t.u.). Annihilation event: a +1/2 and a -1/2 closer than 3 cells at frame n, neither present (no defect of that sign within 2 cells) at frame n+1.
Each partner is followed backwards by nearest-neighbour linking (< 2 cells per frame) while the pair stays the mutual nearest opposite-charge pair.
Passive theory: d(r^2)/dt = -(mu_+ + mu_-) pi K  (attraction pi K / (2 r) for charges +-1/2, one-constant K). Fit r^2 = s (t_a - t) over 3 < r < 12 cells
-> mean mobility mu = s / (2 pi K), K = 32.5 kBT (mean of K1, K3). Also the asymmetry: displacement of each partner towards the annihilation point.
Usage: pair_mobility.py RUN [RUN ...] -> pair_mobility.txt"""
import sys, numpy as np
from scipy.spatial import cKDTree
K = 32.5
def read(run):
    fr, cur, t = [], None, None
    for line in open(run + "/defects.dat") if not run.endswith(".gz") else []:
        p = line.split()
        if len(p) == 2 and p[0].replace(".", "").isdigit():
            if cur is not None: fr.append((t, np.array(cur).reshape(-1, 3)))
            t = float(p[0]); cur = []
        elif len(p) == 4 and cur is not None:
            try: cur.append([float(p[0]), float(p[1]), float(p[2])])
            except ValueError: pass
    if cur is not None: fr.append((t, np.array(cur).reshape(-1, 3)))
    return fr
def pdist(a, b, L): d = np.abs(a - b); d = np.minimum(d, L - d); return np.hypot(d[..., 0], d[..., 1])
def analyse(run, L=96.0):
    import gzip, os
    if not os.path.exists(run + "/defects.dat") and os.path.exists(run + "/defects.dat.gz"):
        txt = gzip.open(run + "/defects.dat.gz", "rt").read(); open("/tmp/claude-1000/_defects_tmp.dat", "w").write(txt); fr = read_text(txt)
    else: fr = read(run)
    T = np.array([t for t, _ in fr]); taus, r2s, n_ev = [], [], 0
    for n in range(10, len(fr) - 1):
        t, f = fr[n]; t1, g = fr[n + 1]
        if not len(f): continue
        P, M = f[f[:, 2] > 0, :2], f[f[:, 2] < 0, :2]
        if not len(P) or not len(M): continue
        d = pdist(P[:, None, :], M[None, :, :], L); i, j = np.unravel_index(np.argmin(d), d.shape) if d.size else (None, None)
        # all close pairs at this frame
        for i in range(len(P)):
            j = np.argmin(d[i]); 
            if d[i, j] > 3 or np.argmin(d[:, j]) != i: continue
            gp, gm = g[g[:, 2] > 0, :2], g[g[:, 2] < 0, :2]
            if (len(gp) and pdist(gp, P[i], L).min() < 2) or (len(gm) and pdist(gm, M[j], L).min() < 2): continue
            n_ev += 1; p, m = P[i].copy(), M[j].copy(); ta = t + 0.5
            for k in range(n, max(n - 300, 0), -1):
                tk, fk = fr[k]; Pk, Mk = fk[fk[:, 2] > 0, :2], fk[fk[:, 2] < 0, :2]
                if not len(Pk) or not len(Mk): break
                dp = pdist(Pk, p, L); dm = pdist(Mk, m, L); ip, im = np.argmin(dp), np.argmin(dm)
                if dp[ip] > 2 or dm[im] > 2: break
                p, m = Pk[ip], Mk[im]; r = pdist(p, m, L)
                # still mutual nearest opposite charges?
                if pdist(Mk, p, L).min() < r - 1e-9 or pdist(Pk, m, L).min() < r - 1e-9: break
                taus.append(ta - tk); r2s.append(r * r)
                if r > 15: break
    taus, r2s = np.array(taus), np.array(r2s); m = (r2s > 9) & (r2s < 144)
    s = np.sum(taus[m] * r2s[m]) / np.sum(taus[m] ** 2) if m.any() else np.nan
    bins = [(3, 5), (5, 8), (8, 12)]; local = []
    for lo, hi in bins:
        k = (r2s > lo * lo) & (r2s < hi * hi); local.append(np.sum(taus[k] * r2s[k]) / np.sum(taus[k] ** 2) if k.any() else np.nan)
    return n_ev, s, local, len(taus)
def read_text(txt):
    fr, cur, t = [], None, None
    for line in txt.splitlines():
        p = line.split()
        if len(p) == 2 and p[0].replace(".", "").isdigit():
            if cur is not None: fr.append((t, np.array(cur).reshape(-1, 3)))
            t = float(p[0]); cur = []
        elif len(p) == 4 and cur is not None:
            try: cur.append([float(p[0]), float(p[1]), float(p[2])])
            except ValueError: pass
    if cur is not None: fr.append((t, np.array(cur).reshape(-1, 3)))
    return fr
if __name__ == "__main__":
    out = []
    for run in sys.argv[1:]:
        n, s, loc, ns = analyse(run.rstrip("/"))
        out.append(f"{run}: {n} annihilations, {ns} samples; r^2 = s (t_a - t): s = {s:.3f} cells^2/t.u. -> mean mobility mu = s/(2 pi K) = {s / (2 * np.pi * K):.4f}"
                   f"  (local s at r 3-5 / 5-8 / 8-12: " + " / ".join(f"{x:.3f}" for x in loc) + ")")
        print(out[-1], flush=True)
    out.append("wall-defect mobility (passive circles, mode relaxation): mu = 0.02; free-defect estimate from gamma1 ~ 50: mu ~ 4/(pi gamma1 ln(r/a)) ~ 0.008-0.01")
    open("pair_mobility.txt", "w").write("\n".join(out) + "\n")
