#!/usr/bin/env python3
"""Per-defect escape model vs bulk charge (user, 2026-10-12): does the jet mechanism (push = 0.40 cos tau g(d), elastic drift ~const) explain why exits slow as the bulk charges up,
and where would the bulk charge saturate? Second half of each run, every frame (2 t.u.).
Per wall +1/2 (< 4 cells) and frame: depth d, cos tau, push P = 0.40 cos tau g(d) (g from jetmap_ff12), Q_b of the frame, run mean of Q_b, n_wall; EXIT = its chain reaches > 6 cells within
10 frames and stays > 4 for 10 more (as turnover.py). Logistic fit P(exit in 20 t.u.) = 1/(1 + exp(-(a + b P))) pooled over runs; then predicted vs measured hazard by dQ = Q_b - run mean.
Entries per frame (bulk +1/2 crossing into the band) are counted as in turnover.py for the flux balance. Usage: perdef.py TAG RUNDIR [...] (env NPROC) -> perdef_<TAG>.npy/.txt"""
import os, sys, numpy as np
from multiprocessing import Pool
from scipy.spatial import cKDTree
from scipy.optimize import minimize
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
import energetics as En  # noqa
GX = np.array([0.5, 1.5, 2.5, 3.5, 4.5, 6]); GY = np.array([0.30, 0.33, 0.58, 0.86, 0.97, 1.0])
def one(run):
    try:
        J, bc, k, t, fr, ev = En.events(run, max_events=0); eps = float(bc.get("aInv", [1])[0] or 1); lam = J["species"][0].get("tumble", 2.0)
        N0 = 2 * abs(k) - 2; sg = 1 if k > 0 else -1; W, s, P = En.wall_polyline(bc); tree = cKDTree(W); t = np.array(t); h = len(fr) // 2
        Wc = W[:, 0] + 1j * W[:, 1]; Tn = np.roll(Wc, -1) - np.roll(Wc, 1); Nn = -1j * Tn / np.abs(Tn)
        if np.mean((Nn * np.conj(Wc)).real) < 0: Nn = -Nn
        D = [tree.query(f[:, :2]) if len(f) else (np.zeros(0), np.zeros(0, int)) for f in fr]
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
        rows = []; ent = []
        for n in range(h, len(fr) - 22):
            f = fr[n]
            if not len(f): continue
            d, j = D[n]; w = (np.sign(f[:, 2]) == sg) & (d < 4); qb = 0.5 * np.sign(f[d >= 4, 2]).sum(); nw = w.sum()
            for i in np.flatnonzero(w):
                ct = np.cos(np.angle(np.exp(1j * (f[i, 3] - np.angle(Nn[j[i]])))))
                ex = crosses(n, f[i, :2], True); b = np.flatnonzero(d > 4); zi = f[i, 0] + 1j * f[i, 1]
                press = sum(0.5 * np.sign(f[jj, 2]) * ((zi - (f[jj, 0] + 1j * f[jj, 1])) / abs(zi - (f[jj, 0] + 1j * f[jj, 1])) ** 2 * np.conj(Nn[j[i]])).real for jj in b)   # bulk push towards the wall (pi K units / q_i)
                rows.append((n, qb, nw, d[i], ct, 0.40 * ct * np.interp(d[i], GX, GY), float(ex), press))
            ne = sum(1 for i in np.flatnonzero((np.sign(f[:, 2]) == sg) & (d > 6)) if crosses(n, f[i, :2], False))
            ent.append((n, qb, nw, ne))
        return dict(run=run, eps=eps, lam=lam, N0=N0, rows=np.array(rows, float), ent=np.array(ent, float))
    except Exception as ex:
        return dict(run=run, err=str(ex))
if __name__ == "__main__":
    tag = sys.argv[1]; runs = [os.path.abspath(r) for r in sys.argv[2:]]
    with Pool(int(os.environ.get("NPROC", 4))) as pool: res = pool.map(one, runs, chunksize=1)
    ok = [r for r in res if "err" not in r and r["eps"] >= 1.3 and len(r["rows"])]; np.save(os.path.join(H, f"perdef_{tag}.npy"), np.array(res, dtype=object), allow_pickle=True)
    R = np.concatenate([np.c_[r["rows"], np.full(len(r["rows"]), r["rows"][:, 1].mean() if len(r["rows"]) else 0), np.full(len(r["rows"]), i)] for i, r in enumerate(ok)])
    qb, nw, d, ct, Pp, ex, pr, qm = R[:, 1], R[:, 2], R[:, 3], R[:, 4], R[:, 5], R[:, 6], R[:, 7], R[:, 8]; dq = qb - qm
    # an exit is counted once per crossing at its first frame; the per-frame flag marks 'exits within the next 20 t.u.' -> hazard per frame ~ flag rate / 10 frames
    def nll(p): z = p[0] + p[1] * Pp; return np.sum(np.logaddexp(0, z) - ex * z)
    fit = minimize(nll, [-4, 5]).x; pred = 1 / (1 + np.exp(-(fit[0] + fit[1] * Pp)))
    def nll2(p): z = p[0] + p[1] * Pp + p[2] * pr; return np.sum(np.logaddexp(0, z) - ex * z)
    fit2 = minimize(nll2, [-4, 5, 0]).x; pred2 = 1 / (1 + np.exp(-(fit2[0] + fit2[1] * Pp + fit2[2] * pr)))
    def nll3(p): z = p[0] + p[1] * Pp + p[2] * pr + p[3] * d; return np.sum(np.logaddexp(0, z) - ex * z)
    fit3 = minimize(nll3, [-4, 5, 0, 0]).x; pred3 = 1 / (1 + np.exp(-(fit3[0] + fit3[1] * Pp + fit3[2] * pr + fit3[3] * d)))
    L = [f"{len(ok)} elongated runs, {len(R)} wall-defect frames, exit flag rate {ex.mean():.4f}; logistic fit of P(exit within 20 t.u.) on the push P = 0.40 cos tau g(d): a {fit[0]:+.2f}, b {fit[1]:+.2f}",
         f"with the bulk pressure (push of the bulk charge towards the wall): a {fit2[0]:+.2f}, b_jet {fit2[1]:+.2f}, b_press {fit2[2]:+.2f};  + depth: b_depth {fit3[3]:+.2f} (b_jet {fit3[1]:+.2f}, b_press {fit3[2]:+.2f})",
         "dQ_bulk | n | depth | cos tau | push P | bulk pressure | measured exit-flag rate | predicted: push / push+pressure / +depth | ratios normalised to dQ=0: meas, push, +press, +depth"]
    B = np.arange(-1.75, 1.8, 0.5); base_m = base_p = None; tab = []
    for i in range(len(B) - 1):
        m = (dq >= B[i]) & (dq < B[i + 1])
        if m.sum() < 300: continue
        tab.append((0.5 * (B[i] + B[i + 1]), m.sum(), d[m].mean(), ct[m].mean(), Pp[m].mean(), pr[m].mean(), ex[m].mean(), pred[m].mean(), pred2[m].mean(), pred3[m].mean()))
    t0 = min(tab, key=lambda x: abs(x[0]))
    for x in tab: L.append(f"  {x[0]:+5.2f} | {x[1]:6d} | {x[2]:.2f} | {x[3]:+.3f} | {x[4]:+.4f} | {x[5]:+.3f} | {x[6]:.4f} | {x[7]:.4f} {x[8]:.4f} {x[9]:.4f} | {x[6] / t0[6]:.2f} {x[7] / t0[7]:.2f} {x[8] / t0[8]:.2f} {x[9] / t0[9]:.2f}")
    # linear trends of depth and cos tau with dQ -> extrapolated predicted hazard
    cd = np.polyfit(dq, d, 1); cc = np.polyfit(dq, ct, 1)
    cp = np.polyfit(dq, pr, 1); L.append(f"\ntrends per unit bulk charge: depth {cd[0]:+.3f} cells, cos tau {cc[0]:+.3f}, bulk pressure {cp[0]:+.4f}")
    # flux balance per run-frame: exits per frame vs entries per frame by dQ (totals, not per defect)
    E = np.concatenate([np.c_[r["ent"], np.full(len(r["ent"]), r["ent"][:, 1].mean())] for r in ok]); edq = E[:, 1] - E[:, 4]
    L.append("\nflux balance (per frame, totals): dQ | exits per frame (sum of exit flags/10) | entries per frame | net to the bulk")
    for i in range(len(B) - 1):
        m = (dq >= B[i]) & (dq < B[i + 1]); me = (edq >= B[i]) & (edq < B[i + 1])
        if me.sum() < 100: continue
        xf = ex[m].sum() / 10 / me.sum(); ef = E[me, 3].sum() / me.sum()
        L.append(f"  {0.5 * (B[i] + B[i + 1]):+5.2f} | {xf:.4f} | {ef:.4f} | {xf - ef:+.4f}")
    print("\n".join(L)); open(os.path.join(H, f"perdef_{tag}.txt"), "w").write("\n".join(L) + "\n")
