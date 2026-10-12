#!/usr/bin/env python3
"""Force model for wall +1/2 (user, 2026-10-11): v_in = alpha*(u_core - u_ring) + beta*u_ring + mu_n*F0 + c, all along the inward wall normal.
Rows: depleting escapers (energetics.events) at lags -20..+3 frames; resident controls (random wall +1/2, same frame and region, NOT escaping) at lags -20, -10, -5, 0; bulk +1/2
(depth > 8; along their head, F0 not computed: alpha_bulk = defect velocity vs own jet). v_in from the chain positions at m-1, m+1 projected on the inward normal at m.
u_core = mean inward flow within 2 cells, u_ring = 4-8 cells; F0 = elastic force on the probe (Frank energy of the measured configuration, 0.5-cell inward displacement; < 0 pulls to the wall).
Usage: forcefit.py TAG RUNDIR [...] -> forcefit_<TAG>.npy (rows) and forcefit_<TAG>.txt (fits, bootstrap over events)"""
import os, sys, numpy as np
from multiprocessing import Pool
from scipy.spatial import cKDTree
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
import energetics as En  # noqa
from selfflow import load_flow  # noqa
En.DELTAS = (0, 0.5)
def one(run):
    J, bc, k, t, fr, ev = En.events(run, max_events=50); Q = bc["Q"]; sg = 1 if k > 0 else -1
    eps = bc["aInv"][0] if "aInv" in bc else 1.0; E = En.bd.Ellipse(float(eps))
    W, s, P = En.wall_polyline(bc); tree = cKDTree(W); Wc = W[:, 0] + 1j * W[:, 1]; Tn = np.roll(Wc, -1) - np.roll(Wc, 1); Nn = -1j * Tn / np.abs(Tn)
    if np.mean((Nn * np.conj(Wc)).real) < 0: Nn = -Nn
    A = np.abs(W[:, 0]).max(); rows, ft, idx = load_flow(run); X = rows[:, 1] + 0.5 - Q[0]; Y = rows[:, 2] + 0.5 - Q[1]; rng = np.random.default_rng(5); out = []
    def flows(tt, p, d):
        ii = idx[ft[np.argmin(np.abs(ft - tt))]]; r = np.hypot(X[ii] - p[0], Y[ii] - p[1]); c = r < 2; g = (r > 4) & (r < 8)
        if not c.any() or not g.any(): return np.nan, np.nan
        vx, vy = rows[ii, 4], rows[ii, 5]; return vx[c].mean() * d.real + vy[c].mean() * d.imag, vx[g].mean() * d.real + vy[g].mean() * d.imag
    def near(m, p, sgn=sg, lim=1.5):
        h = fr[m]; hp = h[np.sign(h[:, 2]) == sgn, :2] if len(h) else np.zeros((0, 2))
        if not len(hp): return None
        dd = np.hypot(*(hp - p).T); kk = np.argmin(dd); return hp[kk] if dd[kk] < lim else None
    def row(kind, q, lag, m, p, pm, pp):
        d, j = tree.query(p); nin = -Nn[j]; uc, ur = flows(t[m], p, nin); L = En.landscape(E, fr[m], p, k)
        if L is None or not np.isfinite(uc): return
        v = ((pp[0] - pm[0]) * nin.real + (pp[1] - pm[1]) * nin.imag) / 4.0          # frames are 2 t.u.: (m+1)-(m-1) = 4 t.u.
        f = fr[m]; i = int(np.argmin(np.hypot(*(f[:, :2] - p).T))); ct = np.cos(np.angle(np.exp(1j * (f[i, 3] - np.angle(Nn[j])))))
        out.append((os.path.basename(os.path.dirname(run)) + "/" + os.path.basename(run), kind, q, lag, d, ct, uc, ur, L["F0"], v))
    for q, e in enumerate(ev):
        n = e["n"]; ch = dict(e["chain"]); c = ch[n]
        for m in range(n + 1, n + 5):
            c2 = near(m, c, lim=2.5)
            if c2 is None: break
            c = c2; ch[m] = c
        for lag in range(-20, 4):
            m = n + lag
            if all(x in ch for x in (m - 1, m, m + 1)): row("esc", q, lag, m, ch[m], ch[m - 1], ch[m + 1])
            if lag in (-20, -10, -5, 0) and m in ch:
                f = fr[m]; dd = tree.query(f[:, :2])[0] if len(f) else np.zeros(0)
                res = [jj for jj in np.flatnonzero((np.sign(f[:, 2]) == sg) & (dd < 4)) if np.hypot(*(f[jj, :2] - ch[m])) > 3 and abs(abs(f[jj, 0]) - abs(ch[m][0])) / A < 0.15]
                if res:
                    pc = f[rng.choice(res), :2]; pm, pp = near(m - 1, pc), near(m + 1, pc)
                    if pm is not None and pp is not None: row("ctl", q, lag, m, pc, pm, pp)
    for m in range(25, len(fr) - 2, 7):                                   # bulk +1/2: own jet ridden?
        f = fr[m]
        if not len(f): continue
        dd = tree.query(f[:, :2])[0]
        for i in np.flatnonzero((f[:, 2] > 0) & (dd > 8))[:3]:
            p = f[i, :2]; pm, pp = near(m - 1, p, 1), near(m + 1, p, 1)
            if pm is None or pp is None: continue
            hd = -np.exp(1j * f[i, 3]); uc, ur = flows(t[m], p, hd)
            if np.isfinite(uc): out.append((os.path.basename(run), "bulk", -1, 0, dd[i], 1.0, uc, ur, 0.0, ((pp[0] - pm[0]) * hd.real + (pp[1] - pm[1]) * hd.imag) / 4.0))
    return out
def fit(R, cols, nb=300):
    y = R["v"].astype(float); Xm = np.c_[[R[c].astype(float) for c in cols]].T if cols else np.zeros((len(y), 0)); Xm = np.c_[Xm, np.ones(len(y))]
    ok = np.all(np.isfinite(Xm), 1) & np.isfinite(y); Xm, y, g = Xm[ok], y[ok], R["ev"][ok]
    b = np.linalg.lstsq(Xm, y, rcond=None)[0]; r2 = 1 - np.sum((y - Xm @ b) ** 2) / np.sum((y - y.mean()) ** 2)
    G = np.unique(g); rng = np.random.default_rng(0); B = []
    for _ in range(nb):
        pick = rng.choice(G, len(G)); ii = np.concatenate([np.flatnonzero(g == z) for z in pick]); B.append(np.linalg.lstsq(Xm[ii], y[ii], rcond=None)[0])
    return b, np.std(B, 0), r2, ok.sum()
if __name__ == "__main__":
    tag = sys.argv[1]; runs = [os.path.abspath(r) for r in sys.argv[2:]]
    with Pool(min(12, len(runs))) as pool: res = pool.map(one, runs, chunksize=1)
    rows = sum(res, []); np.save(os.path.join(H, f"forcefit_{tag}.npy"), np.array(rows, dtype=object), allow_pickle=True)
    a = np.array(rows, dtype=object); R = dict(run=a[:, 0], kind=a[:, 1].astype(str), ev=np.array([f"{x}|{y}" for x, y in zip(a[:, 0], a[:, 2])]), lag=a[:, 3].astype(int),
             d=a[:, 4].astype(float), ct=a[:, 5].astype(float), uc=a[:, 6].astype(float), ur=a[:, 7].astype(float), F0=a[:, 8].astype(float), v=a[:, 9].astype(float))
    R["jet"] = R["uc"] - R["ur"]; L = []
    def sub(mask): return {kk: vv[mask] for kk, vv in R.items()}
    def say(name, S, cols):
        b, e, r2, n = fit(S, cols); L.append(f"{name:52s} n {n:5d}  R2 {r2:.2f}  " + "  ".join(f"{c} {x:+.4f}±{y:.4f}" for c, x, y in zip(cols + ['const'], b, e)))
    for kd, nm in (("esc", "escapers"), ("ctl", "residents"), ("bulk", "bulk +1/2 (along head)")):
        S = sub(R["kind"] == kd); L.append(f"{nm}: n {len(S['v'])}, mean v {np.nanmean(S['v']):+.3f}, jet {np.nanmean(S['jet']):+.3f}, ring {np.nanmean(S['ur']):+.3f}, F0 {np.nanmean(S['F0']):+.1f}")
    L.append("")
    say("bulk +1/2: v = alpha*jet + beta*ring + c", sub(R["kind"] == "bulk"), ["jet", "ur"])
    W_ = (R["kind"] != "bulk")
    say("escapers: v = alpha*jet + beta*ring + mu*F0 + c", sub(R["kind"] == "esc"), ["jet", "ur", "F0"])
    say("escapers, depth < 4 (in the band)", sub((R["kind"] == "esc") & (R["d"] < 4)), ["jet", "ur", "F0"])
    say("escapers, depth >= 4 (leaving)", sub((R["kind"] == "esc") & (R["d"] >= 4)), ["jet", "ur", "F0"])
    say("escapers + residents", sub(W_), ["jet", "ur", "F0"])
    say("residents only", sub(R["kind"] == "ctl"), ["jet", "ur", "F0"])
    R["jetct"] = 0.32 * R["ct"]
    say("escapers + residents, jet from tilt (0.32 cos tau)", sub(W_), ["jetct", "ur", "F0"])
    L.append("\nby lag (escapers): lag(t.u.)  v_meas  model(esc fit)  alpha*jet  beta*ring  mu*F0")
    b, _, _, _ = fit(sub(R["kind"] == "esc"), ["jet", "ur", "F0"])
    for lg in (-20, -15, -10, -5, -2, 0, 1, 2, 3):
        S = sub((R["kind"] == "esc") & (R["lag"] == lg))
        if not len(S["v"]): continue
        m = lambda x: np.nanmean(x); L.append(f"  {2*lg:+4d}   {m(S['v']):+.3f}   {b[0]*m(S['jet'])+b[1]*m(S['ur'])+b[2]*m(S['F0'])+b[3]:+.3f}   {b[0]*m(S['jet']):+.3f}   {b[1]*m(S['ur']):+.3f}   {b[2]*m(S['F0']):+.3f}")
    print("\n".join(L)); open(os.path.join(H, f"forcefit_{tag}.txt"), "w").write("\n".join(L) + "\n")
