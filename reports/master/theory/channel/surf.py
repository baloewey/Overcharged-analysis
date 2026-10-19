#!/usr/bin/env python3
"""'Surfing' +1/2 in channels (user, 2026-10-19): fast fluid pockets that carry +1/2 defects. Are they the defects' own active jets or an advecting stream?
Per run (window T0..T0+TL, flow and defects every 2 t.u.): (1) fast cells = top 2% of |v| inside the channel (excluding 3 cells at each wall): fraction within 3 cells
of a bulk +1/2, of a bulk -1/2, and expected by chance (area fraction of those disks); (2) for bulk +1/2 (>= 4 cells from the walls): local flow u (mean over 2 cells),
head direction h = tail angle + pi; <cos(u, h)>, |u| vs the channel rms; component of u along the channel mean flow; (3) defect velocity w over 2 t.u. (nearest
same-sign link < 2.5 cells): <cos(w, h)>, <cos(w, u)>, |w|/|u|; same for -1/2 (no head; w vs u). -> surf.txt"""
import os, sys, json, gzip, numpy as np
from scipy.spatial import cKDTree
sys.path.insert(0, "/home/baloewe/Physics/Overcharged/runs/cavity"); from edge_criteria import read_with_angle  # noqa
R_ = "/home/baloewe/Physics/Overcharged/runs"; H = os.path.dirname(os.path.abspath(__file__))
def analyse(r, T0=3000.0, TL=1000.0):
    run = f"{R_}/{r}"; J = json.load(open(run + "/input.json")); L, Hd = J["domain"][:2]; y0, y1 = J["BC"][2]["Q"][1], J["BC"][3]["Q"][1]
    F = {}
    with gzip.open(run + "/flowfield.dat.gz", "rt") as fh:
        for l in fh:
            p = l.split()
            if len(p) != 7 or not p[0][0].isdigit(): continue
            tt = float(p[0])
            if tt < T0: continue
            if tt > T0 + TL: break
            F.setdefault(tt, []).append([float(v) for v in p[1:6]])
    t, fr = read_with_angle(run, np.array([0.0, 0.0]), 0.0); t = np.array(t); fi = {round(x, 3): i for i, x in enumerate(t)}
    fast_p = fast_m = fast_chance_p = fast_chance_m = nfast = 0; P, Mn = [], []
    meanflow = np.mean([np.mean(np.array(v)[:, 3]) for v in F.values()]); sgn = np.sign(meanflow) or 1
    for tt in sorted(F):
        if round(tt, 3) not in fi or round(tt + 2, 3) not in fi: continue
        a = np.array(F[tt]); ix, iy = a[:, 0].astype(int), a[:, 1].astype(int); U = np.zeros((L, Hd)); V = np.zeros((L, Hd)); U[ix, iy] = a[:, 3]; V[ix, iy] = a[:, 4]
        f = fr[fi[round(tt, 3)]]; g = fr[fi[round(tt + 2, 3)]]
        if not len(f): continue
        bulk = (f[:, 1] - y0 >= 4) & (y1 - f[:, 1] >= 4); pp, mm = f[bulk & (f[:, 2] > 0)], f[bulk & (f[:, 2] < 0)]
        sp = np.hypot(U, V); yc = np.arange(Hd) + 0.5; inner = (yc > y0 + 3) & (yc < y1 - 3); S = sp[:, inner]
        thr = np.percentile(S, 98); fx, fy = np.nonzero(sp * inner[None, :] >= thr); cx, cy = fx + 0.5, fy + 0.5; nfast += len(cx)
        def near(D):
            if not len(D): return np.zeros(len(cx), bool)
            dx = (cx[:, None] - D[None, :, 0] + L / 2) % L - L / 2; return (np.hypot(dx, cy[:, None] - D[None, :, 1]) < 3).any(1)
        fast_p += near(pp).sum(); fast_m += near(mm).sum()
        area = L * inner.sum(); fast_chance_p += len(cx) * min(1, len(pp) * np.pi * 9 / area); fast_chance_m += len(cx) * min(1, len(mm) * np.pi * 9 / area)
        def local(x, y):
            m = (np.abs(((np.arange(L)[:, None] + 0.5 - x + L / 2) % L - L / 2)) < 2) & (np.abs(np.arange(Hd)[None, :] + 0.5 - y) < 2)
            return U[m].mean(), V[m].mean()
        gs = {s: g[np.sign(g[:, 2]) == s] for s in (1, -1)} if len(g) else {1: np.zeros((0, 4)), -1: np.zeros((0, 4))}
        for D, s, out in ((pp, 1, P), (mm, -1, Mn)):
            for x, y, q, ang in D:
                u = np.array(local(x, y)); G = gs[s]
                if not len(G): continue
                dx = (G[:, 0] - x + L / 2) % L - L / 2; dd = np.hypot(dx, G[:, 1] - y); j = np.argmin(dd)
                if dd[j] > 2.5: continue
                w = np.array([dx[j], G[j, 1] - y]) / 2.0; h = np.array([np.cos(ang + np.pi), np.sin(ang + np.pi)]) if s > 0 else None
                out.append((u[0], u[1], w[0], w[1], *(h if h is not None else (np.nan, np.nan))))
    P, Mn = np.array(P), np.array(Mn); rms = np.sqrt(np.mean([np.mean(np.array(v)[:, 3] ** 2 + np.array(v)[:, 4] ** 2) for v in F.values()]))
    cos = lambda a, b: np.sum(a * b, 1) / np.maximum(np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1), 1e-12)
    uP, wP, hP = P[:, 0:2], P[:, 2:4], P[:, 4:6]; uM, wM = Mn[:, 0:2], Mn[:, 2:4]
    return [f"{r}  (window {T0:.0f}-{T0+TL:.0f} t.u.; mean channel flow {meanflow:+.3f}, rms speed {rms:.3f})",
            f"  fast cells (top 2% speed, inner channel) within 3 cells of a bulk +1/2: {fast_p/nfast:.2f} (chance {fast_chance_p/nfast:.2f}); of a bulk -1/2: {fast_m/nfast:.2f} (chance {fast_chance_m/nfast:.2f})",
            f"  bulk +1/2 (n {len(P)}): local |u| {np.linalg.norm(uP,axis=1).mean():.3f} ({np.linalg.norm(uP,axis=1).mean()/rms:.2f} x rms); <cos(u, head)> {cos(uP,hP).mean():+.2f}; u along mean flow {np.mean(uP[:,0]*sgn):+.3f}",
            f"     defect velocity |w| {np.linalg.norm(wP,axis=1).mean():.3f}; <cos(w, head)> {cos(wP,hP).mean():+.2f}; <cos(w, u)> {cos(wP,uP).mean():+.2f}; w along mean flow {np.mean(wP[:,0]*sgn):+.3f}",
            f"  bulk -1/2 (n {len(Mn)}): local |u| {np.linalg.norm(uM,axis=1).mean():.3f} ({np.linalg.norm(uM,axis=1).mean()/rms:.2f} x rms); |w| {np.linalg.norm(wM,axis=1).mean():.3f}; <cos(w, u)> {cos(wM,uM).mean():+.2f}; w along mean flow {np.mean(wM[:,0]*sgn):+.3f}"]
if __name__ == "__main__":
    from multiprocessing import Pool
    runs = ["channel_current/W24_sym_s5", "channel_current/W24_mixed_s3", "channel_pressure/sym_g+0.001_s2", "channel_pressure/mixed_g+0.003_s1"]
    with Pool(4) as p: res = p.map(analyse, runs)
    L = ["SURFING +1/2? fast pockets vs defects (surf.py)"] + [l for r in res for l in r]
    open(os.path.join(H, "surf.txt"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
