#!/usr/bin/env python3
"""Is the untilting before an escape just rotational diffusion of the tilt? (user, 2026-10-12). Shape sweep (defects every 2 t.u.; ellipses + ovals, aspect >= 1.35), second half.
Variable x = cos tau (tau = tail - outward normal; x = 1: head into the bulk).
(1) Kramers-Moyal from RESIDENT tracks (wall +1/2 followed frame to frame, nearest same-sign within 1.5 cells, while < 4 cells deep; started every 25 frames, up to 30 frames; not
    within 3 cells of an escape chain): drift A(x) = <dx>/dt, diffusion B(x) = <dx^2>/(2 dt) at dt = 2 t.u. (and 6 t.u. for A), binned in x.
    ESCAPERS (energetics.events chains, frames -20..-1 before departure): their <dx>/dt in the same x bins -> excess drift = escaper - resident (with s.e.).
(2) Excursions: in resident tracks, first upward crossing of x >= 0.7 after x <= 0.4; profile x(lag) for 20 frames before the crossing; outcome = leaves the band (> 6 cells within
    10 frames after) or returns. Escapers: same alignment on their last upward crossing of 0.7 before departure. Compare profiles.
(3) P(escape | excursion to x >= 0.7) for residents' excursions; depth at the crossing, escaping vs returning.
Usage: rotdiff.py [NPROC] -> rotdiff.npy, rotdiff.txt"""
import os, sys, csv, numpy as np
from multiprocessing import Pool
from scipy.spatial import cKDTree
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
import energetics as En  # noqa
SW = "/home/baloewe/Physics/Overcharged/runs/shapes_server"; XB = np.linspace(-1, 1, 11)
def one(rel):
    try:
        run = os.path.join(SW, rel); J, bc, k, t, fr, ev = En.events(run, max_events=10 ** 6); sg = 1 if k > 0 else -1
        W, s_w, P = En.wall_polyline(bc); tree = cKDTree(W); t = np.array(t); h = len(fr) // 2
        Wc = W[:, 0] + 1j * W[:, 1]; Tn = np.roll(Wc, -1) - np.roll(Wc, 1); Nn = -1j * Tn / np.abs(Tn)
        if np.mean((Nn * np.conj(Wc)).real) < 0: Nn = -Nn
        def xd(m, p):
            f = fr[m]; i = int(np.argmin(np.hypot(*(f[:, :2] - p).T))); d, j = tree.query(f[i, :2]); return np.cos(np.angle(np.exp(1j * (f[i, 3] - np.angle(Nn[j]))))), d
        def nxt(m, c, lim=1.5):
            g = fr[m]; gp = g[np.sign(g[:, 2]) == sg, :2] if len(g) else np.zeros((0, 2))
            if not len(gp): return None
            dd = np.hypot(*(gp - c).T); kk = np.argmin(dd); return gp[kk] if dd[kk] < lim else None
        # escaper chains (second half)
        esc_pts = {}; esc_series = []
        for e in ev:
            n = e["n"]
            if t[n] < t[h]: continue
            ch = e["chain"]; ser = [(m - n, *xd(m, ch[m])) for m in sorted(ch) if m <= n]
            esc_series.append(ser)
            for m, p in ch.items(): esc_pts.setdefault(m, []).append(p)
        # resident tracks
        tracks = []
        for n0 in range(h, len(fr) - 32, 25):
            f = fr[n0]
            if not len(f): continue
            d, j = tree.query(f[:, :2])
            for i in np.flatnonzero((np.sign(f[:, 2]) == sg) & (d < 4)):
                p = f[i, :2]
                if any(np.hypot(*(np.array(esc_pts.get(m, [])) - p).T).min() < 3 for m in range(n0, n0 + 30) if esc_pts.get(m)): continue
                ser = []; c = p
                for m in range(n0, n0 + 30):
                    if m > n0: c = nxt(m, c)
                    if c is None: break
                    x, dep = xd(m, c); ser.append((m - n0, x, dep))
                    if dep > 6: break
                tracks.append(ser)
        return dict(rel=rel, esc=esc_series, res=tracks)
    except Exception as ex:
        return dict(rel=rel, err=str(ex))
def km(pairs):
    x0 = np.array([p[0] for p in pairs]); dx = np.array([p[1] for p in pairs]); dt = np.array([p[2] for p in pairs]); b = np.digitize(x0, XB) - 1; out = []
    for i in range(len(XB) - 1):
        m = b == i
        if m.sum() < 30: out.append((np.nan, np.nan, np.nan, 0)); continue
        v = dx[m] / dt[m]; out.append((v.mean(), v.std() / np.sqrt(m.sum()), np.mean(dx[m] ** 2 / (2 * dt[m])), m.sum()))
    return out
if __name__ == "__main__":
    man = {r["dir"]: r for r in csv.DictReader(open(SW + "/manifest_shape_sweep.tsv"), delimiter="\t")}
    rels = [r for r in man if man[r]["family"] in ("ellipse", "oval") and float(man[r]["aspect"]) >= 1.35 and os.path.exists(os.path.join(SW, r, "DONE"))]
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 3) as pool: res = pool.map(one, rels, chunksize=2)
    ok = [r for r in res if "err" not in r]; np.save(os.path.join(H, "rotdiff.npy"), np.array(ok, dtype=object), allow_pickle=True)
    # (1) Kramers-Moyal
    rp2, rp6, ep2 = [], [], []
    for r in ok:
        for ser in r["res"]:
            inb = [s for s in ser if s[2] < 4]
            for a, b in zip(inb[:-1], inb[1:]):
                if b[0] - a[0] == 1: rp2.append((a[1], b[1] - a[1], 2.0))
            for a, b in zip(inb[:-3], inb[3:]):
                if b[0] - a[0] == 3: rp6.append((a[1], b[1] - a[1], 6.0))
        for ser in r["esc"]:
            s = [q for q in ser if -20 <= q[0] <= 0]
            for a, b in zip(s[:-1], s[1:]):
                if b[0] - a[0] == 1 and a[0] < 0: ep2.append((a[1], b[1] - a[1], 2.0))
    R2, R6, E2 = km(rp2), km(rp6), km(ep2)
    L = [f"{len(ok)} runs ({len(res) - len(ok)} errors); resident step pairs {len(rp2)}, escaper step pairs (-40..0 t.u.) {len(ep2)}",
         "(1) x = cos tau bin | resident drift A (per t.u., dt 2) ± s.e. | resident drift (dt 6) | resident diffusion B (dt 2) | escaper drift (dt 2) ± s.e. | excess (sigma) | n_res n_esc"]
    for i in range(len(XB) - 1):
        a, b6, e = R2[i], R6[i], E2[i]; ex = (e[0] - a[0]) / np.sqrt(a[1] ** 2 + e[1] ** 2) if np.isfinite(e[0]) and np.isfinite(a[0]) else np.nan
        L.append(f"  [{XB[i]:+.1f},{XB[i+1]:+.1f}) | {a[0]:+.4f} ± {a[1]:.4f} | {b6[0]:+.4f} | {a[2]:.4f} | {e[0]:+.4f} ± {e[1]:.4f} | {ex:+.1f} | {a[3]} {e[3]}")
    # (2)-(3) excursions
    prof_ret, prof_esc_res, depth_ret, depth_esc = [], [], [], []
    for r in ok:
        for ser in r["res"]:
            xs = np.array([q[1] for q in ser]); ds = np.array([q[2] for q in ser]); lo = np.flatnonzero(xs <= 0.4)
            if not len(lo): continue
            up = [q for q in range(lo[0] + 1, len(xs)) if xs[q] >= 0.7 and ds[q] < 4]
            if not up or up[0] < 8: continue
            c = up[0]; prof = np.full(21, np.nan); seg = xs[max(0, c - 20):c + 1]; prof[21 - len(seg):] = seg
            left = (ds[c:c + 11] > 6).any()
            (prof_esc_res if left else prof_ret).append(prof); (depth_esc if left else depth_ret).append(ds[c])
    prof_esc = []
    for r in ok:
        for ser in r["esc"]:
            lg = np.array([q[0] for q in ser]); xs = np.array([q[1] for q in ser]); ups = [q for q in range(1, len(xs)) if xs[q] >= 0.7 > xs[q - 1]]
            if not ups: continue
            c = ups[-1]; prof = np.full(21, np.nan); seg = xs[max(0, c - 20):c + 1]; prof[21 - len(seg):] = seg; prof_esc.append(prof)
    L.append(f"\n(2) excursion profiles aligned on the upward crossing of cos tau = 0.7 (t.u. before crossing): escapers (n {len(prof_esc)}) | resident excursions that return (n {len(prof_ret)}) | resident excursions that escape (n {len(prof_esc_res)})")
    for q in (0, 5, 10, 13, 15, 17, 19, 20):
        f_ = lambda P_: np.nanmean(np.array(P_)[:, q]) if len(P_) else np.nan
        L.append(f"  {2 * (q - 20):+d}: {f_(prof_esc):+.3f} | {f_(prof_ret):+.3f} | {f_(prof_esc_res):+.3f}")
    ntot = len(prof_ret) + len(prof_esc_res)
    L.append(f"\n(3) resident excursions to cos tau >= 0.7: {ntot}; escape within 20 t.u. after the crossing: {len(prof_esc_res)} ({len(prof_esc_res) / max(ntot, 1):.1%}); depth at crossing: escaping {np.mean(depth_esc):.2f} vs returning {np.mean(depth_ret):.2f}")
    print("\n".join(L)); open(os.path.join(H, "rotdiff.txt"), "w").write("\n".join(L) + "\n")
