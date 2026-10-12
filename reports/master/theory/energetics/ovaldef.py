#!/usr/bin/env python3
"""Ovals vs ellipses, defect-only escape analyses (user, 2026-10-12). Shape sweep (defects every 2 t.u.), same code for both families; region by |x|/x_max of the wall site
(tip > 0.866, flank < 0.5, mid otherwise; for ellipses this equals the angle-parameter rule). Per run, second half:
 (a) depleting escapes (energetics.events): count by region; escaper cos tau at -40, -20, -10, 0 t.u.; resident controls (same frame, |x|/x_max within 0.15) at the same lags;
 (b) resident wall +1/2: occupancy and cos tau by region (every 5th frame) -> escape rate per wall defect per region;
 (c) crowding: exits (turnover definition) per wall defect vs n_wall - (2k-2), sampled every 5th frame.
Usage: ovaldef.py [NPROC] -> ovaldef.npy, ovaldef.txt"""
import os, sys, csv, numpy as np
from multiprocessing import Pool
from scipy.spatial import cKDTree
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
import energetics as En  # noqa
SW = "/home/baloewe/Physics/Overcharged/runs/shapes_server"; LAGS = (-20, -10, -5, 0)
def reg(x, xm): u = abs(x) / xm; return 0 if u > 0.866 else (2 if u < 0.5 else 1)
def one(rel):
    try:
        run = os.path.join(SW, rel); J, bc, k, t, fr, ev = En.events(run, max_events=10 ** 6); sg = 1 if k > 0 else -1; N0 = 2 * abs(k) - 2
        W, s, P = En.wall_polyline(bc); tree = cKDTree(W); xm = np.abs(W[:, 0]).max(); t = np.array(t); h = len(fr) // 2
        Wc = W[:, 0] + 1j * W[:, 1]; Tn = np.roll(Wc, -1) - np.roll(Wc, 1); Nn = -1j * Tn / np.abs(Tn)
        if np.mean((Nn * np.conj(Wc)).real) < 0: Nn = -Nn
        def ct_of(f, i, j): return np.cos(np.angle(np.exp(1j * (f[i, 3] - np.angle(Nn[j])))))
        rng = np.random.default_rng(3); esc = []
        for e in ev:
            n = e["n"]
            if t[n] < t[h]: continue
            ch = e["chain"]; g = reg(ch[n][0], xm); row = [g]
            for lg in LAGS:
                m = n + lg
                if m not in ch: row += [np.nan, np.nan]; continue
                f = fr[m]; d, j = tree.query(f[:, :2]); i = int(np.argmin(np.hypot(*(f[:, :2] - ch[m]).T))); row.append(ct_of(f, i, j[i]))
                res = [q for q in np.flatnonzero((np.sign(f[:, 2]) == sg) & (d < 4)) if np.hypot(*(f[q, :2] - ch[m])) > 3 and abs(abs(f[q, 0]) - abs(ch[m][0])) / xm < 0.15]
                row.append(ct_of(f, (q := rng.choice(res)), j[q]) if res else np.nan)
            esc.append(row)
        occ = np.zeros(3); cs = [[], [], []]; nfr = 0
        def nxt(m, c, lim=2.5):
            g_ = fr[m]; gp = g_[np.sign(g_[:, 2]) == sg, :2] if len(g_) else np.zeros((0, 2))
            if not len(gp): return None
            dd = np.hypot(*(gp - c).T); kk = np.argmin(dd); return gp[kk] if dd[kk] < lim else None
        def exits(n, c):
            for m in range(n + 1, min(n + 11, len(fr))):
                c = nxt(m, c)
                if c is None: return False
                if tree.query(c)[0] > 6:
                    for mm in range(m + 1, min(m + 11, len(fr))):
                        c = nxt(mm, c)
                        if c is None or tree.query(c)[0] < 4: return False
                    return True
            return False
        crowd = {}
        for n in range(h, len(fr) - 22, 5):
            f = fr[n]
            if not len(f): continue
            d, j = tree.query(f[:, :2]); w = np.flatnonzero((np.sign(f[:, 2]) == sg) & (d < 4)); nfr += 1; c = len(w) - N0
            for i in w: g = reg(f[i, 0], xm); occ[g] += 1; cs[g].append(ct_of(f, i, j[i]))
            ex = sum(exits(n, f[i, :2]) for i in w); a_ = crowd.setdefault(c, [0, 0]); a_[0] += len(w); a_[1] += ex
        T = t[-1] - t[h]
        return dict(rel=rel, k=k, esc=np.array(esc).reshape(-1, 1 + 2 * len(LAGS)), occ=occ / max(nfr, 1), cos=[np.mean(x) if x else np.nan for x in cs], T=T, crowd=crowd)
    except Exception as ex:
        return dict(rel=rel, err=str(ex))
if __name__ == "__main__":
    man = {r["dir"]: r for r in csv.DictReader(open(SW + "/manifest_shape_sweep.tsv"), delimiter="\t")}
    rels = [r for r in man if man[r]["family"] in ("ellipse", "oval") and os.path.exists(os.path.join(SW, r, "DONE"))]
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 6) as pool: res = pool.map(one, rels, chunksize=2)
    np.save(os.path.join(H, "ovaldef.npy"), np.array(res, dtype=object), allow_pickle=True)
    ok = [r for r in res if "err" not in r]; L = [f"{len(ok)} runs ({len(res) - len(ok)} errors)"]; RN = ("tip", "mid", "flank")
    for fam in ("ellipse", "oval"):
        for asp in sorted(set(round(float(man[r["rel"]]["aspect"]), 1) for r in ok if man[r["rel"]]["family"] == fam)):
            S = [r for r in ok if man[r["rel"]]["family"] == fam and round(float(man[r["rel"]]["aspect"]), 1) == asp]
            E = np.concatenate([r["esc"] for r in S]) if any(len(r["esc"]) for r in S) else np.zeros((0, 9))
            occT = sum(r["occ"] * r["T"] for r in S); nE = np.array([(E[:, 0] == g).sum() for g in range(3)]); rate = 1000 * nE / np.maximum(occT, 1e-9)
            cosr = np.nanmean([r["cos"] for r in S], 0)
            pre = " ".join(f"{np.nanmean(E[:, 1 + 2 * q]):+.2f}/{np.nanmean(E[:, 2 + 2 * q]):+.2f}" for q in range(len(LAGS))) if len(E) else "-"
            cr = {}
            for r in S:
                for c, (nw, ne) in r["crowd"].items(): a_ = cr.setdefault(c, [0, 0]); a_[0] += nw; a_[1] += ne
            crs = " ".join(f"{c:+d}:{1000 * cr[c][1] / cr[c][0]:.0f}" for c in sorted(cr) if cr[c][0] > 2000 and -5 <= c <= 2)
            L.append(f"{fam:7s} aspect {asp:.1f} ({len(S)} runs): depleting escapes tip/mid/flank {nE.tolist()}; rate per wall defect per 1000 t.u. {np.round(rate, 3).tolist()}; resident cos tau {np.round(cosr, 2).tolist()}")
            L.append(f"        escaper/control cos tau at -40/-20/-10/0 t.u.: {pre}")
            L.append(f"        exits per 1000 wall-defect samples vs (n_wall - (2k-2)): {crs}")
    print("\n".join(L)); open(os.path.join(H, "ovaldef.txt"), "w").write("\n".join(L) + "\n")
