#!/usr/bin/env python3
"""Precursors of the untilting (user, 2026-10-12): position/orientation correlations between a depleting escaper and the other wall and bulk defects, aligned on departure,
escaper vs resident control (same frame, |x|/x_max within 0.15, > 3 cells away). Shape sweep (defects every 2 t.u.; ellipses and ovals, aspect >= 1.35), second half.
Lags (frames of 2 t.u.): -30 -25 -20 -15 -10 -5 -3 0. Per defect D and lag:
  WALL: arclength gaps to the left/right wall +1/2 neighbours (units of P/n_wall), smaller and larger gap; neighbours' cos tau (mean); tilt-sense agreement sign(tau_D) sign(tau_nb) (mean);
        |tau_D - tau_nb| (deg); layer tilt coherence |<sign tau>|;
  BULK (> 4 cells from the wall): distance to the nearest +1/2 and -1/2; number of bulk defects within 8 cells; for the nearest bulk +1/2: bearing (cos of the angle between D -> it and the
        inward normal at D: > 0 in front), heading towards D (cos between its head and the direction to D), and |orientation difference| of the comet axes (deg).
Usage: precursor.py [NPROC] -> precursor.npy, precursor.txt"""
import os, sys, csv, numpy as np
from multiprocessing import Pool
from scipy.spatial import cKDTree
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
import energetics as En  # noqa
SW = "/home/baloewe/Physics/Overcharged/runs/shapes_server"; LAGS = (-30, -25, -20, -15, -10, -5, -3, 0)
NAMES = ["gap_min", "gap_max", "nb_cos", "nb_same_sense", "nb_dtau", "layer_coh", "d_plus", "d_minus", "n_bulk8", "plus_bearing", "plus_heading", "plus_dori"]
def feats(f, i, tree, Nn, s_w, P, sg):
    d, j = tree.query(f[:, :2]); wall = np.flatnonzero((np.sign(f[:, 2]) == sg) & (d < 4)); bulk = np.flatnonzero(d >= 4)
    tau = {q: np.angle(np.exp(1j * (f[q, 3] - np.angle(Nn[j[q]])))) for q in wall}
    if i not in tau or len(wall) < 3: return None
    sp = np.array([s_w[j[q]] for q in wall]); si = s_w[j[i]]; ds = (sp - si + P / 2) % P - P / 2; others = [q for q in wall if q != i]
    dso = np.array([ds[list(wall).index(q)] for q in others]); left = others[int(np.argmin(np.where(dso < 0, -dso, np.inf)))] if (dso < 0).any() else None
    right = others[int(np.argmin(np.where(dso > 0, dso, np.inf)))] if (dso > 0).any() else None; mean_sp = P / len(wall)
    gaps = [abs(ds[list(wall).index(q)]) / mean_sp for q in (left, right) if q is not None]; nb = [q for q in (left, right) if q is not None]
    out = [min(gaps), max(gaps), np.mean([np.cos(tau[q]) for q in nb]), np.mean([np.sign(tau[i]) * np.sign(tau[q]) for q in nb]),
           np.degrees(np.mean([abs(np.angle(np.exp(1j * (tau[i] - tau[q])))) for q in nb])), abs(np.mean([np.sign(t) for t in tau.values()]))]
    zi = f[i, 0] + 1j * f[i, 1]; nin = -Nn[j[i]]
    bp = [q for q in bulk if f[q, 2] > 0]; bm = [q for q in bulk if f[q, 2] < 0]
    dz = lambda q: (f[q, 0] + 1j * f[q, 1]) - zi
    out += [min(abs(dz(q)) for q in bp) if bp else np.nan, min(abs(dz(q)) for q in bm) if bm else np.nan, sum(1 for q in bulk if abs(dz(q)) < 8)]
    if bp:
        q = min(bp, key=lambda q: abs(dz(q))); u = dz(q) / abs(dz(q)); head = -np.exp(1j * f[q, 3])
        out += [(u * np.conj(nin)).real, (head * np.conj(-u)).real, np.degrees(abs(np.angle(np.exp(1j * (f[q, 3] - f[i, 3])))))]
    else: out += [np.nan] * 3
    return out
def one(rel):
    try:
        run = os.path.join(SW, rel); J, bc, k, t, fr, ev = En.events(run, max_events=10 ** 6); sg = 1 if k > 0 else -1
        W, s_w, P = En.wall_polyline(bc); tree = cKDTree(W); xm = np.abs(W[:, 0]).max(); t = np.array(t); h = len(fr) // 2
        Wc = W[:, 0] + 1j * W[:, 1]; Tn = np.roll(Wc, -1) - np.roll(Wc, 1); Nn = -1j * Tn / np.abs(Tn)
        if np.mean((Nn * np.conj(Wc)).real) < 0: Nn = -Nn
        rng = np.random.default_rng(5); out = []
        for q_ev, e in enumerate(ev):
            n = e["n"]
            if t[n] < t[h]: continue
            ch = e["chain"]; ctrl = None
            for lg in LAGS:
                m = n + lg
                if m not in ch or not len(fr[m]): continue
                f = fr[m]; i = int(np.argmin(np.hypot(*(f[:, :2] - ch[m]).T))); fe = feats(f, i, tree, Nn, s_w, P, sg)
                if fe is not None: out.append(["esc", q_ev, lg] + fe)
                d = tree.query(f[:, :2])[0]
                if ctrl is None or lg == LAGS[0]:
                    res = [q for q in np.flatnonzero((np.sign(f[:, 2]) == sg) & (d < 4)) if np.hypot(*(f[q, :2] - ch[m])) > 3 and abs(abs(f[q, 0]) - abs(ch[m][0])) / xm < 0.15]
                    ctrl = f[rng.choice(res), :2] if res else None
                if ctrl is not None:
                    dd = np.hypot(*(f[:, :2] - ctrl).T); c = int(np.argmin(dd))
                    if dd[c] < 2.5 and d[c] < 4:
                        ctrl = f[c, :2]; fc = feats(f, c, tree, Nn, s_w, P, sg)
                        if fc is not None: out.append(["ctl", q_ev, lg] + fc)
        return rel, out
    except Exception as ex:
        return rel, [["ERR", str(ex)]]
if __name__ == "__main__":
    man = {r["dir"]: r for r in csv.DictReader(open(SW + "/manifest_shape_sweep.tsv"), delimiter="\t")}
    rels = [r for r in man if man[r]["family"] in ("ellipse", "oval") and float(man[r]["aspect"]) >= 1.35 and os.path.exists(os.path.join(SW, r, "DONE"))]
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 3) as pool: res = pool.map(one, rels, chunksize=2)
    rows = [[rel] + r for rel, o in res for r in o if r[0] != "ERR"]; errs = [(rel, o[0][1]) for rel, o in res if o and o[0][0] == "ERR"]
    np.save(os.path.join(H, "precursor.npy"), np.array(rows, dtype=object), allow_pickle=True)
    A = np.array(rows, dtype=object); kind = A[:, 1].astype(str); lag = A[:, 3].astype(int); F = A[:, 4:].astype(float)
    L = [f"{len(rels)} runs (aspect >= 1.35, ellipses + ovals), {len(errs)} errors; {int(((kind == 'esc') & (lag == 0)).sum())} escapes at departure. Escaper vs control: mean ± s.e."]
    for c, nm in enumerate(NAMES):
        row = []
        for lg in LAGS:
            e = F[(kind == "esc") & (lag == lg), c]; q = F[(kind == "ctl") & (lag == lg), c]; e, q = e[np.isfinite(e)], q[np.isfinite(q)]
            diff = e.mean() - q.mean(); se = np.sqrt(e.var() / max(len(e), 1) + q.var() / max(len(q), 1))
            row.append(f"{2*lg:+d}: {e.mean():.3f}/{q.mean():.3f} ({diff / se:+.1f}σ)")
        L.append(f"{nm:14s} " + "  ".join(row))
    L += [f"ERR {r}: {m}" for r, m in errs[:5]]
    print("\n".join(L)); open(os.path.join(H, "precursor.txt"), "w").write("\n".join(L) + "\n")
