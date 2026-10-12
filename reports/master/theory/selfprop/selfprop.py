#!/usr/bin/env python3
"""Do wall +1/2 in elongated cavities sit deeper and point their heads more into the bulk (stronger inward self-propulsion)? (B. Loewe, 2026-10-09)
Shape sweep (675 runs), second half, frames every 10 t.u. For every wall +1/2 (depth < 4): depth d (distance to the wall polyline), tilt tau = tail angle - outward normal
(defects.dat angle = tail; head = -tail, so the inward component of the self-propulsion is cos tau), region (|x|/A: tip > 1/2, flank < 1/4), and escape within 20 t.u.
(nearest-neighbour chain < 6 cells per 10 t.u., ends deeper than 6). Also the radial drift of each wall defect over 10 t.u. (inward > 0) vs cos tau.
Writes selfprop_events.npz (pooled) and selfprop.txt."""
import csv, glob, json, os, sys
import numpy as np
from multiprocessing import Pool
from scipy.spatial import cKDTree
ROOT = "/home/baloewe/Physics/Overcharged/runs"; sys.path.insert(0, ROOT + "/cavity")
from edge_criteria import read_with_angle  # noqa
from wall_fourier import wall_polyline  # noqa
H = os.path.dirname(os.path.abspath(__file__)); SW = ROOT + "/shapes_server"
def one(rel):
    run = os.path.join(SW, rel)
    try:
        bc = json.load(open(run + "/input.json"))["BC"][0]; Q = bc["Q"]
        W, s, P = wall_polyline(bc); tree = cKDTree(W); A = np.abs(W[:, 0]).max(); Wc = W[:, 0] + 1j * W[:, 1]
        Tn = np.roll(Wc, -1) - np.roll(Wc, 1); Nn = -1j * Tn / np.abs(Tn)
        if np.mean((Nn * np.conj(Wc)).real) < 0: Nn = -Nn                                   # make the normal point outward
        t, fr = read_with_angle(run, np.array(Q[:2], float), 0.5); fr = fr[::5]
        PL = [f[f[:, 2] > 0] if len(f) else np.zeros((0, 4)) for f in fr]; rows = []
        for a in range(len(PL) - 2):
            p = PL[a]
            if not len(p): continue
            d, j = tree.query(p[:, :2]); w = d < 4
            if not w.any(): continue
            pos = p[w, :2]; dd = d[w]; jj = j[w]; tau = np.angle(np.exp(1j * (p[w, 3] - np.angle(Nn[jj])))); reg = np.where(np.abs(pos[:, 0]) / A > 0.5, 0, np.where(np.abs(pos[:, 0]) / A < 0.25, 2, 1))
            alive = np.ones(len(pos), bool); cur = pos.copy(); drift = np.full(len(pos), np.nan); dfin = np.zeros(len(pos))
            for b in (1, 2):
                if not len(PL[a + b]): alive[:] = False; break
                q, k = cKDTree(PL[a + b][:, :2]).query(cur); alive &= q < 6; cur = PL[a + b][k, :2]; dn, _ = tree.query(cur)
                if b == 1: drift = np.where(q < 6, dn - dd, np.nan)                           # depth change in 10 t.u. (> 0 = inward)
                dfin = dn
            esc = alive & (dfin > 6)
            for i in range(len(pos)): rows.append((dd[i], np.cos(tau[i]), reg[i], esc[i], drift[i]))
        return rel, np.array(rows)
    except Exception as ex:
        return rel, None
if __name__ == "__main__":
    man = {r["dir"]: r for r in csv.DictReader(open(SW + "/manifest_shape_sweep.tsv"), delimiter="\t")}
    runs = sorted(os.path.relpath(os.path.dirname(f), SW) for f in glob.glob(SW + "/*/*/*/*/s*/DONE"))
    with Pool(8) as pool: res = pool.map(one, runs, chunksize=4)
    D, C, RG, E, DR, EL, RQ = [], [], [], [], [], [], []
    for rel, a in res:
        if a is None or not len(a): continue
        m = man[rel]; D.append(a[:, 0]); C.append(a[:, 1]); RG.append(a[:, 2]); E.append(a[:, 3]); DR.append(a[:, 4])
        EL.append(np.full(len(a), float(m["aspect"]) >= 1.25)); RQ.append(np.full(len(a), float(m["Req_la"])))
    D, C, RG, E, DR, EL, RQ = map(np.concatenate, (D, C, RG, E, DR, EL, RQ)); EL = EL.astype(bool); E = E.astype(bool)
    np.savez_compressed(os.path.join(H, "selfprop_events.npz"), depth=D, cos_tau=C, region=RG, esc=E, drift=DR, elong=EL, Req=RQ)
    L = []; say = lambda x: (L.append(x), print(x, flush=True))
    say(f"{len(D)} wall-defect samples")
    say("\n== by class and R_eq: mean depth, mean cos(tau) (inward head component), P(escape in 20 t.u.) ==")
    for el in (False, True):
        for R in (1.0, 1.5, 2.0):
            m = (EL == el) & (RQ == R)
            say(f"  {'elongated' if el else 'circle-like':11s} R_eq {R}: n {m.sum():7d}  depth {D[m].mean():.2f}  cos tau {C[m].mean():+.3f}  P(esc) {E[m].mean():.4f}  "
                f"inward drift {np.nanmean(DR[m]):+.3f}")
    say("\n== elongated, by region (tip / middle / flank) ==")
    for g in (0, 1, 2):
        m = EL & (RG == g); say(f"  region {['tip', 'middle', 'flank'][g]:6s}: depth {D[m].mean():.2f}  cos tau {C[m].mean():+.3f}  P(esc) {E[m].mean():.4f}")
    say("\n== P(escape) binned by depth, within class ==")
    for el in (False, True):
        line = f"  {'elongated' if el else 'circle-like':11s}: "
        for lo, hi in ((0, 1), (1, 1.5), (1.5, 2), (2, 2.5), (2.5, 3), (3, 4)):
            m = (EL == el) & (D >= lo) & (D < hi); line += f"{lo}-{hi}: {E[m].mean():.4f} ({m.mean() if False else m.sum()})  "
        say(line)
    say("\n== P(escape) binned by cos(tau), within class ==")
    for el in (False, True):
        line = f"  {'elongated' if el else 'circle-like':11s}: "
        for lo, hi in ((-1, 0), (0, 0.25), (0.25, 0.5), (0.5, 0.75), (0.75, 0.9), (0.9, 1.01)):
            m = (EL == el) & (C >= lo) & (C < hi); line += f"[{lo},{hi}): {E[m].mean():.4f} (n {m.sum()})  "
        say(line)
    say("\n== inward drift in 10 t.u. vs cos(tau), within class (self-propulsion test) ==")
    for el in (False, True):
        line = f"  {'elongated' if el else 'circle-like':11s}: "
        for lo, hi in ((-1, 0), (0, 0.5), (0.5, 0.9), (0.9, 1.01)):
            m = (EL == el) & (C >= lo) & (C < hi); line += f"[{lo},{hi}): {np.nanmean(DR[m]):+.3f}  "
        say(line)
    say("\n== fraction of the escape-rate ratio (elongated/circle-like) explained by the depth and cos(tau) distributions ==")
    # reweight: circle-like conditional P(esc | depth bin, cos bin) applied to the elongated distribution
    db = np.digitize(D, [1, 1.5, 2, 2.5, 3]); cb = np.digitize(C, [0, 0.25, 0.5, 0.75, 0.9]); key = db * 10 + cb
    pc = {k: E[(~EL) & (key == k)].mean() for k in np.unique(key[~EL]) if ((~EL) & (key == k)).sum() > 200}
    pred = np.mean([pc.get(k, np.nan) for k in key[EL]]); obs_e = E[EL].mean(); obs_c = E[~EL].mean()
    say(f"  circle-like P(esc) {obs_c:.4f}; elongated observed {obs_e:.4f}; elongated predicted from its depth x cos(tau) distribution with circle-like rates {pred:.4f}"
        f"  -> explains {(pred - obs_c) / (obs_e - obs_c):.0%} of the excess")
    open(os.path.join(H, "selfprop.txt"), "w").write("\n".join(L) + "\n")
