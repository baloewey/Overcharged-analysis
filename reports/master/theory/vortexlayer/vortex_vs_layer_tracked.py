#!/usr/bin/env python3
"""Is the second vortex the driver of the wall-layer instability in ellipses? (user's hypothesis, 2026-10-07)
Active ellipses W3 (eps 1.2, 1.5, 2.0; s01, s02; flow averaged over 100 t.u.). Per flow interval: vortex state (halves x<0 / x>0 co- or
counter-rotating, as before). Per defect frame (dt = 2): wall +1/2 (depth < 4), escapes = a wall +1/2 found deeper than 6 cells in the next
frame (nearest-neighbour match < 3 cells), entries = the reverse; wall count deficit (vacancies relative to 2k-2 = 10). Positions classified
by the nearest wall point: 'meeting' (|x| < A/4: where the two vortices meet the wall in a left/right two-vortex state) vs 'tips' (|x| > A/2).
Reported: escape rate per wall defect per 1000 t.u. and mean vacancies, in one- vs two-vortex intervals, total and by location; lagged
correlation between the two-vortex indicator and the escape count per interval (does the second vortex come first?). Writes vortex_vs_layer_tracked.txt."""
import gzip, json, os, sys
import numpy as np
from scipy.spatial import cKDTree
ROOT = "/home/baloewe/Physics/Overcharged/runs"; sys.path.insert(0, ROOT + "/cavity")
from edge_criteria import read_with_angle  # noqa
from wall_fourier import wall_polyline, project  # noqa
H = os.path.dirname(os.path.abspath(__file__)); L = []
def say(s): print(s, flush=True); L.append(s)
for eps in ("1.2", "1.5", "2.0"):
    agg = {0: dict(esc=0, esc_meet=0, esc_tip=0, wt=0.0, wt_meet=0.0, wt_tip=0.0, vac=[]), 1: dict(esc=0, esc_meet=0, esc_tip=0, wt=0.0, wt_meet=0.0, wt_tip=0.0, vac=[])}
    lagged = []
    for sd in ("s01", "s02"):
        run = f"{ROOT}/cassini/W3_ellipse/eps{eps}_k06/{sd}"; bc = json.load(open(run + "/input.json"))["BC"][0]; Q = bc["Q"]
        W, s, P = wall_polyline(bc); A = np.abs(W[:, 0]).max(); B = np.abs(W[:, 1]).max()
        rows = np.loadtxt((l for l in gzip.open(run + "/flowfield.dat.gz", "rt") if len(l.split()) == 7 and l.split()[0][0].isdigit()))
        tf = np.unique(rows[:, 0]); tf = tf[tf > tf.max() * 0.2]; state = {}
        for tt in tf:
            r_ = rows[rows[:, 0] == tt]; x = r_[:, 1] + 0.5 - Q[0]; y = r_[:, 2] + 0.5 - Q[1]; m = (x / A) ** 2 + (y / B) ** 2 < 0.9
            def hr(sel):
                cx, cy = x[sel].mean(), y[sel].mean(); dx, dy = x[sel] - cx, y[sel] - cy
                return np.mean((dx * r_[sel, 5] - dy * r_[sel, 4]) / np.maximum(np.hypot(dx, dy), 0.5))
            state[tt] = int(hr(m & (x < 0)) * hr(m & (x > 0)) < 0)
        t, fr = read_with_angle(run, np.array(Q[:2], float), 0.0); t = np.array(t)
        per_int = {tt: [0, 0.0] for tt in tf}
        # precompute per-frame plus defects and depths
        PL, DEP = [], []
        for f in fr:
            pl = f[f[:, 2] > 0]; dep, _ = project(W, s, pl[:, :2]) if len(pl) else (np.array([]), None); PL.append(pl); DEP.append(dep)
        NF = 10                                                                       # follow for 10 frames = 20 t.u.
        for i in range(len(fr) - NF):
            if t[i] <= tf[0] - 100 or t[i] > tf[-1]: continue
            j_ = np.searchsorted(tf, t[i]); iv = tf[j_] if j_ < len(tf) else None
            if iv is None: continue
            st = state[iv]; w = DEP[i] < 4
            agg[st]["vac"].append(10 - w.sum())
            if not w.any(): continue
            pos = PL[i][w, :2].copy(); alive = np.ones(len(pos), bool); x0 = pos[:, 0].copy()
            for q in range(1, NF + 1):
                if len(PL[i + q]) == 0: alive[:] = False; break
                d, jj = cKDTree(PL[i + q][:, :2]).query(pos); alive &= d < 3; pos = PL[i + q][jj, :2]; dfin = DEP[i + q][jj]
            esc = alive & (dfin > 6); meet = np.abs(x0) < A / 4; tip = np.abs(x0) > A / 2; dt = (t[i + 1] - t[i])
            # each frame is an independent start: rate = escapes per wall-defect start, per 20 t.u. -> per 1000 t.u.
            a = agg[st]; a["esc"] += esc.sum(); a["esc_meet"] += (esc & meet).sum(); a["esc_tip"] += (esc & tip).sum()
            a["wt"] += w.sum() * 20.0; a["wt_meet"] += meet.sum() * 20.0; a["wt_tip"] += tip.sum() * 20.0
            per_int[iv][0] += esc.sum(); per_int[iv][1] += w.sum() * 20.0
        seq_two = np.array([state[tt] for tt in tf]); seq_esc = np.array([per_int[tt][0] / max(per_int[tt][1], 1) * 1000 for tt in tf])
        lagged.append((seq_two, seq_esc))
    say(f"\n== eps {eps} (two seeds)")
    for st, lab in ((0, "one vortex "), (1, "two vortices")):
        a = agg[st]
        if a["wt"] == 0: say(f"   {lab}: no intervals"); continue
        say(f"   {lab}: wall-defect time {a['wt']:.0f}; escape probability within 20 t.u. per wall defect {a['esc'] / a['wt'] * 20:.4f} (n = {a['esc']}) "
            f"(meeting region |x|<A/4: {a['esc_meet'] / max(a['wt_meet'], 1) * 20:.4f}; tips |x|>A/2: {a['esc_tip'] / max(a['wt_tip'], 1) * 20:.4f}); "
            f"mean vacancies {np.mean(a['vac']):.2f}")
    # lagged correlation: two-vortex indicator at interval i vs escape rate at interval i + lag (lag in intervals of 100 t.u.)
    line = "   corr(two-vortex at i, escape rate at i+lag): "
    for lag in (-2, -1, 0, 1, 2):
        xs, ys = [], []
        for two, esc in lagged:
            n = len(two)
            for i in range(n):
                if 0 <= i + lag < n: xs.append(two[i]); ys.append(esc[i + lag])
        xs, ys = np.array(xs, float), np.array(ys, float)
        c = np.corrcoef(xs, ys)[0, 1] if xs.std() > 0 and ys.std() > 0 else np.nan
        line += f"lag {lag:+d}: {c:+.2f}  "
    say(line + "(positive lag = escapes after the vortex state)")
open(os.path.join(H, "vortex_vs_layer_tracked.txt"), "w").write("\n".join(L) + "\n")
