"""Restoring force on wall +1/2 in the imposed-swirl runs at full time resolution (frames every 2 t.u.; depth change over 2 t.u.), second half.
Passive vs active at similar coherence. -> restoring_swirl.txt"""
import glob, os, numpy as np
import restoring_coh as RC
from scipy.spatial import cKDTree
def samples2(run):
    import json, sys
    J = json.load(open(run + "/input.json")); bc = J["BC"][0]; Q = bc["Q"]
    if "aInv" not in bc: bc["aInv"] = [1, 1, 1]
    W, s, P = RC.wall_polyline(bc); tree = cKDTree(W); Wc = W[:, 0] + 1j * W[:, 1]; Tn = np.roll(Wc, -1) - np.roll(Wc, 1); Nn = -1j * Tn / np.abs(Tn)
    if np.mean((Nn * np.conj(Wc)).real) < 0: Nn = -Nn
    t, fr = RC.read_with_angle(run, np.array(Q[:2], float), 0.5); PL = [f[f[:, 2] > 0] if len(f) else np.zeros((0, 4)) for f in fr]; out = []; coh = []
    for a in range(len(PL) - 1):
        p, q = PL[a], PL[a + 1]
        if not len(p) or not len(q): continue
        d, j = tree.query(p[:, :2]); w = d < 4
        if not w.any(): continue
        coh.append(abs(np.mean(np.exp(1j * (p[w, 3] - np.angle(Nn[j[w]]))))))
        dq, k = cKDTree(q[:, :2]).query(p[w, :2]); dn, _ = tree.query(q[k, :2]); ok = dq < 3
        out += [(dd, n - dd) for dd, n in zip(d[w][ok], dn[ok])]
    return np.array(out), np.mean(coh)
edges = np.array([0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0])
if __name__ == "__main__":
    L = []
    L.append("force -U'/kT_eff (1/cell) per depth bin " + " ".join(f"[{a},{b})" for a, b in zip(edges[:-1], edges[1:])) + " ; well = U(3.5) - U_min")
    for run in sorted(glob.glob("/home/baloewe/Physics/Overcharged/runs/swirl_ellipse*/*") + glob.glob("/home/baloewe/Physics/Overcharged/runs/swirl_threshold/passive_A0.01_s*")):
        if not os.path.exists(run + "/DONE"): continue
        a, coh = samples2(run); F = []; nn = []
        for lo, hi in zip(edges[:-1], edges[1:]):
            x = a[(a[:, 0] >= lo) & (a[:, 0] < hi), 1]; nn.append(len(x)); F.append(2 * x.mean() / x.var() if len(x) > 40 else np.nan)
        F = np.array(F); w = np.diff(edges); cen = 0.5 * (edges[1:] + edges[:-1]); U = np.concatenate([[0], np.cumsum(-0.5 * (F[:-1] + F[1:]) * np.diff(cen))]); U -= np.nanmin(U)
        L.append(f"{os.path.relpath(run, '/home/baloewe/Physics/Overcharged/runs'):42s} coh {coh:.2f}: F " + " ".join(f"{x:+6.2f}" for x in F) + f" | well {U[-1]:5.2f} | n deep(>2) {sum(nn[4:])}")
    print("\n".join(L)); open("restoring_swirl.txt", "w").write("\n".join(L) + "\n")
