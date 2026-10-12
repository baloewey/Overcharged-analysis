"""Restoring force on WALL -1/2 (negative winding; not self-propelled) in the P6 mirror shapes, same method as restoring_swirl.py (2 t.u. steps, second half),
compared with wall +1/2 in the positive-winding runs of the same family. -> restoring_mirror.txt"""
import glob, json, os, sys, numpy as np
from scipy.spatial import cKDTree
sys.path.insert(0, "/home/baloewe/Physics/Overcharged/runs/cavity")
from edge_criteria import read_with_angle  # noqa
from wall_fourier import wall_polyline  # noqa
edges = np.array([0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0])
def samples(run, sign):
    J = json.load(open(run + "/input.json")); bc = J["BC"][0]; Q = bc["Q"]
    if "aInv" not in bc: bc["aInv"] = [1, 1, 1]
    W, s, P = wall_polyline(bc); tree = cKDTree(W); t, fr = read_with_angle(run, np.array(Q[:2], float), 0.5)
    PL = [f[np.sign(f[:, 2]) == sign] if len(f) else np.zeros((0, 4)) for f in fr]; out = []
    for a in range(len(PL) - 1):
        p, q = PL[a], PL[a + 1]
        if not len(p) or not len(q): continue
        d, _ = tree.query(p[:, :2]); w = d < 4
        if not w.any(): continue
        dq, k = cKDTree(q[:, :2]).query(p[w, :2]); dn, _ = tree.query(q[k, :2]); ok = dq < 3
        out += [(dd, n - dd) for dd, n in zip(d[w][ok], dn[ok])]
    return np.array(out), J.get("species", [{}])[0].get("act", None), bc.get("customAnch")
L = ["force -U'/kT_eff per depth bin " + " ".join(f"[{a},{b})" for a, b in zip(edges[:-1], edges[1:]))]
for run in sorted(glob.glob("/home/baloewe/Physics/Overcharged/runs/cassini/P6_mirror_shapes/*") + glob.glob("/home/baloewe/Physics/Overcharged/runs/cassini/P6_mirror_shapes/*/s*")):
    if not os.path.exists(run + "/input.json"): continue
    try:
        k = json.load(open(run + "/input.json"))["BC"][0].get("customAnch", 0); a, act, _ = samples(run, -1 if k < 0 else 1)
    except Exception as ex: L.append(f"{run}: {ex}"); continue
    F = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        x = a[(a[:, 0] >= lo) & (a[:, 0] < hi), 1] if len(a) else np.array([]); F.append(2 * x.mean() / x.var() if len(x) > 40 else np.nan)
    L.append(f"{os.path.relpath(run, "/home/baloewe/Physics/Overcharged/runs/cassini/P6_mirror_shapes"):30s} k {k:+d} act {act}: F " + " ".join(f"{x:+6.2f}" for x in F) + f"  (n {len(a)})")
print("\n".join(L)); open("restoring_mirror.txt", "w").write("\n".join(L) + "\n")
