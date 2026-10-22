"""Shared helpers for the SI alternative sheets (2026-10-22)."""
import json, sys, numpy as np
from scipy.spatial import cKDTree
sys.path.insert(0, "/home/baloewe/Physics/Overcharged/runs/cavity")
from edge_criteria import read_with_angle  # noqa
from wall_fourier import wall_polyline  # noqa
RUNS = "/home/baloewe/Physics/Overcharged/runs"; LA = 22.21
CASES = [("Circle, $R=1.43\\,\\ell_a$, $k=6$ ($\\rho=0.06$)", RUNS + "/cavity/server_main/main_act0.07/rho0.06/k06/s05"),
         ("Circle, $R=2.87\\,\\ell_a$, $k=6$ ($\\rho=0.03$)", RUNS + "/cavity/server_main/main_act0.07/rho0.03/k06/s05"),
         ("Ellipse 2, $k=5$ ($\\rho=0.04$)", RUNS + "/shapes_server/ellipse/eps2.00/Req1.5/k05/s01"),
         ("Ellipse 2, $k=16$ ($\\rho=0.14$)", RUNS + "/shapes_server/ellipse/eps2.00/Req1.5/k16/s01")]
def load(run, tfrac=0.5):
    bc = json.load(open(run + "/input.json"))["BC"][0]; W, s, P = wall_polyline(bc); tree = cKDTree(W); Q = np.array(bc["Q"][:2], float)
    t, fr = read_with_angle(run, Q, tfrac); s = np.asarray(s)
    Wc = W[:, 0] + 1j * W[:, 1]; T = np.roll(Wc, -1) - np.roll(Wc, 1); T /= np.abs(T); Nn = -1j * T
    if np.mean((Nn * np.conj(Wc - Wc.mean())).real) < 0: Nn = -Nn
    return dict(bc=bc, W=W, s=s, P=P, tree=tree, t=np.array(t), fr=fr, Nn=Nn)
def tracks(D, kind=+1, RL=3.5, GAP=2):
    """+1/2 (kind > 0) tracks: dict id -> list of (frame, x, y, d, s, tilt)."""
    act = {}; out = {}; nid = 0
    for n, f in enumerate(D["fr"]):
        p = f[np.sign(f[:, 2]) == kind] if len(f) else np.zeros((0, 4)); ids = np.full(len(p), -1)
        live = [(i, xy) for i, (xy, ln) in act.items() if n - ln <= 1 + GAP]
        if live and len(p):
            LX = np.array([xy for _, xy in live]); LI = np.array([i for i, _ in live]); dd, jj = cKDTree(LX).query(p[:, :2]); _, kk = cKDTree(p[:, :2]).query(LX)
            ok = (dd < RL) & (kk[jj] == np.arange(len(p))); ids[ok] = LI[jj[ok]]
        new = ids < 0; ids[new] = nid + np.arange(new.sum()); nid += new.sum()
        if len(p):
            d, j = D["tree"].query(p[:, :2]); tl = np.degrees(np.angle(np.exp(1j * (p[:, 3] - np.angle(D["Nn"][j])))))
            for i, (x, y), dv, jv, tv in zip(ids, p[:, :2], d, j, tl): out.setdefault(int(i), []).append((n, x, y, dv, D["s"][jv], tv)); act[int(i)] = ((x, y), n)
        act = {i: v for i, v in act.items() if n - v[1] <= 1 + GAP}
    return out
