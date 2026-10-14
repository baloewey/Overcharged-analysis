#!/usr/bin/env python3
"""'Duck-duck-goose' in circular cavities (user, 2026-10-14, from movies): a bulk +1/2 orbits along the wall layer, finds a spot, enters, and another wall defect
is pushed into the bulk. Test with tracked identities (kickout/kick_identity.track; frames every 2 t.u.; second half).
Events (hysteresis): ENTRY = a +1/2 track that was deeper than 6 cells reaches the wall band (depth < 4); EXIT = a wall +1/2 that leaves the band and gets deeper than 6.
For every entry: the first exit that follows (any defect) -> delay, whether it is the SAME defect (a bounce) or ANOTHER one (a 'goose'), and the angular distance
between the entry and exit sites in units of the site spacing 2 pi/(2k-2). Null model: entry times shifted by a random offset (same exit sequence).
Orbiting: for each entering +1/2, its azimuthal travel in the 4-8 cell ring over the 40 t.u. before entry (signed by the cavity circulation sense estimated from
the wall-defect tilt sign; reported as |travel| in site spacings) vs the same for bulk +1/2 in the ring that do not enter within 40 t.u.
Usage: ddg.py -> ddg.txt"""
import os, sys, glob, json, numpy as np
from multiprocessing import Pool
from scipy.spatial import cKDTree
R_ = "/home/baloewe/Physics/Overcharged/runs"; H = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, R_ + "/cavity"); sys.path.insert(0, os.path.join(H, "..", "kickout")); sys.path.insert(0, os.path.join(H, "..", "energetics"))
from edge_criteria import read_with_angle  # noqa
from kick_identity import track  # noqa
import energetics as En  # noqa
rng = np.random.default_rng(5)


def one(run):
    J = json.load(open(run + "/input.json")); bc = J["BC"][0]; k = bc["customAnch"]; Q = np.array(bc["Q"][:2], float)
    W, s, P = En.wall_polyline(bc); tree = cKDTree(W)
    t, fr = read_with_angle(run, Q, 0.5); t = np.array(t); dt = float(np.median(np.diff(t)))
    D = [tree.query(f[:, :2])[0] if len(f) else np.zeros(0) for f in fr]
    ids, info = track(fr, D, t)
    traj = {}
    for n, (f, idn) in enumerate(zip(fr, ids)):
        for i, kk in enumerate(idn):
            if f[i, 2] > 0: traj.setdefault(kk, []).append((n, D[n][i], np.arctan2(f[i, 1], f[i, 0])))
    entries, exits, orbit_in, orbit_ctl = [], [], [], []
    site = 2 * np.pi / max(2 * abs(k) - 2, 1)
    for kk, tr in traj.items():
        tr = np.array(tr); state = None
        for j in range(len(tr)):
            d = tr[j, 1]
            if d > 6:
                if state == "wall": exits.append((tr[j, 0], tr[j, 2], kk))
                state = "bulk"
            elif d < 4:
                if state == "bulk":
                    entries.append((tr[j, 0], tr[j, 2], kk))
                    w = tr[(tr[:, 0] >= tr[j, 0] - 20) & (tr[:, 0] < tr[j, 0]) & (tr[:, 1] >= 4) & (tr[:, 1] <= 8)]
                    if len(w) >= 3: orbit_in.append(abs(np.unwrap(w[:, 2])[-1] - np.unwrap(w[:, 2])[0]) / site)
                state = "wall"
        ring = tr[(tr[:, 1] >= 4) & (tr[:, 1] <= 8)]
        if len(ring) > 20 and rng.random() < 0.3:
            a = rng.integers(0, len(ring) - 20); w = ring[a:a + 20]
            if np.all(np.diff(w[:, 0]) <= 3): orbit_ctl.append(abs(np.unwrap(w[:, 2])[-1] - np.unwrap(w[:, 2])[0]) / site)
    exits.sort(); ex_t = np.array([e[0] for e in exits]); out = []
    for (n, ph, kk) in entries:
        j = np.searchsorted(ex_t, n, side="right")
        if j < len(exits):
            ne, phe, ke = exits[j]; dphi = abs(np.angle(np.exp(1j * (phe - ph)))) / site
            out.append(((ne - n) * dt, ke == kk, dphi))
    null = []
    for (n, ph, kk) in entries:
        n2 = rng.integers(0, len(fr)); j = np.searchsorted(ex_t, n2, side="right")
        if j < len(exits): null.append((exits[j][0] - n2) * dt)
    return dict(run=os.path.relpath(run, R_), k=k, T=len(fr) * dt, nent=len(entries), nex=len(exits), pairs=out, null=null, oin=orbit_in, octl=orbit_ctl)


GROUPS = {"circles R_eq 1.0": f"{R_}/shapes_server/circle/circle/Req1.0/k*/s*", "circles R_eq 1.5": f"{R_}/shapes_server/circle/circle/Req1.5/k*/s*",
          "circles R_eq 2.0": f"{R_}/shapes_server/circle/circle/Req2.0/k*/s*", "ellipses 2.0, R_eq 1.5-2": f"{R_}/shapes_server/ellipse/eps2.00/Req[12].[05]/k*/s*"}

if __name__ == "__main__":
    runs = {g: sorted(glob.glob(p)) for g, p in GROUPS.items()}; allr = [r for v in runs.values() for r in v]
    with Pool(int(os.environ.get("NPROC", 6))) as pool: res = dict(zip(allr, pool.map(one, allr, chunksize=1)))
    L = ["DUCK-DUCK-GOOSE test (tracked +1/2; entry = bulk (>6) -> band (<4), exit = band -> >6; second half)",
         "group | runs | entries, exits per wall defect per 1000 t.u. | after an entry, next exit: median delay (t.u.) [null] | P(exit within 20 t.u.) [null] | next exit by ANOTHER defect | its distance from the entry site (site spacings): median, P(<= 1.5) [uniform] | orbit before entry (site spacings, 40 t.u.) vs ring control"]
    for g, rr in runs.items():
        R = [res[r] for r in rr]
        if not R: continue
        P = [p for x in R for p in x["pairs"]]; Nl = np.array([v for x in R for v in x["null"]])
        if not P: continue
        dl = np.array([p[0] for p in P]); same = np.array([p[1] for p in P]); dph = np.array([p[2] for p in P])
        nw = np.mean([2 * abs(x["k"]) - 2 for x in R]); Tt = sum(x["T"] for x in R)
        ent = 1000 * sum(x["nent"] for x in R) / (Tt * nw); exr = 1000 * sum(x["nex"] for x in R) / (Tt * nw)
        oth = ~same; maxd = np.mean([(2 * abs(x["k"]) - 2) / 2 for x in R])
        oin = np.array([v for x in R for v in x["oin"]]); octl = np.array([v for x in R for v in x["octl"]])
        L.append(f"{g:26s} | {len(R):2d} | {ent:.2f}, {exr:.2f} | {np.median(dl):.0f} [{np.median(Nl):.0f}] | {np.mean(dl <= 20):.2f} [{np.mean(Nl <= 20):.2f}] | {np.mean(oth):.2f} | "
                 f"{np.median(dph[oth]) if oth.any() else np.nan:.1f}, {np.mean(dph[oth] <= 1.5) if oth.any() else np.nan:.2f} [{1.5 / maxd:.2f}] | {np.median(oin) if len(oin) else np.nan:.2f} vs {np.median(octl) if len(octl) else np.nan:.2f}")
    print("\n".join(L)); open(os.path.join(H, "ddg.txt"), "w").write("\n".join(L) + "\n")
