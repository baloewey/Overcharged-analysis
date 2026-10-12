#!/usr/bin/env python3
"""Elastic escape barrier of a wall +1/2 when a bulk +1/2 approaches (cf. the SI model of arXiv:2503.19664, here for an anchored wall).
Anchored disk R = 31.83, k = 6, Frank energy relaxed outside rigid natural cores (ejection/tilted_structures.energy), one-constant 32.5 and anisotropic (42, 23.1).
Configuration: ring of 10 +1/2 at depth 2.06 (sites at 2 pi m/10), wall defect j at angle 0; bulk charge +1 = one fixed +1/2 at (-0.45 R, 0) plus the APPROACHER.
Approacher positions: in front of j (radial, distance r = 3, 4, 5, 6, 8, 10, 14 from j's site), to the side (depth 6, half a site spacing along the wall), far (0, 0.45 R; reference).
For each: energy of j over a grid (depth 1-9, tangential offset -4..4 cells), the rest of the ring fixed; barrier = minimax path cost from the rest site to depth >= 8
(4-neighbour grid moves), relative to the rest-site energy. Writes elastic_kick.txt."""
import heapq, os, sys, time
import numpy as np
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(os.path.dirname(H), "ejection"))
from tilted_structures import energy  # noqa
R = 31.830989; N = 10; d0 = 2.06; rs = R - d0
depths = np.arange(1.0, 9.01, 1.0); offs = np.arange(-4.0, 4.01, 2.0)
ring = [rs * np.exp(2j * np.pi * m / N) for m in range(N)]; fixed = -0.45 * R + 0j
def zpos(depth, off): return (R - depth) * np.exp(1j * off / (R - depth))
def barrier(app, K1, K3):
    Emap = np.full((len(depths), len(offs)), np.inf)
    for a, dp in enumerate(depths):
        for b, of in enumerate(offs):
            zj = zpos(dp, of)
            if app is not None and abs(zj - app) < 1.5: continue
            zs = [zj] + ring[1:] + [fixed] + ([app] if app is not None else [0.45j * R])
            Emap[a, b] = energy(zs, [0.5] * len(zs), [0.0] * len(zs), R, K1, K3)
    i0 = int(np.argmin(np.abs(depths - d0))); j0 = int(np.argmin(np.abs(offs))); E0 = Emap[i0, j0]
    best = np.full(Emap.shape, np.inf); best[i0, j0] = E0; pq = [(E0, i0, j0)]
    while pq:
        c, i, j = heapq.heappop(pq)
        if c > best[i, j]: continue
        for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            a, b = i + di, j + dj
            if 0 <= a < len(depths) and 0 <= b < len(offs) and np.isfinite(Emap[a, b]):
                nc = max(c, Emap[a, b])
                if nc < best[a, b]: best[a, b] = nc; heapq.heappush(pq, (nc, a, b))
    out = {ex: best[depths >= ex].min() - E0 for ex in (5, 6, 8)}
    return out, Emap - E0
if __name__ == "__main__":
    L = []
    def say(x): print(x, flush=True); L.append(x)
    cases = [("far (reference)", None)] + [(f"in front, r = {r}", rs - r + 0j) for r in (3, 4, 5, 6, 8, 10, 14)] + \
            [("side, depth 6, half a site along the wall", zpos(6.0, np.pi * (R - 6) / N))]
    for lab, K1, K3 in (("one-constant", 32.5, 32.5), ("anisotropic", 42.0, 23.1)):
        ref = None; maps = {}
        for name, app in cases:
            t0 = time.time(); b, Em = barrier(app, K1, K3); ref = b if ref is None else ref; maps[name] = Em
            say(f"{lab:12s} approacher {name:42s}: barrier to depth 5 / 6 / 8: " + " / ".join(f"{b[ex]:5.1f} ({b[ex] - ref[ex]:+5.1f})" for ex in (5, 6, 8)) + f" kBT  [{time.time() - t0:.0f} s]")
        np.savez(os.path.join(H, f"elastic_kick_maps_{lab}.npz"), depths=depths, offs=offs, **{k.replace(' ', '_').replace(',', '').replace('=', ''): v for k, v in maps.items()})
    open(os.path.join(H, "elastic_kick.txt"), "w").write("\n".join(L) + "\n")
