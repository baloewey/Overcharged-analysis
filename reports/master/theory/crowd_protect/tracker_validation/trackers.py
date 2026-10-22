#!/usr/bin/env python3
"""The three defect linkers used in the project, in one place for validation (2026-10-22). Input: frames, each an array (x, y, charge, angle);
output: list of tracks dict(sg, n=[frame indices], xy=[(x, y)]). Same-sign linking only. Optional periodic x (channels, period L).
  NEW   : predictive (last displacement; v2: cost = min(distance to prediction, distance to last position)), global Hungarian assignment, gate (cells), up to `gap` missed frames (crowd_protect/exchange/exchange_count2.py,
          channel/crossing2/tracker2.py).
  NN    : mutual nearest neighbour within rl cells of the last position, up to `gap` missed frames, no prediction (crowd_protect/exchange/exchange_count.py).
  CHAIN : the escape follow-up of shape_runs.csv / events.npz: frames subsampled every `sub` frames, each start followed by the nearest same-sign defect
          within `rl` cells (not mutual), no prediction (shapes/mirror_shape_analysis.py). Returned as tracks on the subsampled frames (frame indices in
          the original numbering), built greedily from every defect of every subsampled frame (a start-anchored chain, so tracks overlap)."""
import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.spatial import cKDTree
def _d(a, b, L):
    dx = a[:, None, 0] - b[None, :, 0]
    if L: dx = (dx + L / 2) % L - L / 2
    return np.hypot(dx, a[:, None, 1] - b[None, :, 1])
def track_new(fr, gate=9.0, gap=2, L=None, dual=True):
    """dual: cost = min(distance to the predicted position, distance to the last position) (v2, 2026-10-22: slow jittering wall defects match
    where they are, fast ones where they are going; v1 used the prediction alone, which propagated wall swaps along crowded rings)."""
    T = []; act = {+1: [], -1: []}
    for n, f in enumerate(fr):
        for sg in (+1, -1):
            p = f[np.sign(f[:, 2]) == sg][:, :2] if len(f) else np.zeros((0, 2)); A = act[sg]; used = np.zeros(len(p), bool)
            if A and len(p):
                pred = []
                for ti in A:
                    tr = T[ti]; m = n - tr["n"][-1]; v = np.zeros(2)
                    if len(tr["n"]) > 1:
                        v = (np.array(tr["xy"][-1]) - np.array(tr["xy"][-2])) / (tr["n"][-1] - tr["n"][-2])
                        if L: v[0] = (v[0] + L / 2) % L - L / 2
                    pred.append(np.array(tr["xy"][-1]) + v * m)
                C = _d(np.array(pred), p, L)
                if dual: C = np.minimum(C, _d(np.array([T[ti]["xy"][-1] for ti in A]), p, L))
                r, c = linear_sum_assignment(np.where(C < gate, C, 1e6))
                for a, b in zip(r, c):
                    if C[a, b] < gate: T[A[a]]["n"].append(n); T[A[a]]["xy"].append(tuple(p[b])); used[b] = True
            keep = [ti for ti in A if n - T[ti]["n"][-1] <= gap]
            for b in np.flatnonzero(~used): T.append(dict(sg=sg, n=[n], xy=[tuple(p[b])])); keep.append(len(T) - 1)
            act[sg] = keep
    return T
def track_nn(fr, rl=3.5, gap=2, L=None):
    T = []; act = {+1: [], -1: []}
    for n, f in enumerate(fr):
        for sg in (+1, -1):
            p = f[np.sign(f[:, 2]) == sg][:, :2] if len(f) else np.zeros((0, 2)); A = act[sg]; used = np.zeros(len(p), bool)
            if A and len(p):
                last = np.array([T[ti]["xy"][-1] for ti in A]); C = _d(p, last, L); jj = C.argmin(1); dd = C[np.arange(len(p)), jj]; kk = C.argmin(0)
                for b in range(len(p)):
                    if dd[b] < rl and kk[jj[b]] == b: T[A[jj[b]]]["n"].append(n); T[A[jj[b]]]["xy"].append(tuple(p[b])); used[b] = True
            keep = [ti for ti in A if n - T[ti]["n"][-1] <= gap]
            for b in np.flatnonzero(~used): T.append(dict(sg=sg, n=[n], xy=[tuple(p[b])])); keep.append(len(T) - 1)
            act[sg] = keep
    return T
def chain_follow(fr, start_n, xy0, sg, steps=2, sub=5, rl=6.0, L=None):
    """CHAIN: position after `steps` subsampled frames of a start (None if the chain breaks)."""
    pos = np.array(xy0, float)
    for b in range(1, steps + 1):
        n = start_n + b * sub
        if n >= len(fr) or not len(fr[n]): return None
        q = fr[n][np.sign(fr[n][:, 2]) == sg][:, :2]
        if not len(q): return None
        d = _d(pos[None], q, L)[0]; j = d.argmin()
        if d[j] >= rl: return None
        pos = q[j]
    return pos
