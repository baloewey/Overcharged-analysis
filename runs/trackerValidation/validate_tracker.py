#!/usr/bin/env python3
"""
Validate the default (PBC-aware, union-find) defect tracker against
  (a) an independent Python re-implementation working from directorfield.dat, and
  (b) the legacy tracker (same seed -> same dynamics).

Checks
  1. dynamics identical for defectTracker 1 and 0 (directorfield.dat byte-identical)
  2. charge field: new == legacy on every cell where legacy computes it (topochargefield.dat)
  3. new tracker == Python re-implementation, frame by frame (number of defects, charges,
     positions matched within 1e-3 cells, minimum image)
  4. net charge per frame: new vs legacy (PBC: must be 0)
  5. defects near the periodic seam: found by the new tracker, missed by the legacy one

Usage: validate_tracker.py NEW_RUN_DIR LEGACY_RUN_DIR [--tolD 0.01]
"""
import argparse
import filecmp
import json
import os

import numpy as np


def read_director(path, L):
    frames, ts, cur_t, rows = [], [], None, []
    with open(path) as f:
        for line in f:
            tok = line.split()
            if len(tok) != 8:
                continue
            try:
                v = [float(x) for x in tok]
            except ValueError:
                continue
            if cur_t is not None and v[0] != cur_t:
                frames.append(np.array(rows)); ts.append(cur_t); rows = []
            cur_t = v[0]; rows.append(v)
    if rows:
        frames.append(np.array(rows)); ts.append(cur_t)
    out = []
    for t, a in zip(ts, frames):
        if len(a) != L[0] * L[1]:
            continue
        ang = np.zeros(L)
        ang[a[:, 1].astype(int), a[:, 2].astype(int)] = np.arctan2(a[:, 5], a[:, 4])
        out.append((t, ang))
    return out


def read_defects(path):
    t, frames, cur = [], [], None
    for line in open(path):
        tok = line.split()
        try:
            v = [float(x) for x in tok]
        except ValueError:
            continue
        if len(v) == 2:
            if cur is not None:
                frames.append(np.array(cur).reshape(-1, 3))
            t.append(v[0]); cur = []
        elif len(v) == 4 and cur is not None:
            cur.append(v[:3])
    if cur is not None:
        frames.append(np.array(cur).reshape(-1, 3))
    return dict(zip(t, frames))


def read_topofield(path, L):
    out = {}
    for line in open(path):
        tok = line.split()
        if len(tok) != 6:
            continue
        try:
            v = [float(x) for x in tok]
        except ValueError:
            continue
        out.setdefault(v[0], np.zeros(L))[int(v[1]), int(v[2])] = v[4]
    return out


def charge_field(ang, pbc):
    """Winding of the nematic director around the 8 neighbours (same ring as the solver)."""
    di = [1, 1, 0, -1, -1, -1, 0, 1]
    dj = [0, 1, 1, 1, 0, -1, -1, -1]
    ring = [np.roll(np.roll(ang, -a, 0), -b, 1) for a, b in zip(di, dj)]
    tot = np.zeros_like(ang)
    for m in range(8):
        d = ring[(m + 1) % 8] - ring[m]
        d = (d + np.pi / 2) % np.pi - np.pi / 2          # nematic: (-pi/2, pi/2]
        tot += d
    q = tot / (2 * np.pi)
    if not pbc:
        q[0, :] = q[-1, :] = q[:, 0] = q[:, -1] = 0.0
    return q


def cluster(q, tol, pbc):
    L = q.shape
    cand = np.abs(np.abs(q) - 0.5) < tol
    idx = {c: n for n, c in enumerate(zip(*np.nonzero(cand)))}
    parent = list(range(len(idx)))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for (i, j), n in idx.items():
        for a in (-1, 0, 1):
            for b in (-1, 0, 1):
                if a == b == 0:
                    continue
                ii, jj = i + a, j + b
                if pbc:
                    ii %= L[0]; jj %= L[1]
                if (ii, jj) in idx and q[i, j] * q[ii, jj] > 0:
                    ra, rb = find(n), find(idx[(ii, jj)])
                    if ra != rb:
                        parent[ra] = rb
    groups = {}
    for c, n in idx.items():
        groups.setdefault(find(n), []).append(c)
    out = []
    for cells in groups.values():
        c = np.array(cells, float)
        ref = c[0]
        d = c - ref
        if pbc:
            d -= np.array(L) * np.round(d / np.array(L))
        pos = ref + 0.5 + d.mean(0)
        if pbc:
            pos %= np.array(L)
        out.append((pos[0], pos[1], np.mean([q[int(x), int(y)] for x, y in cells])))
    return np.array(out).reshape(-1, 3)


def match(a, b, L, tol=1e-3):
    """Every defect in a has one in b with the same sign within tol (minimum image)."""
    if len(a) != len(b):
        return False
    used = set()
    for x, y, q in a:
        d = b[:, :2] - [x, y]
        d -= np.array(L) * np.round(d / np.array(L))
        dist = np.hypot(d[:, 0], d[:, 1])
        ok = [n for n in np.argsort(dist) if n not in used and np.sign(b[n, 2]) == np.sign(q)]
        if not ok or dist[ok[0]] > tol:
            return False
        used.add(ok[0])
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("new"); ap.add_argument("legacy")
    ap.add_argument("--tolD", type=float, default=0.01)
    a = ap.parse_args()
    inp = json.load(open(os.path.join(a.new, "input.json")))
    L = tuple(inp["domain"][:2])
    pbc = inp.get("domainWalls", 0) == 1 and "BC" not in inp

    same = filecmp.cmp(os.path.join(a.new, "directorfield.dat"),
                       os.path.join(a.legacy, "directorfield.dat"), shallow=False)
    print(f"1. dynamics identical (directorfield.dat byte-identical): {same}")

    tn = read_topofield(os.path.join(a.new, "topochargefield.dat"), L)
    tl = read_topofield(os.path.join(a.legacy, "topochargefield.dat"), L)
    maxdiff = max(np.abs(tn[t][1:-1, 1:-1] - tl[t][1:-1, 1:-1]).max() for t in tn if t in tl)
    edge = sum(int((np.abs(tn[t]) > 0.25).sum() - (np.abs(tn[t][1:-1, 1:-1]) > 0.25).sum()) for t in tn)
    print(f"2. charge field, interior cells: max |new - legacy| = {maxdiff:.3g} over {len(tn)} frames;"
          f" charged cells in the outer layer (new only): {edge}")

    dn = read_defects(os.path.join(a.new, "defects.dat"))
    dl = read_defects(os.path.join(a.legacy, "defects.dat"))
    frames = read_director(os.path.join(a.new, "directorfield.dat"), L)
    nmatch = ncount = 0
    bad = []
    for t, ang in frames:
        if t not in dn:
            continue
        py = cluster(charge_field(ang, pbc), a.tolD, pbc)
        ncount += 1
        if match(dn[t], py, L):
            nmatch += 1
        else:
            bad.append((t, len(dn[t]), len(py)))
    print(f"3. new tracker == Python re-implementation: {nmatch}/{ncount} frames"
          + (f"; first mismatches (t, N_new, N_py): {bad[:5]}" if bad else ""))

    common = sorted(set(dn) & set(dl))
    qn = np.array([0.5 * np.sign(dn[t][:, 2]).sum() if len(dn[t]) else 0 for t in common])
    ql = np.array([0.5 * np.sign(dl[t][:, 2]).sum() if len(dl[t]) else 0 for t in common])
    Nn = np.array([len(dn[t]) for t in common]); Nl = np.array([len(dl[t]) for t in common])
    print(f"4. frames with net charge 0: new {np.mean(np.abs(qn) < .25):.3f}, legacy {np.mean(np.abs(ql) < .25):.3f}"
          f"  (<N> new {Nn.mean():.2f}, legacy {Nl.mean():.2f}; {len(common)} frames)")

    seam = 0
    for t in common:
        if len(dn[t]):
            x, y = dn[t][:, 0], dn[t][:, 1]
            seam += int(((x < 1.5) | (x > L[0] - 1.5) | (y < 1.5) | (y > L[1] - 1.5)).sum())
    print(f"5. defects within 1.5 cells of the box edge found by the new tracker: {seam} "
          f"(legacy can't see these)")


if __name__ == "__main__":
    main()
