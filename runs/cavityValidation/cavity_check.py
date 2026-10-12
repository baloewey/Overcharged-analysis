#!/usr/bin/env python3
"""
Consistency checks for circular-cavity runs with custom (winding k) anchoring.

For each run directory it reads input.json, directorfield.dat and defects.dat:

  From the director field (independent of the solver's defect tracker):
    Q(r)     enclosed topological charge = winding of the director around a
             circle of radius r (Q-tensor interpolated bilinearly, nematic
             angle differences wrapped to (-pi/2, pi/2]).
    align    <cos 2(psi - psi_anch)> just inside the wall, where the imposed
             homeotropic custom anchoring is psi_anch = k*phi (mod pi).
             1 = perfect anchoring, 0 = uncorrelated.

  From defects.dat (the solver's tracker):
    N, N+, N-, net charge q = (N+ - N-)/2, number of defects near the wall.

Expected for winding k: Q(r -> R) = k and net tracker charge = k.

Usage:
    python3 cavity_check.py passive/k*/ [--tfrac 0.5] [--wall 2.0]
"""
import argparse
import json
import os

import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "cavity"))
from defect_audit import audit_frames, audit_summary, open_data  # noqa: E402


def read_input(rundir):
    with open(os.path.join(rundir, "input.json")) as f:
        inp = json.load(f)
    bc = inp["BC"][0]
    k = bc.get("customAnch", 0)
    return {"L": inp["domain"][:2], "C": bc["Q"][:2], "R": bc["R"],
            "k": k if k != 0 else 1, "ca": k, "act": inp["species"][0].get("act", 0)}


def read_director(path, L):
    """Return times and arrays (nframes, Lx, Ly) of Qxx, Qxy (S-weighted) and S."""
    t_list, frames = [], []
    cur_t, buf = None, []
    with open_data(path) as f:
        for line in f:
            tok = line.split()
            if len(tok) != 8:
                continue
            try:
                v = [float(x) for x in tok]
            except ValueError:
                continue
            if cur_t is None or v[0] != cur_t:
                if buf:
                    frames.append(np.array(buf))
                    t_list.append(cur_t)
                cur_t, buf = v[0], []
            buf.append(v)
    if buf:
        frames.append(np.array(buf))
        t_list.append(cur_t)
    out_t, Qxx, Qxy, S = [], [], [], []
    for t, a in zip(t_list, frames):
        if len(a) != L[0] * L[1]:
            continue                                  # truncated frame
        i, j = a[:, 1].astype(int), a[:, 2].astype(int)
        nx, ny, s = a[:, 4], a[:, 5], a[:, 7]
        qxx = np.zeros(L); qxy = np.zeros(L); ss = np.zeros(L)
        qxx[i, j] = s * (nx * nx - ny * ny)           # 2D Q = S(2nn - I)/... up to a factor
        qxy[i, j] = s * (2 * nx * ny)
        ss[i, j] = s
        out_t.append(t); Qxx.append(qxx); Qxy.append(qxy); S.append(ss)
    return np.array(out_t), np.array(Qxx), np.array(Qxy), np.array(S)


def bilinear(field, x, y):
    """Bilinear interpolation of a cell-centred field (cell i centred at i+0.5)."""
    gx, gy = x - 0.5, y - 0.5
    i0, j0 = np.floor(gx).astype(int), np.floor(gy).astype(int)
    fx, fy = gx - i0, gy - j0
    Lx, Ly = field.shape
    i0c, i1c = np.clip(i0, 0, Lx - 1), np.clip(i0 + 1, 0, Lx - 1)
    j0c, j1c = np.clip(j0, 0, Ly - 1), np.clip(j0 + 1, 0, Ly - 1)
    return ((1 - fx) * (1 - fy) * field[i0c, j0c] + fx * (1 - fy) * field[i1c, j0c]
            + (1 - fx) * fy * field[i0c, j1c] + fx * fy * field[i1c, j1c])


def director_angle_on_circle(qxx, qxy, C, r, M):
    phi = 2 * np.pi * np.arange(M) / M
    x, y = C[0] + r * np.cos(phi), C[1] + r * np.sin(phi)
    a, b = bilinear(qxx, x, y), bilinear(qxy, x, y)
    return phi, 0.5 * np.arctan2(b, a)


def winding(psi):
    d = np.diff(np.append(psi, psi[0]))
    d = (d + np.pi / 2) % np.pi - np.pi / 2          # nematic: wrap to [-pi/2, pi/2)
    return d.sum() / (2 * np.pi)


def read_defects(path):
    t, frames, cur = [], [], None
    with open_data(path) as f:
        for line in f:
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
    return np.array(t), frames


def analyse(rundir, tfrac, wall):
    p = read_input(rundir)
    L, C, R, k = p["L"], p["C"], p["R"], p["k"]
    out = {"run": rundir.rstrip("/"), **p}

    # ---- director field
    t, Qxx, Qxy, S = read_director(os.path.join(rundir, "directorfield.dat"), L)
    keep = t >= tfrac * t.max()
    radii = np.arange(1.0, R - 0.49, 1.0)
    Qr = np.array([[winding(director_angle_on_circle(Qxx[f], Qxy[f], C, r,
                                                     max(64, int(16 * r)))[1])
                    for r in radii] for f in np.nonzero(keep)[0]])
    out["radii"], out["Qr_mean"], out["Qr_std"] = radii, Qr.mean(0), Qr.std(0)
    out["Qr_final"] = Qr[-1]
    # anchoring quality just inside the wall
    r_in = R - 0.75
    al = []
    for f in np.nonzero(keep)[0]:
        phi, psi = director_angle_on_circle(Qxx[f], Qxy[f], C, r_in, 256)
        al.append(np.mean(np.cos(2 * (psi - k * phi))))
    out["align"] = np.mean(al)
    out["S_bulk"] = np.mean([S[f][S[f] > 0].mean() for f in np.nonzero(keep)[0]])
    out["dir_frames"] = int(keep.sum())

    # ---- solver defect tracker
    td, fr = read_defects(os.path.join(rundir, "defects.dat"))
    kd = td >= tfrac * td.max()
    N = np.array([len(x) for x in fr])[kd]
    Np = np.array([(x[:, 2] > 0).sum() if len(x) else 0 for x in fr])[kd]
    Nm = N - Np
    q = 0.5 * (Np - Nm)
    rdef = [np.hypot(x[:, 0] - C[0], x[:, 1] - C[1]) if len(x) else np.array([]) for x in fr]
    nwall = np.array([(r > R - wall).sum() for r in rdef])[kd]
    aud = audit_summary(audit_frames([x for x, keep in zip(fr, kd) if keep], k))
    out.update({"miss": aud["frac_missing_imposed"], "unpaired": aud["frac_unpaired"],
                "pairs": aud["N_pairs_mean"]})
    out.update({"N": N.mean(), "Nplus": Np.mean(), "Nminus": Nm.mean(),
                "q_mean": q.mean(), "q_frac_ok": np.mean(np.abs(q - k) < 0.25),
                "N_wall": nwall.mean(), "def_frames": int(kd.sum())})
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--tfrac", type=float, default=0.5,
                    help="use frames with t >= tfrac * t_max (default: second half)")
    ap.add_argument("--wall", type=float, default=2.0,
                    help="defects with r > R - wall count as 'at the wall'")
    args = ap.parse_args()

    res = [analyse(d, args.tfrac, args.wall) for d in args.runs]
    print(f"{'run':<24} {'CA':>3} {'k':>3} {'Q(R-1)':>7} {'align':>6} {'S':>5} | "
          f"{'N':>6} {'N+':>6} {'N-':>6} {'q':>6} {'q=k':>5} {'miss':>5} {'unpr':>5} {'pairs':>5} {'Nwall':>6}")
    for r in res:
        print(f"{r['run']:<24} {r['ca']:>3} {r['k']:>3} {r['Qr_mean'][-1]:7.2f} "
              f"{r['align']:6.3f} {r['S_bulk']:5.2f} | {r['N']:6.2f} {r['Nplus']:6.2f} "
              f"{r['Nminus']:6.2f} {r['q_mean']:6.2f} {r['q_frac_ok']:5.2f} {r['miss']:5.2f} "
              f"{r['unpaired']:5.2f} {r['pairs']:5.2f} {r['N_wall']:6.2f}")
    print("\nQ(r) (mean over frames; r in cells):")
    print(f"{'run':<24} " + " ".join(f"{x:5.0f}" for x in res[0]["radii"]))
    for r in res:
        print(f"{r['run']:<24} " + " ".join(f"{x:5.2f}" for x in r["Qr_mean"]))
    print("\nQ(R-1) from director field and q from the tracker should both equal k."
          "\nq=k: fraction of frames where tracker net charge equals k;"
          f" Nwall: defects within {args.wall} cells of the wall.")


if __name__ == "__main__":
    main()
