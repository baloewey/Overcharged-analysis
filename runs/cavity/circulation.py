#!/usr/bin/env python3
"""
Flow circulation order parameter vs bulk topological charge in circular cavities (H4).

Per run (directory with input.json, flowfield.dat, defects.dat):
  m(t) = < v_phi >_cells          azimuthal flow, averaged over cells with r < R - rcut
  s(t) = < |v| >_cells            mean flow speed (normalisation)
  Phi_persist = | <m>_t | / <s>_t  persistent sense of rotation (user's definition, normalised)
  Phi_inst    = < |m| >_t / <s>_t  instantaneous coherent circulation (insensitive to reversals)
  (raw, un-normalised versions |<m>_t| and <|m|>_t are also reported)
  bulk charge q_b(t) = net charge of tracked defects with r < R - wall
  P1 = fraction of frames with q_b = +1,  <q_b>
  audit (defect_audit.py): fraction of frames with total net charge = k, with
  missing imposed defects, or with unpaired dynamical defects; mean number of
  dynamical ±1/2 pairs. P1_ok / qb_ok: bulk statistics on audit-passing frames only
  (compare with P1 / qb_mean: a difference means the audit filter biases the result).

flowfield.dat frames are averages over the flowOut interval.

Runs in seed directories (…/rho0.10/k05/s01) are grouped into their parameter set
(…/rho0.10/k05); set values are the mean ± standard error over seeds. Any other run
directory is its own set.

Usage:
    python3 circulation.py main_act0.07/rho*/k*/s*/ [--wall 2.5] [--rcut 1.0] [--out circ.csv]
"""
import argparse
import json
import os
import re
from collections import defaultdict

import numpy as np

from defect_audit import audit_frames, audit_summary, open_data, read_frames


def read_input(rundir):
    with open(os.path.join(rundir, "input.json")) as f:
        inp = json.load(f)
    bc = inp["BC"][0]
    ca = bc.get("customAnch", 0)
    R = bc["R"]
    return {"L": inp["domain"][:2], "C": bc["Q"][:2], "R": R, "k": ca if ca else 1,
            "rho": abs(ca if ca else 1) / (np.pi * R), "act": inp["species"][0]["act"]}


def read_flow(path, L):
    """Yield (t, vx[L], vy[L]) per complete frame."""
    t_cur, rows = None, []

    def emit():
        a = np.array(rows)
        if len(a) != L[0] * L[1]:
            return None
        vx = np.zeros(L); vy = np.zeros(L)
        i, j = a[:, 1].astype(int), a[:, 2].astype(int)
        vx[i, j], vy[i, j] = a[:, 4], a[:, 5]
        return t_cur, vx, vy

    with open_data(path) as f:
        for line in f:
            tok = line.split()
            if len(tok) != 7:
                continue
            try:
                v = [float(x) for x in tok]
            except ValueError:
                continue
            if t_cur is not None and v[0] != t_cur:
                fr = emit()
                if fr:
                    yield fr
                rows = []
            t_cur = v[0]
            rows.append(v)
    if rows:
        fr = emit()
        if fr:
            yield fr


def analyse(rundir, wall, rcut, tmin):
    p = read_input(rundir)
    L, C, R = p["L"], p["C"], p["R"]
    ii, jj = np.meshgrid(np.arange(L[0]) + 0.5, np.arange(L[1]) + 0.5, indexing="ij")
    x, y = ii - C[0], jj - C[1]
    r = np.hypot(x, y)
    mask = (r < R - rcut) & (r > 0.5)
    ephx, ephy = -y[mask] / r[mask], x[mask] / r[mask]

    m, s = [], []
    for t, vx, vy in read_flow(os.path.join(rundir, "flowfield.dat"), L):
        if t < tmin:
            continue
        vxm, vym = vx[mask], vy[mask]
        m.append(np.mean(vxm * ephx + vym * ephy))
        s.append(np.mean(np.hypot(vxm, vym)))
    m, s = np.array(m), np.array(s)

    td, frames = read_frames(os.path.join(rundir, "defects.dat"))
    frames = [f for t, f in zip(td, frames) if t >= tmin]
    qb = np.array([0.5 * np.sign(f[np.hypot(f[:, 0] - C[0], f[:, 1] - C[1]) < R - wall, 2]).sum()
                   if len(f) else 0.0 for f in frames])
    aud = audit_frames(frames, p["k"])

    res = {"run": rundir.rstrip("/"), **p, "flow_frames": len(m), "def_frames": len(qb)}
    if len(m):
        res.update({"m_persist": abs(m.mean()), "m_inst": np.abs(m).mean(),
                    "speed": s.mean(), "Phi_persist": abs(m.mean()) / s.mean(),
                    "Phi_inst": np.abs(m).mean() / s.mean()})
    if len(qb):
        ok = aud["ok_net"]
        res.update({"P1": np.mean(np.abs(qb - 1.0) < 0.25), "qb_mean": qb.mean(),
                    "P1_ok": np.mean(np.abs(qb[ok] - 1.0) < 0.25) if ok.any() else np.nan,
                    "qb_ok": qb[ok].mean() if ok.any() else np.nan,
                    **audit_summary(aud)})
    return res


KEYS = ["Phi_persist", "Phi_inst", "m_persist", "m_inst", "speed", "P1", "qb_mean",
        "P1_ok", "qb_ok", "frac_ok", "frac_missing_imposed", "frac_unpaired", "N_pairs_mean"]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--wall", type=float, default=2.5,
                    help="bulk = defects with r < R - wall (default 2.5 cells)")
    ap.add_argument("--rcut", type=float, default=1.0,
                    help="flow average over cells with r < R - rcut")
    ap.add_argument("--tmin", type=float, default=0.0,
                    help="ignore frames before this time (warm-up is already excluded)")
    ap.add_argument("--out", default="circulation.csv")
    a = ap.parse_args()

    runs = [analyse(d, a.wall, a.rcut, a.tmin) for d in a.runs]
    sets = defaultdict(list)
    for r in runs:
        # group seeds (…/k05/s01) into their parameter set; otherwise each run is its own set
        base = os.path.basename(r["run"])
        key = os.path.dirname(r["run"]) if re.fullmatch(r"s\d+", base) else r["run"]
        sets[key].append(r)

    rows = []
    for key, rr in sets.items():
        row = {"set": key, "k": rr[0]["k"], "R": rr[0]["R"], "rho": rr[0]["rho"],
               "act": rr[0]["act"], "nseeds": len(rr)}
        for kname in KEYS:
            v = np.array([r[kname] for r in rr if kname in r])
            row[kname] = v.mean() if len(v) else np.nan
            row[kname + "_err"] = v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else np.nan
        rows.append(row)
    rows.sort(key=lambda r: (r["rho"], r["k"]))

    cols = ["set", "k", "R", "rho", "act", "nseeds"] + [c for k in KEYS for c in (k, k + "_err")]
    with open(a.out, "w") as f:
        f.write(",".join(cols) + "\n")
        for r in rows:
            f.write(",".join(str(r[c]) for c in cols) + "\n")

    print(f"{'set':<34} {'k':>4} {'R':>6} {'rho':>5} {'n':>2} {'Phi_pers':>9} {'Phi_inst':>9} "
          f"{'P(qb=1)':>8} {'<qb>':>6} | {'q=k':>5} {'miss':>5} {'unpr':>5} {'pairs':>5}")
    for r in rows:
        print(f"{r['set']:<34} {r['k']:>4} {r['R']:6.1f} {r['rho']:5.2f} {r['nseeds']:>2} "
              f"{r['Phi_persist']:9.3f} {r['Phi_inst']:9.3f} {r['P1']:8.2f} {r['qb_mean']:6.2f} | "
              f"{r['frac_ok']:5.2f} {r['frac_missing_imposed']:5.2f} {r['frac_unpaired']:5.2f} "
              f"{r['N_pairs_mean']:5.2f}")
    print(f"\nwritten {a.out}  (set values: mean over seeds; *_err = standard error)\n"
          "audit: q=k frames with total charge = k; miss = imposed defects missing; "
          "unpr = unpaired dynamical defects; pairs = mean dynamical ±1/2 pairs")


if __name__ == "__main__":
    main()
