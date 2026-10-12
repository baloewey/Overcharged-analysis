#!/usr/bin/env python3
"""
Active length scale from MPCD defect data.

For each run directory (containing defects.dat) this:
  1. reads the defects (position, charge) in every output frame,
  2. discards the initial transient (automatically, or from --tmin),
  3. time-averages the number of defects N over the stationary part, with a
     statistical error that accounts for time correlations (integrated
     autocorrelation time),
  4. reports the active length  l_a = sqrt(A / <N>).

Because the solver's defect finder can split one defect into several clusters
(single-pass labelling) and misses defects in the outer cell layer, four
estimators are compared:

  all          every stationary frame, raw defect count
  neutral      only frames whose net charge equals --qexpect (charge conservation)
  merged       same-sign defects closer than --rmerge cells joined into one
  merged+neut  merged, then only neutral frames

Agreement within errors means counting artefacts are negligible. For filtered
estimators, <N> of the rejected frames is also reported: if it differs from the
kept frames, the filter biases l_a.

A is the area over which the solver actually measures topological charge.
topoChargeAndDefectsOut() (pout.c) skips the outermost cell layer, so
A = (Lx-2)(Ly-2); the naive Lx*Ly value is also written to the CSV.

The error on l_a comes from error propagation:  dl/l = (1/2) dN/N.

Usage:
    python3 active_length.py L100/act_*/ --L 100 100 [--tmin T] [--rmerge 1.5]
"""
import argparse
import json
import os
import sys

import numpy as np


# ---------------------------------------------------------------- reading
def read_frames(path):
    """Return times t and a list of (x, y, q) arrays, one per frame."""
    t, frames = [], []
    cur, nexp = None, 0
    with open(path) as f:
        for line in f:
            tok = line.split()
            if not tok:
                continue
            try:
                vals = [float(x) for x in tok]
            except ValueError:
                continue                        # banner / column headers
            if len(vals) == 2:                  # frame header: t numDefects
                if cur is not None:
                    frames.append(np.array(cur).reshape(-1, 3))
                t.append(vals[0])
                nexp = int(vals[1])
                cur = []
            elif len(vals) == 4 and cur is not None:   # QX QY charge angle
                cur.append(vals[:3])
    if cur is not None:
        if len(cur) == nexp:
            frames.append(np.array(cur).reshape(-1, 3))
        else:                                   # truncated last frame
            t = t[:-1]
    return np.asarray(t), frames


# ---------------------------------------------------------------- merging
def merge_count(fr, L, rmerge):
    """Number of defects after joining same-sign defects within rmerge
    (minimum-image distance), and the resulting net charge."""
    n = len(fr)
    if n == 0:
        return 0, 0.0
    x, y, s = fr[:, 0], fr[:, 1], np.sign(fr[:, 2])
    parent = np.arange(n)

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    dx = x[:, None] - x[None, :]
    dy = y[:, None] - y[None, :]
    dx -= L[0] * np.round(dx / L[0])
    dy -= L[1] * np.round(dy / L[1])
    close = (dx * dx + dy * dy <= rmerge * rmerge) & (s[:, None] == s[None, :])
    for i, j in zip(*np.nonzero(np.triu(close, 1))):
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[ri] = rj
    roots = {find(i) for i in range(n)}
    q = 0.5 * sum(s[r] for r in roots)
    return len(roots), q


# ---------------------------------------------------------------- statistics
def autocorr(x):
    """Normalised autocorrelation function via FFT."""
    x = np.asarray(x, float) - np.mean(x)
    n = len(x)
    if n < 2 or np.var(x) == 0:
        return np.ones(1)
    f = np.fft.rfft(x, n=2 * n)
    acf = np.fft.irfft(f * np.conj(f))[:n]
    return acf / acf[0]


def integrated_time(x, c=5.0):
    """Integrated autocorrelation time tau_int (in samples), Sokal windowing.
    Statistical inefficiency g = 2 tau_int."""
    rho = autocorr(x)
    if len(rho) == 1:
        return 0.5
    taus = np.cumsum(rho) - 0.5          # tau(M) = 1/2 + sum_{k=1..M} rho_k
    for M in range(1, len(rho)):
        if M >= c * taus[M]:
            return max(taus[M], 0.5)
    return max(taus[-1], 0.5)


def detect_equilibration(x, nskip=None):
    """Start index t0 maximising the number of effectively independent
    samples of x[t0:] (Chodera, JCTC 12, 1799 (2016))."""
    n = len(x)
    if nskip is None:
        nskip = max(1, n // 200)
    best = (0, -1.0)
    for t0 in range(0, n - 20, nskip):
        neff = (n - t0) / (2.0 * integrated_time(x[t0:]))
        if neff > best[1]:
            best = (t0, neff)
    return best[0]


def block_error(x, nblocks=10):
    """Standard error from nblocks contiguous blocks (cross-check)."""
    m = len(x) // nblocks
    if m < 1:
        return np.nan
    b = np.asarray(x[: m * nblocks], float).reshape(nblocks, m).mean(axis=1)
    return b.std(ddof=1) / np.sqrt(nblocks)


def series_stats(x, dt_out):
    """Mean, SEM (correlation-corrected), block SEM, tau (time units), n_eff.
    For filtered series the kept frames are treated as consecutive, which
    slightly underestimates tau when many frames are dropped."""
    x = np.asarray(x, float)
    if len(x) < 20:
        return None
    tau = integrated_time(x)
    neff = len(x) / (2.0 * tau)
    std = x.std(ddof=1)
    return {"N_mean": x.mean(), "N_std": std, "N_sem": std / np.sqrt(neff),
            "N_sem_block": block_error(x), "tau_int_time": tau * dt_out,
            "n_eff": neff}


# ---------------------------------------------------------------- analysis
def read_activity(rundir):
    try:
        with open(os.path.join(rundir, "input.json")) as f:
            return float(json.load(f)["species"][0]["act"])
    except Exception:
        return np.nan


def analyse(rundir, L, tmin, rmerge, qexpect):
    run = rundir.rstrip("/")
    act = read_activity(rundir)
    t, frames = read_frames(os.path.join(rundir, "defects.dat"))
    base = {"run": run, "act": act, "frames": len(frames)}
    if len(frames) < 50:
        return [dict(base, estimator="-", status="too few frames")]

    N = np.array([len(fr) for fr in frames])
    Q = np.array([0.5 * np.sign(fr[:, 2]).sum() if len(fr) else 0.0 for fr in frames])
    NQm = np.array([merge_count(fr, L, rmerge) for fr in frames])
    Nm, Qm = NQm[:, 0], NQm[:, 1]

    i0 = detect_equilibration(N) if tmin is None else int(np.searchsorted(t, tmin))
    sl = slice(i0, None)
    dt_out = np.median(np.diff(t))
    A_eff = (L[0] - 2) * (L[1] - 2)
    A_box = L[0] * L[1]

    neutral = np.abs(Q[sl] - qexpect) < 0.25
    neutral_m = np.abs(Qm[sl] - qexpect) < 0.25
    ests = {
        "all":         (N[sl], np.ones_like(neutral)),
        "neutral":     (N[sl], neutral),
        "merged":      (Nm[sl], np.ones_like(neutral)),
        "merged+neut": (Nm[sl], neutral_m),
    }

    rows = []
    for name, (x, keep) in ests.items():
        r = dict(base, estimator=name, t_start=t[i0], t_end=t[-1],
                 frames_used=int(keep.sum()), frac_kept=keep.mean(),
                 net_charge_mean=(Qm if "merged" in name else Q)[sl].mean(),
                 N_rejected_mean=x[~keep].mean() if (~keep).any() else np.nan)
        st = series_stats(x[keep], dt_out)
        if st is None:
            rows.append(dict(r, status="too few frames kept"))
            continue
        r.update(st)
        m, e = st["N_mean"], st["N_sem"]
        r["l_active"] = np.sqrt(A_eff / m)
        r["l_active_err"] = 0.5 * r["l_active"] * e / m
        r["l_active_fullbox"] = np.sqrt(A_box / m)
        h = len(x[keep]) // 2
        xk = x[keep]
        r["drift_sigma"] = (xk[h:].mean() - xk[:h].mean()) / (np.sqrt(2) * e) if e > 0 else np.nan
        r["status"] = "ok"
        rows.append(r)
    return rows


COLS = ["act", "estimator", "status", "frames", "t_start", "t_end", "frames_used",
        "frac_kept", "N_mean", "N_std", "N_sem", "N_sem_block", "N_rejected_mean",
        "tau_int_time", "n_eff", "net_charge_mean", "l_active", "l_active_err",
        "l_active_fullbox", "drift_sigma", "run"]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+", help="run directories containing defects.dat")
    ap.add_argument("--L", nargs=2, type=int, default=[100, 100], help="box size Lx Ly")
    ap.add_argument("--tmin", type=float, default=None,
                    help="discard t < tmin (MPCD time); default: automatic")
    ap.add_argument("--rmerge", type=float, default=1.5,
                    help="merge same-sign defects closer than this (cells)")
    ap.add_argument("--qexpect", type=float, default=0.0,
                    help="expected net charge (0 for PBC; set by BC otherwise)")
    ap.add_argument("--out", default="active_length.csv", help="output CSV")
    args = ap.parse_args()

    rows = [r for d in args.runs
            for r in analyse(d, args.L, args.tmin, args.rmerge, args.qexpect)]
    rows.sort(key=lambda r: r["act"])   # stable: keeps estimator order

    with open(args.out, "w") as f:
        f.write(",".join(COLS) + "\n")
        for r in rows:
            f.write(",".join(str(r.get(c, "")) for c in COLS) + "\n")

    print(f"{'act':>6} {'estimator':<12} {'t_start':>7} {'kept':>5} {'<N>':>7} {'±':>5} "
          f"{'blk±':>5} {'<N>rej':>7} {'tau':>6} {'n_eff':>5} {'l_a':>6} {'±':>5} {'drift':>6}")
    for r in rows:
        if r["status"] != "ok":
            print(f"{r['act']:6.3f} {r['estimator']:<12} {r['status']} ({r['frames']} frames)")
            continue
        print(f"{r['act']:6.3f} {r['estimator']:<12} {r['t_start']:7.0f} {r['frac_kept']:5.2f} "
              f"{r['N_mean']:7.2f} {r['N_sem']:5.2f} {r['N_sem_block']:5.2f} "
              f"{r['N_rejected_mean']:7.2f} {r['tau_int_time']:6.1f} {r['n_eff']:5.0f} "
              f"{r['l_active']:6.2f} {r['l_active_err']:5.2f} {r['drift_sigma']:6.2f}")
    print(f"\nwritten {args.out}\n"
          "kept = fraction of stationary frames used; <N>rej = mean N of rejected frames\n"
          "tau in MPCD time; drift = 2nd-half minus 1st-half mean in sigma (|drift|>2: not stationary)")


if __name__ == "__main__":
    sys.exit(main())
