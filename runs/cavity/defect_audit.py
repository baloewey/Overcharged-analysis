"""
Topological-charge audit of tracked defects in a cavity with imposed winding k.

With all tracked defects of charge ±1/2 (the tracker keeps only |q| ~ 1/2 clusters):
    N_same = defects with the sign of k,   N_opp = defects of the opposite sign
    net charge q = (N_same - N_opp)/2 * sign(k)

Per frame the audit checks
    ok_net     q == k                                 (charge conservation)
    ok_imposed N_same >= 2|k|                        (all BC-imposed defects present)
    ok_pairs   N_same - 2|k| == N_opp                 (dynamical defects = N - 2|k| is even
                                                      and half of them are opposite sign)
Given ±1/2 charges, ok_pairs <=> ok_net; they're reported separately because the
way a frame fails tells you why:
    missing_imposed   N_same < 2|k|           -> wall defects lost
    unpaired          ok_imposed but not ok_pairs -> a pair member missed / spurious / split

N_pairs = N_opp is the number of dynamically created ±1/2 pairs (physics observable).
"""
import gzip
import os

import numpy as np


def open_data(path):
    """Open a text output file, falling back to path + '.gz' (runs compressed by run_queue.sh)."""
    if os.path.exists(path):
        return open(path)
    if os.path.exists(path + ".gz"):
        return gzip.open(path + ".gz", "rt")
    raise FileNotFoundError(path + " (nor .gz)")


def audit_frames(frames, k):
    """frames: list of (n,3) arrays (x, y, charge). Returns dict of per-frame arrays."""
    sk = 1.0 if k >= 0 else -1.0
    ak = abs(k)
    n_same = np.array([(np.sign(f[:, 2]) == sk).sum() if len(f) else 0 for f in frames])
    n_opp = np.array([(np.sign(f[:, 2]) == -sk).sum() if len(f) else 0 for f in frames])
    q = 0.5 * sk * (n_same - n_opp)
    ok_net = np.abs(q - k) < 0.25
    ok_imposed = n_same >= 2 * ak
    ok_pairs = (n_same - 2 * ak) == n_opp
    return {"n_same": n_same, "n_opp": n_opp, "q": q, "ok_net": ok_net,
            "ok_imposed": ok_imposed, "ok_pairs": ok_pairs,
            "n_dyn": n_same + n_opp - 2 * ak}


def audit_summary(a):
    """Fractions over frames."""
    n = len(a["q"])
    if n == 0:
        return {}
    return {
        "frac_ok": a["ok_net"].mean(),
        "frac_missing_imposed": (~a["ok_imposed"]).mean(),
        "frac_unpaired": (a["ok_imposed"] & ~a["ok_pairs"]).mean(),
        "frac_odd_dyn": (a["n_dyn"] % 2 != 0).mean(),
        "q_mean": a["q"].mean(),
        "N_pairs_mean": a["n_opp"].mean(),
    }


def read_frames(path):
    """Read defects.dat -> (t, list of (n,3) arrays [x, y, charge]); drops a truncated last frame."""
    t, frames, cur, nexp = [], [], None, 0
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
                t.append(v[0]); nexp = int(v[1]); cur = []
            elif len(v) == 4 and cur is not None:
                cur.append(v[:3])
    if cur is not None:
        if len(cur) == nexp:
            frames.append(np.array(cur).reshape(-1, 3))
        else:
            t = t[:-1]
    return np.array(t), frames
