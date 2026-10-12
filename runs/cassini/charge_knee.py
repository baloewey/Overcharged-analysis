#!/usr/bin/env python3
"""
Bulk charge from the knee of Q(c) (user's idea, 2026-10-02): traditional vs less-charged layer vs no layer.

Q(c) = accumulated charge inside level c (lobe_profiles.profile; circle c = r). Smoothed slope S = dQ/d(c/b).
Knee = the local minimum of S closest to the wall, searched between the minimum of Q (the -1/2 core) and c/b = 0.985,
that is clearly flatter than both sides: S_knee <= --left x (max slope within 0.12 inwards) and S_knee <= --right x
(max slope between the knee and the wall), and Q must rise by >= 0.5 over the 0.2 in c/b inside the knee
(a plateau has to be reached; excludes noise minima in a flat interior). Q is smoothed over 0.045 in c/b first, which removes the cell-quantised steps of the layer. Q_knee = bulk charge; k - Q_knee = layer charge; no knee = no layer.
Also reports the plateau contrast (max slope towards the wall / S_knee).

Usage: charge_knee.py RUN [RUN ...] [--labels ...]
"""
import argparse, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from lobe_profiles import profile  # noqa: E402


def knee(p, left=0.7, right=0.2, win=9, minrise=0.5):
    """Smoothing window win bins (profile with 200 shells: 9 bins = 0.045 in c/b) removes the cell-quantised steps of the
    layer (~0.02 wide) but keeps a real plateau (~0.1 wide)."""
    x = p["edges"][1:]; Q = p["Q"]
    ker = np.ones(win) / win
    Qs = np.convolve(np.r_[np.full(win, Q[0]), Q, np.full(win, Q[-1])], ker, "same")[win:-win]
    S = np.gradient(Qs, x)
    i0 = int(np.argmin(Qs)); hi = np.searchsorted(x, 0.985)
    best = None
    for i in range(max(i0 + 2, 3), hi - 2):
        if not (S[i] <= S[i - 1] and S[i] <= S[i + 1]):
            continue
        jl = np.searchsorted(x, x[i] - 0.12)
        sl = S[jl:i].max() if i > jl else np.nan
        sr = S[i + 1:hi + 1].max()
        rise = Qs[i] - Qs[max(np.searchsorted(x, x[i] - 0.2), 0)]
        if S[i] <= left * sl and S[i] <= right * sr and rise >= minrise:
            best = (x[i], Qs[i], sr / max(S[i], 1e-6), S[i])
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+"); ap.add_argument("--labels", nargs="+")
    a = ap.parse_args()
    labels = a.labels or [os.path.normpath(r) for r in a.runs]
    print(f"{'run':<34} {'k':>3} {'knee c/b':>8} {'Q_knee':>6} {'layer':>6} {'contrast':>8}")
    for run, lab in zip(a.runs, labels):
        p = profile(run, 200, 0.2); kn = knee(p)
        if kn is None:
            print(f"{lab:<34} {p['k']:>3}   no knee")
        else:
            print(f"{lab:<34} {p['k']:>3} {kn[0]:8.3f} {kn[1]:6.2f} {p['k'] - kn[1]:6.2f} {kn[2]:8.1f}")


if __name__ == "__main__":
    main()
