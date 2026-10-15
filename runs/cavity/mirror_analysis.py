#!/usr/bin/env python3
"""H5: mirror runs (customAnch = -k). Per run: peak of Q(r) before the wall plunge (r/R < 0.97) and its position, width of
the region within 0.1 of the peak, mean number of -1/2 and +1/2 within 3 cells of the wall per frame.
Bulk +1 with a wall layer of charge -k-1 would need 2|k|+2 wall -1/2. Writes analysis_weekend/mirror.csv."""
import csv, glob, json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "cassini"))
from lobe_profiles import profile  # noqa: E402
from defect_spatial import load    # noqa: E402
out = os.path.join(HERE, "analysis_weekend", "mirror.csv")
w = csv.writer(open(out, "w", newline="")); w.writerow(["run", "k", "R", "rho", "Qpeak", "r_peak", "width_peak", "wall_minus", "wall_plus", "bulk_minus", "bulk_plus"])
for d in sorted(glob.glob(os.path.join(HERE, "server_main", "main_act0.07_neg", "rho*", "km*", "s[01][0-9]", "DONE"))):
    run = os.path.dirname(d)
    try:
        bc = json.load(open(os.path.join(run, "input.json")))["BC"][0]; R = bc["R"]; k = bc["customAnch"]
        p = profile(run, 200, 0.2); x = p["edges"][1:]; Q = p["Q"]
        m = x < 0.97; i = int(np.argmax(np.where(m, Q, -1e9))); qp = Q[i]
        band = np.where(m & (Q > qp - 0.1))[0]; width = (x[band].max() - x[band].min()) if len(band) else 0
        dd = load(run, 0.2); wm = wp = bm = bp = 0
        for f in dd["frames"]:
            s = R - np.hypot(f[:, 0], f[:, 1])
            wm += np.sum((s < 3) & (f[:, 2] < 0)); wp += np.sum((s < 3) & (f[:, 2] > 0))
            bm += np.sum((s >= 3) & (f[:, 2] < 0)); bp += np.sum((s >= 3) & (f[:, 2] > 0))
        n = len(dd["frames"])
        w.writerow([os.path.relpath(run, os.path.join(HERE, "server_main")), k, R, abs(k) / (np.pi * R), qp, x[i], width, wm / n, wp / n, bm / n, bp / n])
    except Exception as e:
        print(run, "FAILED", e, flush=True)
print("written", out)
