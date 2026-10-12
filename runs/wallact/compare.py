"""MPCD wall-band knockout (wallActOff 4) vs the paired shape-sweep runs (same input and seed, activity everywhere): second half; wall +1/2, bulk charge, tilt coherence,
circulation, P(two), and the pull-back force at 2-3 cells (restoring_swirl.samples2, 2 t.u. steps)."""
import glob, os, sys, numpy as np
sys.path.insert(0, "/home/baloewe/Physics/Overcharged/reports/master/theory/selfprop"); sys.path.insert(0, "/home/baloewe/Physics/Overcharged/reports/master/theory/overnight")
import restoring_swirl as RS
import subprocess
edges = RS.edges
def pull(run):
    a, coh = RS.samples2(run); x = a[(a[:, 0] >= 2.0) & (a[:, 0] < 3.0), 1]; return 2 * x.mean() / x.var(), coh
for d in sorted(glob.glob("*_s0*")):
    if not os.path.exists(d + "/DONE"): continue
    rel = d.replace("_", "/")
    ref = "/home/baloewe/Physics/Overcharged/runs/shapes_server/" + rel
    for lab, run in (("knockout", d), ("reference", ref)):
        F, coh = pull(run if run.startswith("/") else os.path.abspath(run))
        out = subprocess.run([sys.executable, "/home/baloewe/Physics/Overcharged/reports/master/theory/overnight/counts.py", os.path.relpath(os.path.dirname(os.path.abspath(run)), "/home/baloewe/Physics/Overcharged/runs")],
                             capture_output=True, text=True).stdout
        line = [l for l in out.splitlines() if l.startswith(os.path.basename(run))]
        sh = line[0].split("second half:")[1].strip() if line else "?"
        print(f"{rel:30s} {lab:9s}: force at 2-3 cells {F:+.2f} (negative = back to the wall), wall-layer coherence {coh:.2f} | {sh}")
