#!/usr/bin/env python3
"""chi scan: vortex order of each circle run (second half of the flow output): circulation C = |<u_theta>| / <|u|> per frame (r < R - 1), mean |C| and its
sign-flip count; rms speed. Checks whether a single coherent vortex survives at chi 0.25 before reading tilts. -> chiscan_circ.txt"""
import glob, gzip, json, os, numpy as np
from multiprocessing import Pool
ROOT = "/home/baloewe/Physics/Overcharged/runs/chiscan"
def one(run):
    bc = json.load(open(run + "/input.json"))["BC"][0]; R = bc["R"]; Q = np.array(bc["Q"][:2]); fr = {}
    with gzip.open(run + "/flowfield.dat.gz", "rt") as fh:
        for line in fh:
            p = line.split()
            if len(p) < 6 or not p[0][0].isdigit(): continue
            fr.setdefault(float(p[0]), []).append([float(p[1]), float(p[2]), float(p[4]), float(p[5])])
    ts = sorted(fr)[len(fr) // 2:]; C = []; V = []
    for t in ts:
        a = np.array(fr[t]); x = a[:, 0] + 0.5 - Q[0]; y = a[:, 1] + 0.5 - Q[1]; r = np.hypot(x, y); m = (r < R - 1) & (r > 1)
        ut = (x[m] * a[m, 3] - y[m] * a[m, 2]) / r[m]; sp = np.hypot(a[m, 2], a[m, 3]); C.append(ut.mean() / sp.mean()); V.append(np.sqrt((sp ** 2).mean()))
    C = np.array(C); return os.path.basename(run), np.abs(C).mean(), int((np.diff(np.sign(C)) != 0).sum()), np.mean(V)
if __name__ == "__main__":
    with Pool(14) as p: res = p.map(one, sorted(glob.glob(ROOT + "/circ_k*_s*")))
    L = ["run | mean |C| | sign flips | rms speed"] + [f"{r} | {c:.2f} | {f} | {v:.4f}" for r, c, f, v in res]
    open("chiscan_circ.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
