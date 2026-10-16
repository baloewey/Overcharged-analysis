#!/usr/bin/env python3
"""Director rotation profile across the wall band (user, 2026-10-16), MPCD circles at R = 31.83, k 3-20, second half, director every 100 t.u.
Per depth bin (R - r, 1-cell bins, 0-14 cells): circulation-signed complex order Z = <S exp(2i(theta - phi - pi/2))> (all cells of the ring) ->
alpha = arg(Z)/2 = director angle from the local TANGENT (flow direction), positive with the circulation; |Z|/<S> = how uniform the ring is (the wall winding
(k-1)phi averages to 0 at the wall and is screened by the layer). Zero-mode picture: the far field sits at the Leslie angle from the flow, theta_L = 30 deg
(lambda 2), and the layer tilt is tau = 2 x (far-field rotation); a soft band would absorb part of the rotation. Also the depth where |Z|/<S> reaches 0.5
(screening depth of the wall winding). -> band_profile.txt, band_profile.npz"""
import os, sys, glob, json, numpy as np
from multiprocessing import Pool
ROOT = "/home/baloewe/Physics/Overcharged"; H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
import bulk_spiral as BS  # noqa
EDG = np.arange(0, 15.01, 1.0)


def one(run):
    bc = json.load(open(run + "/input.json"))["BC"][0]; R = bc["R"]; k = bc["customAnch"]; Q = np.array(bc["Q"][:2], float)
    D = BS.frames(run + "/directorfield.dat.gz", 8); Fl = BS.frames(run + "/flowfield.dat.gz", 7)
    td = np.unique(D[:, 0]); td = td[td >= td[0] + 0.5 * (td[-1] - td[0])]; Z = np.zeros(len(EDG) - 1, complex); Sm = np.zeros(len(EDG) - 1); n = 0
    for tt in td:
        d = D[D[:, 0] == tt]; f = Fl[np.abs(Fl[:, 0] - tt) < 1e-6]
        if not len(f): continue
        xf = f[:, 1] + 0.5 - Q[0]; yf = f[:, 2] + 0.5 - Q[1]; rf = np.hypot(xf, yf); circ = np.sign(((xf * f[:, 5] - yf * f[:, 4]) / np.maximum(rf, 1e-9))[rf < R - 1].mean())
        if circ == 0: continue
        x = d[:, 1] + 0.5 - Q[0]; y = d[:, 2] + 0.5 - Q[1]; r = np.hypot(x, y); ph = np.arctan2(y, x); dep = R - r; m = (dep >= 0) & (dep < EDG[-1]) & (d[:, 7] > 0)
        z = d[m, 7] * np.exp(2j * circ * (np.arctan2(d[m, 5], d[m, 4]) - ph[m] - np.pi / 2))
        b = np.digitize(dep[m], EDG) - 1
        for i in range(len(EDG) - 1):
            q = b == i
            if q.any(): Z[i] += z[q].mean(); Sm[i] += d[m, 7][q].mean()
        n += 1
    return dict(k=k, a=2 * np.pi * R / (2 * k - 2), Z=Z / n, S=Sm / n)


if __name__ == "__main__":
    runs = sorted(os.path.dirname(p) for p in glob.glob(ROOT + "/runs/cavity/server_main/main_act0.07/rho*/k*/s[01][0-9]/directorfield.dat.gz")
                  if abs(json.load(open(os.path.dirname(p) + "/input.json"))["BC"][0]["R"] - 31.83) < 0.05)
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 6) as pool: res = pool.map(one, runs)
    T = {}
    for l in open(os.path.join(H, "bulk_spiral.txt")):
        p = [c.strip() for c in l.split("|")]
        if len(p) > 4 and p[0].isdigit(): T[int(p[0])] = float(p[3])
    ks = sorted(set(x["k"] for x in res)); mid = 0.5 * (EDG[1:] + EDG[:-1])
    L = ["DIRECTOR ACROSS THE WALL BAND (R 31.83; alpha = angle from the tangent, + with the circulation; uniformity |Z|/<S>); Leslie theta_L = 30 deg (lambda 2)",
         "k | spacing | tilt tau | tau/2 | alpha (uniformity) at depth 0.5, 1.5, 2.5, 3.5, 4.5, 6.5, 8.5, 10.5, 12.5 | depth where uniformity > 0.5 | max alpha (depth)"]
    out = {}
    for kk in ks:
        S_ = [x for x in res if x["k"] == kk]; Z = np.mean([x["Z"] for x in S_], 0); S = np.mean([x["S"] for x in S_], 0); u = np.abs(Z) / S; al = np.degrees(np.angle(Z) / 2)
        out[kk] = (al, u); dsc = mid[np.argmax(u > 0.5)] if (u > 0.5).any() else np.nan; j = np.argmax(np.where(u > 0.3, al, -99))
        L.append(f"{kk:2d} | {S_[0]['a']:5.1f} | {T.get(kk, np.nan):5.1f} | {T.get(kk, np.nan)/2:5.1f} | " + "  ".join(f"{al[i]:+5.1f} ({u[i]:.2f})" for i in (0, 1, 2, 3, 4, 6, 8, 10, 12)) + f" | {dsc:.1f} | {al[j]:+.1f} ({mid[j]:.1f})")
    np.savez(os.path.join(H, "band_profile.npz"), ks=np.array(ks), mid=mid, alpha=np.array([out[k][0] for k in ks]), unif=np.array([out[k][1] for k in ks]))
    print("\n".join(L)); open(os.path.join(H, "band_profile.txt"), "w").write("\n".join(L) + "\n")
