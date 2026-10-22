#!/usr/bin/env python3
"""Gauss decomposition of the drag energy (analytic, one constant K 44, point charges; geometry as drag_energy.py: circle R 31.83, ring of N = 2k-2 wall +1/2 at the
depth-law radius, central +1). Probe = wall defect 0 pulled radially inward by delta. The other N-1 wall charges are written exactly as
  (a) a SMEARED ring of all N charges (uniform line, total N/2) : by the 2D shell theorem it exerts no force inside it (with its wall image);
  (b) the HOLE: charge -1/2 at the probe's original site (removes the probe's own share), probe-hole distance regularised with the core radius 1.5;
  (c) the CORRUGATION: discrete N-1 charges minus (a) minus (b) = the periodic part, which decays away from the layer over ~ spacing/2 pi.
Plus the probe's SELF-IMAGE (wall repulsion) and the central +1. E(delta) - E(0) for each part, kBT. -> gauss_split.txt"""
import numpy as np
R = 31.830989; K = 44.0; KS = [4, 6, 9, 12, 14]; DEL = [0, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 12]
A = 1.5                                                                                  # core radius: hole-probe distance regularised as sqrt(d^2 + A^2)
def hole(zp, z0): return 0.5 * (0.5 * np.log(abs(zp - z0) ** 2 + A ** 2) + np.log(abs(R ** 2 - zp * np.conj(z0))))
def pair(zp, zo, qo):     # probe (q 1/2) with charges zo, qo and their images, units of pi K
    zo = np.asarray(zo); qo = np.asarray(qo); return -np.sum(qo * (np.log(np.abs(zp - zo)) + np.log(np.abs(R ** 2 - zp * np.conj(zo)))))
L = ["GAUSS SPLIT of the drag energy (kBT; analytic one constant K 44): columns delta = " + str(DEL)]
for k in KS:
    N = 2 * k - 2; rs = R * ((N - 1 + 4) / (3 * N - 1 + 4)) ** (1 / (2 * N)) - 0.61; z0 = rs + 0j; M = 4000
    zr = rs * np.exp(2j * np.pi * np.arange(N) / N); sm = rs * np.exp(2j * np.pi * (np.arange(M) + 0.5) / M)
    P = {"discrete N-1 (total wall part)": [], "smeared ring of N": [], "hole -1/2": [], "corrugation": [], "self-image": [], "central +1": []}
    for dl in DEL:
        zp = z0 - dl; d = pair(zp, zr[1:], [0.5] * (N - 1)); s = pair(zp, sm, [0.5 * N / M] * M); h = hole(zp, z0)
        P["discrete N-1 (total wall part)"].append(d); P["smeared ring of N"].append(s); P["hole -1/2"].append(h); P["corrugation"].append(d - s - h)
        P["self-image"].append(-0.25 * np.log(R ** 2 - abs(zp) ** 2)); P["central +1"].append(pair(zp, [0j], [1.0]) + 0.5 * np.log(R ** 2))
    L.append(f"\nk {k} (N {N}, rho {2 * k / (2 * np.pi * R):.2f}, spacing {2 * np.pi * rs / N:.1f}, wall distance of the ring {R - rs:.2f})")
    tot = 0
    for nm, v in P.items():
        v = np.pi * K * (np.array(v) - v[0]); L.append(f"  {nm:32s}: " + " ".join(f"{x:7.1f}" for x in v))
        if nm not in ("smeared ring of N", "hole -1/2", "corrugation"): tot = tot + v
    L.append(f"  {'TOTAL':32s}: " + " ".join(f"{x:7.1f}" for x in tot))
    L.append(f"  line-charge slope pi^2 K rho / 2 = {np.pi ** 2 * K * (2 * k / (2 * np.pi * R)) / 2:.1f} kBT/cell (constant force of an infinite uniform line of +1/2 at density rho)")
open("gauss_split.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
