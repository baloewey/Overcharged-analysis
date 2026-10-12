#!/usr/bin/env python3
"""Energetics of flow-induced ejection (user's mechanism), relaxed, with the proper baseline.
  E_A0      : 2k ring, natural core orientations (the passive state without flow; its tilts wind around the ring)
  E_A(tau)  : 2k ring with every core forced to tilt tau (what flow alignment demands)
  E_B(tau)  : 2k - 2 ring tilted by registry + core (tau) and a bulk pair at +-0.45 R (natural orientation)
Reported (kBT): ejection cost without flow E_B(0) - E_A0; frustration cost of a uniform tilt on the 2k ring E_A(tau) - E_A0; ejection balance under
uniform tilt E_B(tau) - E_A(tau); and the tilt above which the ejected structure is cheaper than the natural ring, E_B(tau) < E_A0 is NOT the relevant
criterion once flow imposes the tilt: the relevant comparison is E_B(tau) vs E_A(tau). Core radius a = 1.5 and 2.5 (core dependence).
One-constant (32.5) and anisotropic (42, 23.1); k = 4, 6, 8; R = 31.83. Writes ejection_energetics.txt."""
import os, sys
import numpy as np
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
from tilted_structures import energy, structA, structB, cores_for
L = []
def say(s): print(s, flush=True); L.append(s)
R = 31.830989
for k in (4, 6, 8):
    for lab, K1, K3 in (("one-constant", 32.5, 32.5), ("anisotropic", 42.0, 23.1)):
        for a in (1.5, 2.5):
            zA, qA, nA = structA(k, R)
            EA0 = energy(zA, qA, [0.0] * len(zA), R, K1, K3, a=a)
            row = []
            for td in (0, 30, 60, 90):
                tau = np.radians(td); cA, _ = cores_for(zA, qA, R, nA, tau); zB, qB, nB = structB(k, R, 0.45 * R, tau); cB, _ = cores_for(zB, qB, R, nB, tau)
                EA = energy(zA, qA, cA, R, K1, K3, a=a); EB = energy(zB, qB, cB, R, K1, K3, a=a); row.append((td, EA - EA0, EB - EA0, EB - EA))
            say(f"k={k} {lab:12s} core {a}: ejection cost without flow (natural ring -> B untilted) {row[0][2]:+.0f} kBT | under a uniform tilt tau: "
                + " ".join(f"{td}deg: ring frustration {fa:+.0f}, B-A {d:+.0f};" for td, fa, fb, d in row))
open(os.path.join(H, "ejection_energetics.txt"), "w").write("\n".join(L) + "\n")
