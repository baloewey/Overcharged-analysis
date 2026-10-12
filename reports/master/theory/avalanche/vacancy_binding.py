#!/usr/bin/env python3
"""Test 4 of the cascade question: does a vacancy weaken the binding of its neighbours? (Frank energy, anchored disk, ejection/tilted_structures.energy)
Circle R = 31.83, k = 6, active structure: ring of 2k-2 = 10 +1/2 at the depth-law radius (structB, registry 0, natural cores) + bulk pair at +-0.45 R (y axis).
Pull one ring defect radially inward by delta (others fixed, director relaxed outside rigid cores a = 1.5) and record E(delta) - E(0) (kBT); binding = cost at delta = 6.
  F       : full ring, pull site 0;
  V_gap   : site 1 (36 deg) removed and its +1/2 put at the centre (charge conserved), pull site 0 (neighbour of the vacancy) and site 5 (far);
  V_even  : 9 evenly spaced ring + pair + centre +1/2 (ring re-spaced after the escape), pull site 0.
One-constant (32.5) and anisotropic (42, 23.1). Writes vacancy_binding.txt."""
import os, sys, time
import numpy as np
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(os.path.dirname(H), "ejection"))
from tilted_structures import energy  # noqa
R = 31.830989; N = 10; rs = R * ((N - 1 + 4) / (3 * N - 1 + 4)) ** (1 / (2 * N)) - 0.61; pair = [0.45j * R, -0.45j * R]
deltas = [0, 1, 2, 3, 4, 6, 8, 10]; L = []
def say(x): print(x, flush=True); L.append(x)
def ring(n, skip=()): return [rs * np.exp(2j * np.pi * m / n) for m in range(n) if m not in skip]
def profile(zs, j, K1, K3):
    out = []
    for d in deltas:
        z = list(zs); z[j] = z[j] * (rs - d) / rs; out.append(energy(z, [0.5] * len(z), [0.0] * len(z), R, K1, K3))
    return np.array(out) - out[0]
say(f"ring radius {rs:.2f} (depth {R - rs:.2f}); deltas {deltas}")
for lab, K1, K3 in (("one-constant", 32.5, 32.5), ("anisotropic", 42.0, 23.1)):
    t0 = time.time()
    F = ring(N) + pair; Vg = ring(N, skip=(1,)) + pair + [0j]; Ve = ring(N - 1) + pair + [0j]
    cases = {"F site 0 (full ring)": (F, 0), "V_gap site 0 (next to vacancy)": (Vg, 0), "V_gap far site (5)": (Vg, 4), "V_even site 0 (9 re-spaced)": (Ve, 0)}
    res = {}
    for name, (zs, j) in cases.items():
        res[name] = p = profile(zs, j, K1, K3); say(f"{lab:12s} {name:32s} E(delta)-E(0): " + " ".join(f"{d}:{e:+.1f}" for d, e in zip(deltas, p)))
    b = {n: p[deltas.index(6)] for n, p in res.items()}
    say(f"{lab:12s} binding (delta 6): full {b['F site 0 (full ring)']:.1f}, vacancy neighbour {b['V_gap site 0 (next to vacancy)']:.1f}, far {b['V_gap far site (5)']:.1f}, "
        f"re-spaced {b['V_even site 0 (9 re-spaced)']:.1f} kBT  ({time.time() - t0:.0f} s)")
open(os.path.join(H, "vacancy_binding.txt"), "w").write("\n".join(L) + "\n")
# robustness: where the escaped +1/2 sits (centre / just inside the vacancy / opposite side), and a ring without the extra bulk +1/2 (charge not conserved; pure vacancy)
L2 = []
for lab, K1, K3 in (("one-constant", 32.5, 32.5), ("anisotropic", 42.0, 23.1)):
    for where, ze in (("centre", [0j]), ("near vacancy 0.6R@36deg", [0.6 * R * np.exp(1j * np.pi / 5)]), ("opposite 0.6R@216deg", [0.6 * R * np.exp(1j * 6 * np.pi / 5)]), ("removed", [])):
        Vg = ring(N, skip=(1,)) + pair + ze
        bn = profile(Vg, 0, K1, K3)[deltas.index(6)]; bf = profile(Vg, 4, K1, K3)[deltas.index(6)]; bo = profile(Vg, 1, K1, K3)[deltas.index(6)]
        s_ = f"{lab:12s} escaped +1/2 {where:26s}: binding next to vacancy {bn:.1f} (site 0), {bo:.1f} (site 2), far {bf:.1f} kBT"; print(s_, flush=True); L2.append(s_)
open(os.path.join(H, "vacancy_binding.txt"), "a").write("\n" + "\n".join(L2) + "\n")
