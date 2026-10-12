#!/usr/bin/env python3
"""Per-site misfit of the tilted wall layer in ellipses (user's request, 2026-10-08).
eMode-2 anchoring: director angle k 2 pi s/P (s arclength from (A, 0), ccw); relative to the outward normal phi_n(s): psi_rel(s) = k 2 pi s/P - phi_n(s).
Homeotropic points: psi_rel = 0 (mod pi). Local lock (generalised): tail tilt of a wall +1/2 at s = 2 psi_rel(s) + const; local lock slope d tilt/ds = 2 psi_rel'(s)
= 2 (k 2 pi/P - kappa(s)).
Measured per wall +1/2 (depth < 4), W3 ellipses (2 seeds) and the circle reference: arclength s, depth, tail tilt tau (tail - normal at the nearest wall point),
predicted lock tilt 2 psi_rel(s) (mod 2 pi). Per region (tip |x| > A/2, middle, flank |x| < A/4):
  (1) lock: circular correlation and fitted offset of tau vs 2 psi_rel(s) (lock holds if tau - 2 psi_rel(s) is concentrated, resultant length near 1);
  (2) displacement from the nearest homeotropic point, |ds| (cells), and |tau| (deg); (3) misfit: m = (tau_flow - tau)/2, the rotation a core would need to
      reach the common tilt tau_flow demanded by the flow (taken as the circle value, 63 deg, with the local sign of tau); its rms per region;
  (4) predicted binding loss from the forced-rotation cost: dU = (K/2) [2 pi/arccosh(d/a)] m^2 (K = 23 kBT, a = 1.5, d = depth + 0.61), in kBT, per region.
Compared with the measured well depths (measured_binding.txt). Writes site_misfit.txt."""
import json, os, sys
import numpy as np
ROOT = "/home/baloewe/Physics/Overcharged/runs"; sys.path.insert(0, ROOT + "/cavity")
from edge_criteria import read_with_angle  # noqa
H = os.path.dirname(os.path.abspath(__file__)); L = []
def say(s): print(s, flush=True); L.append(s)
K, a = 23.0, 1.5; TFLOW = np.radians(63.0)
def wall(A, B, k, n=20000):
    t = np.linspace(0, 2 * np.pi, n, endpoint=False); W = A * np.cos(t) + 1j * B * np.sin(t)
    seg = np.abs(np.roll(W, -1) - W); s = np.concatenate([[0], np.cumsum(seg)[:-1]]); P = seg.sum()
    Nn = (B * np.cos(t) + 1j * A * np.sin(t)); Nn /= np.abs(Nn); phin = np.unwrap(np.angle(Nn))
    psi = k * 2 * np.pi * s / P - phin                                         # anchoring relative to the normal
    return W, s, P, Nn, psi
for lab, runs in (("circle", [f"{ROOT}/cavity/Rscan_R31.83_seeds/{x}/k06" for x in ("s02", "s03")]),
                  ("eps 1.2", [f"{ROOT}/cassini/W3_ellipse/eps1.2_k06/{x}" for x in ("s01", "s02")]),
                  ("eps 1.5", [f"{ROOT}/cassini/W3_ellipse/eps1.5_k06/{x}" for x in ("s01", "s02")]),
                  ("eps 2.0", [f"{ROOT}/cassini/W3_ellipse/eps2.0_k06/{x}" for x in ("s01", "s02")])):
    rec = {"tip": [], "middle": [], "flank": []}
    for run in runs:
        bc = json.load(open(run + "/input.json"))["BC"][0]; k = int(bc["customAnch"])
        A = bc["R"] * bc.get("aInv", [1, 1])[0]; B = bc["R"] * bc.get("aInv", [1, 1])[1]
        W, s, P, Nn, psi = wall(A, B, k)
        t, fr = read_with_angle(run, np.array(bc["Q"][:2], float), 0.2)
        for f in fr[::3]:
            p = f[f[:, 2] > 0]
            for x, y, q, ang in p:
                z = x + 1j * y; i = np.argmin(np.abs(W - z)); d = abs(W[i] - z)
                if d >= 4: continue
                tau = np.angle(np.exp(1j * (ang - np.angle(Nn[i]))))
                lock = np.angle(np.exp(2j * psi[i]))                        # predicted tilt up to a constant
                # offset from the nearest homeotropic point along the wall: psi(s_h) = 0 mod pi
                dpsi = np.angle(np.exp(2j * psi[i])) / 2; slope = (psi[(i + 1) % len(s)] - psi[i - 1]) / (2 * (s[1] - s[0]))
                ds = -dpsi / slope
                reg = "tip" if abs(x) > A / 2 else ("flank" if abs(x) < A / 4 else "middle")
                if lab == "circle": reg = "tip" if abs(x) > A / 2 else ("flank" if abs(x) < A / 4 else "middle")
                rec[reg].append((tau, lock, ds, d, slope))
    say(f"\n== {lab}")
    for reg in ("tip", "middle", "flank"):
        r = np.array(rec[reg])
        if len(r) < 100: continue
        tau, lock, ds, d, slope = r.T
        off = np.angle(np.mean(np.exp(1j * (tau - lock)))); conc = abs(np.mean(np.exp(1j * (tau - lock))))
        sgn = np.sign(np.sin(tau)); m = (sgn * TFLOW - tau) / 2; m = np.angle(np.exp(2j * m)) / 2
        C = 2 * np.pi / np.arccosh(np.maximum((d + 0.61) / a, 1.0001)); dU = 0.5 * K * C * m ** 2
        say(f"  {reg:6s}: n {len(r):6d}; lock slope 2 psi' {2 * np.mean(slope):.3f} rad/cell; lock: tau - 2 psi_rel concentration {conc:.2f} (offset {np.degrees(off):+.0f} deg); "
            f"|ds| from homeotropic point {np.mean(np.abs(ds)):.2f} cells; |tau| {np.degrees(np.mean(np.abs(tau))):.0f} deg; misfit rms {np.degrees(np.sqrt(np.mean(m ** 2))):.0f} deg; "
            f"forced-rotation cost median {np.median(dU):.0f} kBT (mean {np.mean(dU):.0f}); depth {np.mean(d):.2f}")
open(os.path.join(H, "site_misfit.txt"), "w").write("\n".join(L) + "\n")
