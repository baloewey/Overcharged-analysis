"""MPCD: director angle relative to the wall tangent, alpha(depth) = arg<S e^{2i(theta - phi - 90)}>/2, and S(depth), azimuthally averaged over a window
with steady imposed flow (before the flip, or the second half). Also the wall-defect tilt in the same window (defects.dat, tail angle - normal, coherent mean).
Prediction: tilt = arccos(1/lam_eff), lam_eff = lam (1 + S4)/(2 S2) (2D moment of Jeffery's equation), S4 from a Maier-Saupe (von Mises) distribution with <cos 2psi> = S2."""
import gzip, json, os, sys
import numpy as np
from scipy.special import i0e, i1e
from scipy.optimize import brentq
from scipy.spatial import cKDTree
ROOT = "/home/baloewe/Physics/Overcharged/runs"; sys.path.insert(0, ROOT + "/cavity")
from edge_criteria import read_with_angle  # noqa
def S4_of(S2):
    k = brentq(lambda x: i1e(x) / i0e(x) - S2, 1e-6, 500); return 1 - 2 * S2 / k
def lam_eff(lam, S2): return lam * (1 + S4_of(S2)) / (2 * S2)
def analyse(rel, t0, t1, stride=10):
    run = os.path.join(ROOT, rel); J = json.load(open(run + "/input.json")); bc = J["BC"][0]; Qc = bc["Q"]; R = bc["R"]; lam = J["species"][0].get("tumble", 2)
    acc = {}; nfr = 0; times = set()
    fn = run + "/directorfield.dat"; op = open
    if not os.path.exists(fn): fn += ".gz"; op = gzip.open
    with op(fn, "rt") as f:
        for line in f:
            p = line.split()
            if len(p) != 8 or not p[0][0].isdigit(): continue
            t = float(p[0])
            if t < t0 or t > t1: continue
            times.add(t)
            if (len(times) - 1) % stride: continue
            S = float(p[7])
            if S <= 0: continue
            x = float(p[1]) + 0.5 - Qc[0]; y = float(p[2]) + 0.5 - Qc[1]; r = np.hypot(x, y)
            if r >= R: continue
            th = np.arctan2(float(p[5]), float(p[4])); al2 = 2 * (th - np.arctan2(y, x) - np.pi / 2); d = int((R - r) / 0.5)
            a = acc.setdefault(d, [0j, 0.0, 0]); a[0] += S * np.exp(1j * al2); a[1] += S; a[2] += 1
    prof = [((d + 0.5) * 0.5, np.degrees(np.angle(a[0])) / 2, a[1] / a[2], abs(a[0]) / a[2]) for d, a in sorted(acc.items())]
    t, fr = read_with_angle(run, np.array(Qc[:2], float), 0.0); zs = []
    for tt, f in zip(t, fr):
        if t0 <= tt <= t1 and len(f):
            r = np.hypot(f[:, 0], f[:, 1]); w = (r > R - 4) & (f[:, 2] > 0)
            if w.any(): zs.append(np.mean(np.exp(1j * (f[w, 3] - np.arctan2(f[w, 1], f[w, 0])))))
    z = np.mean(zs); return lam, prof, np.degrees(np.angle(z)), abs(z)
if __name__ == "__main__":
    runs = [("swirl_tumble/passive_lam0.5_A0.01", 1000, 1500), ("swirl_tumble/passive_lam1.2_A0.01", 1000, 1500), ("swirl_threshold/passive_A0.01_s2", 1500, 3000),
            ("swirl_threshold/passive_A0.01_s3", 1500, 3000), ("swirl_tumble/passive_lam4.0_A0.01", 1000, 1500)]
    L = []
    for rel, t0, t1 in runs:
        lam, prof, tilt, coh = analyse(rel, t0, t1)
        L.append(f"{rel}: lam {lam}, measured wall tilt {tilt:+.1f} (coherence {coh:.2f})")
        for dep, al, S, o in prof:
            if dep in (0.25, 0.75, 1.25, 1.75, 2.25, 2.75, 3.75, 4.75, 6.25, 8.25, 12.25, 16.25, 20.25):
                le = lam_eff(lam, min(S, 0.999)) if S > 0.05 else np.nan
                L.append(f"   depth {dep:5.2f}: alpha {al:+6.1f}  2alpha {2 * al:+6.1f}  S {S:.3f}  order {o:.2f}  lam_eff {le:.2f}  arccos(1/lam_eff) "
                         f"{np.degrees(np.arccos(1 / le)) if le > 1 else float('nan'):.1f}")
        print("\n".join(L[-14:]), flush=True)
    open("mpcd_profile.txt", "w").write("\n".join(L) + "\n")
