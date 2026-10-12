"""Solvability (zero-mode projection) for the flat-wall tilt correction.
Static one-constant texture theta_0 (relaxed flatwall state, s = 0), rotated along its zero mode so that the far field is at theta_L:
  theta_L-texture  T(x, y) = theta_0(x - D, y) + q D,  q = pi / a_s.
Zero mode shape  m(x, y) = dT/dD = q - d_x T  (0 at the wall, q in the far field); normalised weight w = m / q.
Steady torque density (one-constant Ericksen-Leslie, simple shear v_x = s y):  g = gamma1 [ (s/2)(lam cos 2T - 1) - s y d_x T ].
Projected torque per period  P = int int w g dx dy  (far field contributes 0 because T -> theta_L there).
The outer region (y >~ a_s) is a 1D alignment boundary layer: K delta'' = gamma1 s sqrt(lam^2 - 1) delta, restoring torque per unit length K delta0 / l_E,
l_E = sqrt(K / (gamma1 s sqrt(lam^2 - 1))). Solvability: P / a_s + K delta0 / l_E = 0 -> delta0 = -P l_E / (K a_s); tilt correction = 2 delta0.
Compared with (i) the flat-wall simulations, (ii) the band ansatz P = -gamma1 s a_s d / 2. Model units gamma1 = K = 1."""
import json, numpy as np
import flatwall as F
np.seterr(all="ignore")
def static(a_s=20.0, H=30.0):
    # reuse flatwall.run internals: relax without flow, return q1, q2 and grid
    import types
    S0, h = F.S0, F.h; nx, ny = int(round(a_s / h)), int(round(H / h)) + 2
    x = (np.arange(nx) + 0.5) * h; y = (np.arange(ny) - 1) * h; X, Y = np.meshgrid(x, y, indexing="ij")
    out = {}
    def capture(q1, q2, X_, Y_, a):
        out["q"] = (q1.copy(), q2.copy()); return F.measure.__wrapped__(q1, q2, X_, Y_, a) if hasattr(F.measure, "__wrapped__") else orig(q1, q2, X_, Y_, a)
    orig = F.measure; F.measure = capture
    try: F.run(2.0, 0.0, a_s=a_s, H=H, T0=1500, T1=0)
    finally: F.measure = orig
    return out["q"], X, Y
def projected(lam, a_s=20.0, s=1.0):
    (q1, q2), X, Y = static(a_s); h = F.h; q = np.pi / a_s
    th2 = np.unwrap(np.unwrap(np.arctan2(q2, q1), axis=1), axis=0)                    # 2 theta, unwrapped
    dthx = 0.5 * (np.roll(th2, -1, 0) - np.roll(th2, 1, 0)); dthx = np.angle(np.exp(1j * dthx)) / (2 * h)   # d_x theta (wrap-safe)
    th = 0.5 * np.arctan2(q2, q1); far = np.angle(np.mean(np.exp(2j * th[:, -12:-2]))) / 2; thL = 0.5 * np.arccos(1 / lam)
    T = th + (thL - far)                                                                    # rotate along the zero mode (translation does not change the integrand)
    w = 1 - dthx / q; inn = Y > 0
    g = (s / 2) * (lam * np.cos(2 * T) - 1) - s * np.maximum(Y, 0) * dthx
    P = np.sum((w * g)[inn]) * h * h; Pv = np.sum((w * (s / 2) * (lam * np.cos(2 * T) - 1))[inn]) * h * h
    S = 2 * np.hypot(q1, q2); i = np.unravel_index(np.argmin(np.where(inn & (Y > 0.3), S, 9)), S.shape); d = Y[i]
    return P, Pv, d
if __name__ == "__main__":
    R = json.load(open("flatwall_scan.json")); L = []
    cache = {}
    for r in R:
        if r["a_s"] != 20.0 or r["H"] != 30.0 or "tilt" not in r: continue
        lam, s = r["lam"], r["s"]
        if lam not in cache: cache[lam] = projected(lam)
        P1, Pv1, d0 = cache[lam]; P, Pv = P1 * s, Pv1 * s                                  # torque is linear in s at fixed texture
        lE = np.sqrt(1 / (s * np.sqrt(lam ** 2 - 1))); corr = np.degrees(-2 * P * lE / 20.0); corr_v = np.degrees(-2 * Pv * lE / 20.0)
        meas = r["tilt"] - np.degrees(np.arccos(1 / lam)); ans = np.degrees(r["depth"] * np.sqrt(s) * (lam ** 2 - 1) ** -0.25)
        L.append(f"lam {lam:3.1f} s {s:.3f}: measured excess {meas:+5.1f} | solvability (static texture, full) {corr:+5.1f}, vorticity+strain only {corr_v:+5.1f} | band ansatz {ans:+5.1f}"
                 f"   [P/(s a_s) = {P1 / 20:+.3f}, static depth {d0:.2f}]")
    print("\n".join(L)); open("solvability.txt", "w").write("\n".join(L) + "\n")
