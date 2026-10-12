#!/usr/bin/env python3
"""Rotational viscosity of the MPCD nematic from the decay of equilibrium director modes (passive periodic box, director along x; runs of frank_constant.py).
theta_q(t) as in frank_constant.py; normalised autocorrelation C_q(t) = Re<theta_q(t0+t) theta_q*(t0)> / <|theta_q|^2>, averaged over t0 and over +-q.
Decay rate Gamma_q from a fit of log C over C > 0.2 (and from the 1/e time). Mode along y (q = (0, q)): splay (K1); along x (q = (q, 0)): bend (K3).
gamma_eff = K q^2 / Gamma_q with K1 = 42.0, K3 = 23.1 (this box's measured values). Writes rot_visc.txt."""
import sys, numpy as np
from frank_constant import read_frames
K1, K3 = 42.0, 23.1
def thetas(run, t0=200.0):
    fr = read_frames(run, t0); L = int(fr[0][:, 0].max()) + 1; out = []
    for f in fr:
        th = np.zeros((L, L)); ix, iy = f[:, 0].astype(int), f[:, 1].astype(int); ang = np.arctan2(f[:, 3], f[:, 2])
        th0 = 0.5 * np.arctan2(np.mean(np.sin(2 * ang)), np.mean(np.cos(2 * ang))); d = np.angle(np.exp(2j * (ang - th0))) / 2
        th[ix, iy] = d; out.append(np.fft.fft2(th))
    return np.array(out), L
def acf(series, nlag):
    x = series - series.mean(); n = len(x); c = np.array([np.mean((x[l:] * np.conj(x[:n - l])).real) for l in range(nlag)]); return c / c[0]
L_ = []
for run in sys.argv[1:]:
    TQ, L = thetas(run); dt = 2.0; L_.append(f"{run}: L {L}, frames {len(TQ)}")
    for m in (1, 2, 3, 4, 6):
        q = 2 * np.pi * m / L
        for lab, idxs, K in (("splay (q || y)", [(0, m), (0, -m)], K1), ("bend  (q || x)", [(m, 0), (-m, 0)], K3)):
            C = np.mean([acf(TQ[:, i % L, j % L], 60) for i, j in idxs], axis=0); t = np.arange(len(C)) * dt
            sel = (C > 0.2) & (t > 0); G = -np.polyfit(t[sel], np.log(C[sel]), 1)[0] if sel.sum() >= 2 else np.nan
            te = np.interp(np.exp(-1), C[::-1], t[::-1]) if (C < np.exp(-1)).any() else np.nan
            L_.append(f"   m {m} q {q:.3f} {lab}: Gamma {G:.4f} (1/e time {te:5.1f}) -> gamma_eff = K q^2/Gamma = {K * q * q / G:6.2f}   C(2,4,8,16) {np.round(C[[1, 2, 4, 8]], 3)}")
print("\n".join(L_)); open("rot_visc.txt", "w").write("\n".join(L_) + "\n")
