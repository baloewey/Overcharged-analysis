#!/usr/bin/env python3
"""Are the vortices of the long winding channel pinned, possibly with sign flips? (user, 2026-10-15: the movie suggests pinned vortices; in a vortex lattice the
vorticity can flip). Sign-blind tests on the mid-band vorticity omega(x, t) (vortex_row.kymo):
  (1) autocorrelation in time at FIXED position (no drift allowed) of omega (signed) and of |omega| (sign-blind), averaged over x;
  (2) time-averaged |omega|(x) profile: pinned vortex positions give a modulation along x; its contrast compared with the same quantity from the two halves of the run
      (reproducible peaks = pinning) and with a surrogate where each frame is shifted by a random amount (no pinning);
  (3) at fixed x, the sign of omega: mean duration of a sign (flip time) and whether flips are periodic (peak in the spectrum of omega(x, t) in time).
-> vortex_pinning.txt, fig_vortex_pinning.png (+pdf)"""
import os, sys, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
from vortex_row import kymo  # noqa
R_ = "/home/baloewe/Physics/Overcharged/runs/channel_long"; rng = np.random.default_rng(4)
L_ = ["VORTEX PINNING TEST (mid-band vorticity, second half, every 10 t.u.)"]
fig, ax = plt.subplots(1, 3, figsize=(16, 4.4))
for name, col, lab in (("W40_m24_L768", "#c2452d", "winding walls"), ("W40_m0_L768", "0.4", "homeotropic walls")):
    ts, K, J, VX = kymo(os.path.join(R_, name)); dt = np.median(np.diff(ts)); n, L = K.shape
    def acf(X):
        Z = (X - X.mean(0)) / (X.std(0) + 1e-12); return np.array([np.mean(Z[:n - l] * Z[l:]) for l in range(0, 40)])
    cs, ca = acf(K), acf(np.abs(K))
    prof = np.abs(K).mean(0); h1, h2 = np.abs(K[: n // 2]).mean(0), np.abs(K[n // 2:]).mean(0)
    rep = np.corrcoef(h1 - h1.mean(), h2 - h2.mean())[0, 1]
    sur = np.array([np.roll(np.abs(K[i]), rng.integers(L)) for i in range(n)]).mean(0)
    contrast, contrast_sur = prof.std() / prof.mean(), sur.std() / sur.mean()
    sgn = np.sign(K); flips = np.sum(sgn[1:] != sgn[:-1], 0); tflip = (n - 1) * dt / np.maximum(flips, 1)
    spec = np.mean(np.abs(np.fft.rfft(K - K.mean(0), axis=0)) ** 2, 1); fq = np.fft.rfftfreq(n, dt); pk = np.argmax(spec[1:]) + 1
    tdec_s = np.argmax(cs < np.exp(-1)) * dt; tdec_a = np.argmax(ca < np.exp(-1)) * dt if (ca < np.exp(-1)).any() else np.inf
    L_.append(f"{name}: fixed-x autocorrelation of omega: 1/e after {tdec_s:.0f} t.u., minimum {cs.min():+.2f} at {np.argmin(cs)*dt:.0f} t.u.; of |omega|: 1/e after {tdec_a:.0f} t.u., "
              f"value at 100/200/390 t.u. {ca[10]:.2f}/{ca[20]:.2f}/{ca[39]:.2f}; time-averaged |omega|(x): contrast {contrast:.3f} (shuffled {contrast_sur:.3f}), "
              f"halves correlation {rep:+.2f}; sign at fixed x lasts {np.median(tflip):.0f} t.u. (median); temporal spectrum peak at period {1/fq[pk]:.0f} t.u. (peak/mean {spec[pk]/spec[1:].mean():.1f})")
    lag = np.arange(40) * dt; ax[0].plot(lag, cs, "-", color=col, label=f"{lab}: signed $\\omega$"); ax[0].plot(lag, ca, "--", color=col, label=f"{lab}: $|\\omega|$")
    if name == "W40_m24_L768":
        x = np.arange(L) + 0.5; ax[1].plot(x, h1, color="#c2452d", lw=1, label="first half"); ax[1].plot(x, h2, color="#1f5fa8", lw=1, label="second half"); ax[1].plot(x, sur, color="0.6", lw=1, label="shuffled positions")
    ax[2].semilogy(fq[1:], spec[1:] / spec[1:].mean(), color=col, label=lab)
ax[0].axhline(np.exp(-1), color="0.7", ls=":"); ax[0].axhline(0, color="0.85"); ax[0].set_xlabel("time lag (t.u.)"); ax[0].set_ylabel("autocorrelation at fixed position")
ax[0].set_title("(a) Memory of the pattern, no drift allowed", fontsize=10); ax[0].legend(fontsize=7.5, frameon=False)
ax[1].set_xlabel("position along the channel x (cells)"); ax[1].set_ylabel("time-averaged $|\\omega|$ (1/t.u.)"); ax[1].set_title("(b) Winding walls: are vortex positions fixed?", fontsize=10); ax[1].legend(fontsize=7.5, frameon=False)
ax[2].set_xlabel("frequency (1/t.u.)"); ax[2].set_ylabel("temporal spectrum of $\\omega$ (normalised)"); ax[2].set_title("(c) Periodic sign flips?", fontsize=10); ax[2].legend(fontsize=7.5, frameon=False)
for a_ in ax: [a_.spines[s_].set_visible(False) for s_ in ("top", "right")]
fig.tight_layout(); fig.savefig(os.path.join(H, "fig_vortex_pinning.png"), dpi=180)
print("\n".join(L_)); open(os.path.join(H, "vortex_pinning.txt"), "w").write("\n".join(L_) + "\n")
