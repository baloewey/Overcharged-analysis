#!/usr/bin/env python3
"""
Frank constants of the 2D MPCD nematic from equilibrium director fluctuations.

Input: a passive (zeta = 0), defect-free, periodic L x L box, density 20, mean-field potential U = 10 kBT, collision operator 21,
initially aligned along x, with the director field written every 2 t.u. (directorfield.dat, one line per cell: t, cell x, cell y,
cell z, n_x, n_y, n_z, S). kBT = 1, cell size = 1.

Theory (small rotations of a director whose mean points along e_par):
    write n = (cos theta, sin theta) with theta measured from the mean director. To linear order
        splay  div n  = d_perp theta      (derivative perpendicular to the mean director)
        bend   curl n = d_par  theta      (derivative along the mean director)
    so the Frank energy F = 1/2 int [K1 (div n)^2 + K3 (curl n)^2] dA is, per Fourier mode theta_q = sum_r theta(r) e^{-i q.r},
        F = (1/2A) sum_q (K1 q_perp^2 + K3 q_par^2) |theta_q|^2            (A = L^2)
    and equipartition (each mode carries kBT/2) gives
        <|theta_q|^2> = A kBT / (K1 q_perp^2 + K3 q_par^2).
    Modes with q along the director are pure bend, modes with q across it pure splay; oblique modes mix the two.

Why the director frame matters: the box is aligned along x, but its mean director wanders slowly (up to ~16 deg). The original version
of this script used the lab axes (q_par = q_x, q_perp = q_y); with a mean tilt phi the "bend" axis then measures K3 cos^2 phi + K1 sin^2 phi.
Here each frame's wavevectors are rotated into that frame's director frame before fitting. The lab-frame fit is kept for comparison.

Fits (all modes with 0 < |q| < QMAX; frames with t >= T0):
  lab_lsq   original estimator: frame-averaged P(q) = <|theta_q|^2>, least squares of A/P = K1 q_y^2 + K3 q_x^2 + c over the modes.
  lab_ml    maximum likelihood, lab axes (see below).
  dir_ml    maximum likelihood, director-frame axes (the result we quote).
  dir_lsq   cross-check with the original estimator in the director frame: frames grouped in 2-deg bins of the mean angle, P averaged within
            each bin, least squares over (bin, mode) with the bin's mean angle used to rotate q.
Maximum likelihood: in equilibrium theta_q is a complex Gaussian, so the per-frame power P = |theta_q|^2 is exponentially distributed with
mean A/D(q), D = K1 q_perp^2 + K3 q_par^2 + c. Minimising sum over frames and modes of [ -log(D/A) + P D/A ] (convex, since D is linear in
the parameters) handles a different rotation in every frame without averaging P first. Frames are correlated in time, so this gives the
estimate but not its error; seed-to-seed and box-size spread give the error.

The constant c: an empirical offset in the inverse spectrum that absorbs short-wavelength corrections (cell-level estimation noise, lattice
effects). Note it is not literally white noise, which would add to P rather than to 1/P; it matters only at small q, where c/q^2 is
comparable to K (see K_eff below).

K_eff(q) for plotting: modes within +-ANG_WIN of the director (bend) or of its normal (splay), binned in |q|, K_eff = A kBT/(q^2 <P>).

Usage: frank_constant.py RUN [RUN ...] [--qmax 0.8] [--t0 200]
"""
import argparse, gzip, os
import numpy as np
from scipy.optimize import minimize

ANG_WIN = np.radians(10.0)   # half-width of the angular window that counts a mode as on-axis (bend or splay) for K_eff(q)
BIN_DEG = 2.0                # bin width (deg) of the mean director angle for the dir_lsq cross-check
MIN_BIN_FRAMES = 20          # dir_lsq: ignore angle bins with fewer frames (1/mean of few exponentials is biased)


def read_frames(run, t0):
    """Return a list of frames (t >= t0); each frame is an array of rows (cell x, cell y, n_x, n_y, S)."""
    f = os.path.join(run, "directorfield.dat")
    op = gzip.open(f + ".gz", "rt") if os.path.exists(f + ".gz") else open(f)
    frames, cur, tcur = [], [], None
    with op as fh:
        for line in fh:
            p = line.split()
            if len(p) != 8 or not p[0][0].isdigit(): continue      # skip the banner and header lines
            t = float(p[0])
            if tcur is not None and t != tcur:                       # a new time stamp closes the previous frame
                if tcur >= t0: frames.append(np.array(cur))
                cur = []
            tcur = t; cur.append((int(p[1]), int(p[2]), float(p[4]), float(p[5]), float(p[7])))
    if cur and tcur >= t0: frames.append(np.array(cur))
    return frames


def spectra(frames):
    """Per frame: mean director angle theta0 and the power |theta_q|^2 of the deviation field. Returns L, theta0[nf], P[nf, L, L], <S>."""
    L = int(frames[0][:, 0].max()) + 1
    P = np.empty((len(frames), L, L)); th0s = np.empty(len(frames)); S = []
    for i, f in enumerate(frames):
        ix, iy = f[:, 0].astype(int), f[:, 1].astype(int)
        # every frame must hold each of the L^2 cells exactly once: a missing cell would stay 0 in th (read as perfectly aligned)
        # and bias the spectrum, e.g. a frame cut off mid-write at the end of a run
        assert len(f) == L * L and len(np.unique(ix * L + iy)) == L * L and ix.max() < L and iy.max() < L, \
            f"incomplete frame {i}: {len(f)} rows, {len(np.unique(ix * L + iy))} distinct cells, expected {L * L}"
        ang = np.arctan2(f[:, 3], f[:, 2])                           # director angle of each cell (defined mod pi)
        # mean director of the frame from the cell-averaged Q tensor (Q ~ (cos 2a, sin 2a)), so +-n count the same
        th0 = 0.5 * np.arctan2(np.mean(np.sin(2 * ang)), np.mean(np.cos(2 * ang)))
        d = (ang - th0 + np.pi / 2) % np.pi - np.pi / 2              # deviation from the mean, wrapped into [-pi/2, pi/2)
        th = np.zeros((L, L)); th[ix, iy] = d
        P[i] = np.abs(np.fft.fft2(th)) ** 2                          # |theta_q|^2 with theta_q = sum_r theta(r) e^{-iqr}
        th0s[i] = th0; S.append(f[:, 4].mean())
    return L, th0s, P, float(np.mean(S))


def rotate(QX, QY, th0):
    """Wavevector components along (par) and across (perp) a director at angle th0."""
    c, s = np.cos(th0), np.sin(th0)
    return QX * c + QY * s, -QX * s + QY * c


def fit_ml(qpar2, qperp2, P, A):
    """Maximum-likelihood K1, K3, c for exponential P with mean A/D, D = K1 qperp2 + K3 qpar2 + c. Arrays are flat over (frame, mode)."""
    def nll(p):
        D = p[0] * qperp2 + p[1] * qpar2 + p[2]
        lam = D / A                                                  # rate of the exponential = 1 / mean
        val = np.sum(-np.log(lam) + lam * P)
        r = -1.0 / lam + P                                           # d nll / d lam
        grad = np.array([np.sum(r * qperp2), np.sum(r * qpar2), np.sum(r)]) / A
        return val, grad
    res = minimize(nll, x0=[30.0, 30.0, 0.05], jac=True, method="L-BFGS-B", bounds=[(1e-3, None), (1e-3, None), (0.0, None)])
    return res.x


def fit_lsq(qpar2, qperp2, y):
    """Least squares of y = A/<P> against K1 qperp2 + K3 qpar2 + c (the original estimator)."""
    X = np.column_stack([qperp2, qpar2, np.ones(len(y))])
    return np.linalg.lstsq(X, y, rcond=None)[0]


def analyse(run, qmax=0.8, t0=200.0):
    frames = read_frames(run, t0)
    L, th0, P, S = spectra(frames); nf = len(frames); A = L * L
    q = 2 * np.pi * np.fft.fftfreq(L)                                # FFT wavenumbers, spacing 2 pi / L
    QX, QY = np.meshgrid(q, q, indexing="ij")
    Q2 = QX ** 2 + QY ** 2
    m = (Q2 > 0) & (Q2 < qmax ** 2)                                  # modes used in every fit
    qx, qy = QX[m], QY[m]; Pm = P[:, m]                              # Pm[frame, mode]
    out = dict(L=L, frames=nf, S=S, th0_mean=np.degrees(th0.mean()), th0_sd=np.degrees(th0.std()), th0_max=np.degrees(np.abs(th0).max()))

    # 1. original: lab axes, frame-averaged spectrum, least squares
    out["lab_lsq"] = fit_lsq(qx ** 2, qy ** 2, A / Pm.mean(0))
    # 2. maximum likelihood, lab axes (every frame uses q_par = q_x, q_perp = q_y)
    ones = np.ones((nf, 1))
    out["lab_ml"] = fit_ml((ones * qx ** 2).ravel(), (ones * qy ** 2).ravel(), Pm.ravel(), A)
    # 3. maximum likelihood, director frame: rotate each frame's wavevectors by that frame's theta0
    qpar, qperp = rotate(qx[None, :], qy[None, :], th0[:, None])     # shape (frame, mode)
    out["dir_ml"] = fit_ml((qpar ** 2).ravel(), (qperp ** 2).ravel(), Pm.ravel(), A)
    # 4. cross-check: director frame with the original estimator, frames grouped by their mean angle
    b = np.floor(np.degrees(th0) / BIN_DEG).astype(int); X1, X3, Y = [], [], []
    for k in np.unique(b):
        sel = b == k
        if sel.sum() < MIN_BIN_FRAMES: continue
        qp, qn = rotate(qx, qy, th0[sel].mean())
        X1.append(qn ** 2); X3.append(qp ** 2); Y.append(A / Pm[sel].mean(0))
    out["dir_lsq"] = fit_lsq(np.concatenate(X3), np.concatenate(X1), np.concatenate(Y))
    out["dir_lsq_frames"] = int(sum((b == k).sum() for k in np.unique(b) if (b == k).sum() >= MIN_BIN_FRAMES))

    # K_eff(q) along the director (bend) and across it (splay), director frame, binned in |q|
    ang = np.abs(np.arctan2(qperp, qpar))                            # angle of q from the director, in [0, pi]
    ang = np.minimum(ang, np.pi - ang)                               # fold: q and -q are the same direction
    qa = np.sqrt(qpar ** 2 + qperp ** 2); dq = 2 * np.pi / L
    for key, on in (("bend", ang < ANG_WIN), ("splay", np.abs(ang - np.pi / 2) < ANG_WIN)):
        rows = []
        for j in range(1, int(qmax / dq) + 1):
            s = on & (np.abs(qa - j * dq) < dq / 2)
            if s.any(): rows.append((j * dq, A / ((j * dq) ** 2 * Pm[s].mean())))
        out[key] = rows
    # keys used by the figure script: the quoted (director-frame ML) constants
    out["K1"], out["K3"], out["c"] = out["dir_ml"]
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("runs", nargs="+"); ap.add_argument("--qmax", type=float, default=0.8)
    ap.add_argument("--t0", type=float, default=200.0); a = ap.parse_args()
    for r in a.runs:
        d = analyse(r, a.qmax, a.t0)
        print(f"{r}: L {d['L']}, {d['frames']} frames, <S> {d['S']:.3f}; mean director angle {d['th0_mean']:+.1f} deg, sd {d['th0_sd']:.1f}, max |.| {d['th0_max']:.1f}")
        for k in ("lab_lsq", "lab_ml", "dir_lsq", "dir_ml"):
            K1, K3, c = d[k]
            print(f"   {k:8s} K1 (splay) = {K1:6.2f}  K3 (bend) = {K3:6.2f}  K1/K3 = {K1 / K3:.3f}  c = {c:.3g}" + (f"  ({d['dir_lsq_frames']} frames in bins)" if k == "dir_lsq" else ""))
        print("   K_eff along the director (bend):", " ".join(f"{qq:.2f}:{k:.2f}" for qq, k in d["bend"][:6]))
        print("   K_eff across the director (splay):", " ".join(f"{qq:.2f}:{k:.2f}" for qq, k in d["splay"][:6]))


if __name__ == "__main__":
    main()
