#!/usr/bin/env python3
"""
Frank constants of the MPCD nematic (passive, density 20, mfpot 10, collOp 21, kBT = 1) from equilibrium director fluctuations.

Director field per cell (dirSOut) of a defect-free, initially x-aligned periodic box. Per frame: mean director angle theta0 from the
cell-averaged Q tensor; deviation dtheta = angle(n) - theta0 wrapped into (-pi/2, pi/2]. 2D FFT (cell area 1): theta_q = sum_r dtheta e^{-iqr}.
Equipartition of F = 1/2 int (K1 (d_y theta)^2 + K3 (d_x theta)^2) dA for a director along x:
    A kBT / <|theta_q|^2> = K1 q_y^2 + K3 q_x^2          (A = L^2)
Fitted, with a constant for the white noise of the cell estimate, over 0 < |q| < QMAX; also K_eff(q) = A kBT/(q^2 <|theta_q|^2>) along
the two axes. Frames: t >= T0 (after warm-up and relaxation). Usage: frank_constant.py RUN [RUN ...] [--qmax 0.8] [--t0 200]
"""
import argparse, gzip, os
import numpy as np


def read_frames(run, t0):
    f = os.path.join(run, "directorfield.dat")
    op = gzip.open(f + ".gz", "rt") if os.path.exists(f + ".gz") else open(f)
    frames, cur, tcur = [], [], None
    with op as fh:
        for line in fh:
            p = line.split()
            if len(p) != 8 or not p[0][0].isdigit(): continue
            t = float(p[0])
            if tcur is not None and t != tcur:
                if tcur >= t0: frames.append(np.array(cur))
                cur = []
            tcur = t; cur.append((int(p[1]), int(p[2]), float(p[4]), float(p[5]), float(p[7])))
    if cur and tcur >= t0: frames.append(np.array(cur))
    return frames


def analyse(run, qmax, t0):
    frames = read_frames(run, t0)
    L = int(max(f[:, 0].max() for f in frames[:1])) + 1
    acc = np.zeros((L, L)); S = []; nf = 0
    for f in frames:
        th = np.zeros((L, L)); ix, iy = f[:, 0].astype(int), f[:, 1].astype(int)
        ang = np.arctan2(f[:, 3], f[:, 2])
        c2, s2 = np.mean(np.cos(2 * ang)), np.mean(np.sin(2 * ang)); th0 = 0.5 * np.arctan2(s2, c2)
        d = (ang - th0 + np.pi / 2) % np.pi - np.pi / 2
        th[ix, iy] = d
        acc += np.abs(np.fft.fft2(th)) ** 2; S.append(f[:, 4].mean()); nf += 1
    P = acc / nf
    q = 2 * np.pi * np.fft.fftfreq(L)
    QX, QY = np.meshgrid(q, q, indexing="ij")
    m = (QX ** 2 + QY ** 2 > 0) & (QX ** 2 + QY ** 2 < qmax ** 2)
    A = L * L
    y = A / P[m]
    X = np.column_stack([QY[m] ** 2, QX[m] ** 2, np.ones(m.sum())])
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    # along the axes
    ax_bend = [(qq, A / (qq * qq * P[i, 0])) for i, qq in enumerate(q) if 0 < qq < qmax]
    ax_splay = [(qq, A / (qq * qq * P[0, j])) for j, qq in enumerate(q) if 0 < qq < qmax]
    return dict(L=L, frames=nf, S=np.mean(S), K1=coef[0], K3=coef[1], c=coef[2], bend=ax_bend, splay=ax_splay)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("runs", nargs="+"); ap.add_argument("--qmax", type=float, default=0.8)
    ap.add_argument("--t0", type=float, default=200.0); a = ap.parse_args()
    for r in a.runs:
        d = analyse(r, a.qmax, a.t0)
        print(f"{r}: L {d['L']}, {d['frames']} frames, <S> {d['S']:.3f};  K1 (splay) = {d['K1']:.2f} kBT,  K3 (bend) = {d['K3']:.2f} kBT,"
              f"  mean {0.5 * (d['K1'] + d['K3']):.2f} kBT;  noise constant {d['c']:.3g}")
        print("   K_eff along q_x (bend):", " ".join(f"{qq:.2f}:{k:.2f}" for qq, k in d["bend"][:6]))
        print("   K_eff along q_y (splay):", " ".join(f"{qq:.2f}:{k:.2f}" for qq, k in d["splay"][:6]))


if __name__ == "__main__":
    main()
