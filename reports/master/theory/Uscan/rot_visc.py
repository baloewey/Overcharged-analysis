#!/usr/bin/env python3
"""Effective rotational viscosity of the passive MPCD nematic vs mfpot (2026-10-27, user: is changing U a clean change of K?).
Data: the passive L64 Frank boxes (runs/Kcalib/L64_s1, s2 at U 10; runs/Uscan/Kcalib_U5_s*, Kcalib_U20_s* ), director every 2 t.u., t >= 200.
Overdamped director modes: theta_q(t) relaxes as exp(-Gamma_q t) with Gamma_q = K_mode q^2 / gamma and static power <|theta_q|^2> = A kBT / (K_mode q^2), so
    gamma = A kBT / (Gamma_q <|theta_q|^2>)          (mode by mode; needs neither K nor the fit offset)
and K_mode = A kBT / (q^2 <|theta_q|^2>). With backflow on, gamma is the flow-renormalised (effective) rotational viscosity of that mode; splay (q across the
director) and bend (q along it) are reported separately. Gamma_q: slope of -ln C(tau) for C > 0.35 (C = normalised Re autocorrelation of theta_q).
Modes: 0.09 < |q| < 0.6, within 20 deg of the director (bend) or of its normal (splay), director = run-mean angle (frames with |theta0 - mean| > 10 deg dropped).
-> rot_visc.txt"""
import os, sys, numpy as np
sys.path.insert(0, "/home/baloewe/Physics/Overcharged/runs/Kcalib")
import frank_constant as F  # noqa
R = "/home/baloewe/Physics/Overcharged/runs"
RUNS = [(5, R + "/Uscan/Kcalib_U5_s1"), (5, R + "/Uscan/Kcalib_U5_s2"), (10, R + "/Kcalib/L64_s1"), (10, R + "/Kcalib/L64_s2"),
        (20, R + "/Uscan/Kcalib_U20_s1"), (20, R + "/Uscan/Kcalib_U20_s2")]
DT = 2.0; WIN = np.radians(20); CUT = 0.35


def one(run):
    fr = F.read_frames(run, 200.0); L = int(fr[0][:, 0].max()) + 1; A = L * L; th0 = []; TQ = []
    for f in fr:
        ix, iy = f[:, 0].astype(int), f[:, 1].astype(int); ang = np.arctan2(f[:, 3], f[:, 2])
        t0 = 0.5 * np.arctan2(np.mean(np.sin(2 * ang)), np.mean(np.cos(2 * ang))); d = (ang - t0 + np.pi / 2) % np.pi - np.pi / 2
        th = np.zeros((L, L)); th[ix, iy] = d; TQ.append(np.fft.fft2(th)); th0.append(t0)
    th0 = np.array(th0); TQ = np.array(TQ); mean0 = 0.5 * np.angle(np.mean(np.exp(2j * th0)))
    keep = np.abs((th0 - mean0 + np.pi / 2) % np.pi - np.pi / 2) < np.radians(10)
    # contiguous stretches only (autocorrelation needs consecutive frames): use the longest kept run of frames
    best, cur, start = (0, 0), 0, 0
    for i, k in enumerate(keep):
        cur = cur + 1 if k else 0
        if cur and cur > best[1] - best[0]: best = (i - cur + 1, i + 1)
    TQ = TQ[best[0]:best[1]]; nf = len(TQ)
    q = 2 * np.pi * np.fft.fftfreq(L); QX, QY = np.meshgrid(q, q, indexing="ij"); qa = np.hypot(QX, QY)
    qpar = QX * np.cos(mean0) + QY * np.sin(mean0); qang = np.arccos(np.clip(np.abs(qpar) / np.maximum(qa, 1e-12), 0, 1))
    out = {"S": np.mean([f[:, 4].mean() for f in fr]), "frames": nf, "drop": int((~keep).sum())}
    for key, sel in (("bend", qang < WIN), ("splay", qang > np.pi / 2 - WIN)):
        m = sel & (qa > 0.09) & (qa < 0.6); g, Km, qq = [], [], []
        for i, j in zip(*np.nonzero(m)):
            x = TQ[:, i, j]; P = np.mean(np.abs(x) ** 2); maxlag = min(nf // 4, 200)
            C = np.array([np.mean((x[l:] * np.conj(x[:nf - l])).real) for l in range(maxlag)]) / P
            n = np.argmax(C < CUT) if (C < CUT).any() else maxlag
            if n < 3: continue                                    # decays within < 2 frames: not resolved
            lag = np.arange(n) * DT; G = -np.sum(lag * np.log(C[:n])) / np.sum(lag * lag)
            g.append(A / (G * P)); Km.append(A / (qa[i, j] ** 2 * P)); qq.append(qa[i, j])
        out[key] = (np.median(g), np.percentile(g, 25), np.percentile(g, 75), np.median(Km), len(g))
    return out


if __name__ == "__main__":
    L = ["EFFECTIVE ROTATIONAL VISCOSITY vs mfpot (passive L64, gamma = A kBT / (Gamma_q <|theta_q|^2>), median [IQR] over modes 0.09 < q < 0.6)",
         "U | run | <S> | frames used | gamma_splay [IQR] | K_splay (mode median) | n | gamma_bend [IQR] | K_bend | n | K/gamma splay, bend (diffusivity)"]
    for U, run in RUNS:
        o = one(run); s, b = o["splay"], o["bend"]
        L.append(f"{U:2d} | {os.path.basename(run)} | {o['S']:.3f} | {o['frames']} (dropped {o['drop']}) | {s[0]:6.1f} [{s[1]:.1f}-{s[2]:.1f}] | {s[3]:5.1f} | {s[4]} | "
                 f"{b[0]:6.1f} [{b[1]:.1f}-{b[2]:.1f}] | {b[3]:5.1f} | {b[4]} | {s[3] / s[0]:.3f}, {b[3] / b[0]:.3f}")
        print(L[-1], flush=True)
    open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "rot_visc.txt"), "w").write("\n".join(L) + "\n")
