#!/usr/bin/env python3
"""Random walks of the passive wall lattice (user's idea, 2026-10-06): does fluid noise explain the excess spread across the wall?
Overdamped Langevin for N = 2k wall +1/2 in a disk (R = 31.83, effective wall Rw = R - 0.61), full electrostatic force (K = K3 = 23 kBT, the value
that fits the tangential spread in mc_exact.py), mobility mu = 0.03 (mobility.py), kBT = 1:
    dz_i = mu F_i dt + sqrt(2 mu kBT dt) xi_i + lam * u_i dt
u_i: advection by the cavity fluid, independent Ornstein-Uhlenbeck noise per defect with tangential/normal variance and correlation time measured
at the defect depth (fluid_noise.json, interpolated in depth); lam scales it (lam = 1: measured fluid).
Scenarios: thermal only; thermal + fluid (lam = 1); thermal + fluid with lam = 2, 4, 8 (how much would be needed); fluid only.
Output: tangential (internal) and normal rms with positions rounded to the cell grid (as the detector), mean depth, and the 1/e time of the
internal tangential modes m = 2..N/2 (median), against the measurement. Writes langevin.txt."""
import json, os
import numpy as np
H = os.path.dirname(os.path.abspath(__file__)); FL = os.path.dirname(H)
M = json.load(open(os.path.join(FL, "fluct_measured.json"))); FN = json.load(open(os.environ.get("FLUIDNOISE", os.path.join(H, "fluid_noise.json"))))
mobtxt = open(os.path.join(FL, "mobility.txt")).read()
R = 31.830989; Rw = R - 0.61; K = 23.0; MU = 0.03; Q = 33.5; q2 = 0.25
dd = np.array(sorted(float(d) for d in FN["depth"])); fv = lambda key: np.array([FN["depth"][str(d) if str(d) in FN["depth"] else f"{d:g}"][key] for d in dd])
VT, VN, TT, TN = fv("var_t"), fv("var_n"), fv("tau_t"), fv("tau_n")
def force(z):
    dz = z[:, None] - z[None, :]; np.fill_diagonal(dz, 1.0)
    pair = np.conj(1 / dz); np.fill_diagonal(pair, 0)
    im = np.conj(-np.conj(z)[None, :] / (Rw ** 2 - z[:, None] * np.conj(z)[None, :])); np.fill_diagonal(im, 0)
    self_ = -2 * z / (Rw ** 2 - np.abs(z) ** 2)
    gradE = -q2 * (2 * pair.sum(1) + 2 * im.sum(1) + self_)
    return -np.pi * K * gradE
def run(N, lam, thermal=True, T=40000.0, dt=0.02, every=1.0, seed=1):
    rng = np.random.default_rng(seed)
    rs = R * ((N - 1) / (3 * N - 1)) ** (1 / (2 * N)) - 0.61
    z = rs * np.exp(2j * np.pi * (np.arange(N) + 0.5) / N); ut = np.zeros(N); un = np.zeros(N)
    nst = int(T / dt); ke = int(every / dt); out = []
    for s in range(nst):
        d = R - np.abs(z); ehat = z / np.abs(z); that = 1j * ehat
        vt, vn = np.interp(d, dd, VT), np.interp(d, dd, VN); tt, tn = np.interp(d, dd, TT), np.interp(d, dd, TN)
        ut += -ut * dt / tt + np.sqrt(2 * vt * dt / tt) * rng.normal(size=N)
        un += -un * dt / tn + np.sqrt(2 * vn * dt / tn) * rng.normal(size=N)
        dz = MU * force(z) * dt + lam * (ut * that + un * ehat) * dt
        if thermal: dz += np.sqrt(2 * MU * dt) * (rng.normal(size=N) + 1j * rng.normal(size=N))
        zn = z + dz
        bad = np.abs(zn) >= Rw - 0.02; zn[bad] = z[bad]; z = zn
        if s % ke == 0 and s * dt > 2000: out.append(z.copy())
    return np.array(out)
def measure(Z, N):
    Zr = (np.round(Z.real + Q) - Q) + 1j * (np.round(Z.imag + Q) - Q)
    ph = np.unwrap(np.sort(np.angle(Zr) % (2 * np.pi), axis=1), axis=0); d = R - np.abs(Zr)
    dphi = ph - ph.mean(0); coll = dphi.mean(1); U = R * (dphi - coll[:, None])
    um = np.fft.fft(U, axis=1); taus = []
    for m in range(2, N // 2 + 1):
        x = um[:, m]; v = np.mean(np.abs(x) ** 2)
        ac = np.array([np.real(np.mean(x[:len(x) - L] * np.conj(x[L:]))) / v for L in range(300)])
        i = np.argmax(ac < np.exp(-1)); taus.append(i * 1.0 if i else np.nan)
    return d.mean(), U.std(), (d - d.mean(0)).std(), np.nanmedian(taus)
L = ["fluid noise (smoothed 1.5 cells) by depth: " + ", ".join(f"d {d:g}: rms_t {np.sqrt(a):.3f} rms_n {np.sqrt(b):.3f} tau_t {c:.1f} tau_n {e:.1f}" for d, a, b, c, e in zip(dd, VT, VN, TT, TN))]
print(L[-1])
for k in (2, 6, 10):
    N = 2 * k; m = M[f"circle_R31.83_k{k:02d}"]
    import re
    mt = re.search(rf"circle_R31.83_k{k:02d}:.*?per mode m: (.*?)  ->", mobtxt)
    tmeas = np.median([float(x.split(",")[1]) for x in mt.group(1).split(";")[1:]]) if mt else np.nan
    L.append(f"k={k:2d} N={N:2d} measured: depth {np.mean(m['depth']):.2f}, tangential {np.sqrt(np.mean(np.square(m['sig_t_internal']))):.2f}, "
             f"normal {np.sqrt(np.mean(np.square(m['sig_n']))):.2f}, mode time (m>=2, median) {tmeas:.0f}")
    print(L[-1])
    for lab, lam, th in (("thermal only", 0.0, True), ("thermal + fluid x1", 1.0, True), ("thermal + fluid x2", 2.0, True),
                         ("thermal + fluid x4", 4.0, True), ("thermal + fluid x8", 8.0, True), ("fluid only x1", 1.0, False)):
        Z = run(N, lam, th, T=20000.0 if N > 12 else 30000.0)
        dm, st, sn, tm = measure(Z, N)
        L.append(f"      {lab:20s}: depth {dm:.2f}, tangential {st:.2f}, normal {sn:.2f}, mode time {tm:.0f}"); print(L[-1])
open(os.path.join(H, "langevin.txt"), "w").write("\n".join(L) + "\n")
