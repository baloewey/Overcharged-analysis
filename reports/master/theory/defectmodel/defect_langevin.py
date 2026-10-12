#!/usr/bin/env python3
"""Particle model of the active cavity (user's idea, 2026-10-07): every +1/2 defect has its own overdamped Langevin equation; defects interact
elastically (exact one-constant disk energy with same-sign images, which carries the anchoring) and hydrodynamically (each +1/2 is a point force
along its head in a no-slip Stokes disk); orientations are slaved to positions through the one-constant director field.
  N = 2k defects of charge 1/2, no creation/annihilation (the observed active state 2k-2 wall + 2 bulk is reachable).
  E/(pi K) = - sum_{i!=j} q^2 ln|z_i - z_j| - sum_{i,j} q^2 ln|R^2 - z_i conj(z_j)|           (kBT = 1)
  theta_rest(z_i) = sum_{j!=i} q [arg(z_i - z_j) + arg(1 - z_i conj(z_j)/R^2)]   (anchoring theta = k phi on the wall, lc.c)
  tail direction = 2 theta_rest (direction of div Q), head p_i = -exp(2 i theta_rest)
  dz_i/dt = mu F_i + v0 p_i + sum_{j!=i} (f/eta) g(z_i/R; z_j/R) . p_j + sqrt(2 mu kBT) xi_i
  g: unit-disk no-slip Green's function (green_disk.wl), scale-free in 2D.
Fluid circulation: mean azimuthal velocity on a ring set of sample points from the sum of Stokeslets.
Usage as a module: simulate(k, R, v0, fh, ...) -> dict of time series."""
import json, os
import numpy as np
from scipy.ndimage import map_coordinates
H = os.path.dirname(os.path.abspath(__file__))
_meta = json.load(open(os.path.join(H, "green_disk_meta.json")))
RS = np.array(_meta["rs"]); NR, NP = _meta["nr"], _meta["nphi"]
G = np.fromfile(os.path.join(H, "green_disk.bin"), dtype="<f8").reshape(len(RS), 2, NR, NP, 2)
def green(zt, zs, p):
    """velocity at targets zt (unit disk) from unit point forces at zs (unit disk) along unit vectors p (complex). Arrays broadcast (M,)."""
    rs = np.abs(zs); a = np.angle(zs); rot = np.exp(-1j * a)
    zt_r = zt * rot; p_r = p * rot                                   # rotate so that the source is on the +x axis
    r = np.clip(np.abs(zt_r), 0, 1); ph = np.angle(zt_r) % (2 * np.pi)
    irs = np.clip(rs / (RS[1] - RS[0]), 0, len(RS) - 1); ir = r * (NR - 1); ip = ph / (2 * np.pi) * (NP - 1)
    u = np.zeros(zt.shape, complex)
    for d, w in ((0, p_r.real), (1, p_r.imag)):
        for c in (0, 1):
            val = np.zeros(zt.shape)
            coords = np.vstack([irs, np.full_like(irs, d), ir, ip, np.full_like(irs, c)])
            val = map_coordinates(G, coords, order=1, mode="nearest")
            u += (w * val) * (1 if c == 0 else 1j)
    return u * np.exp(1j * a)                                        # rotate the velocity back
def elastic(z, R, q=0.5):
    dz = z[:, None] - z[None, :]; np.fill_diagonal(dz, 1.0)
    pair = np.conj(1 / dz); np.fill_diagonal(pair, 0)
    im = np.conj(-np.conj(z)[None, :] / (R ** 2 - z[:, None] * np.conj(z)[None, :])); np.fill_diagonal(im, 0)
    self_ = -2 * z / (R ** 2 - np.abs(z) ** 2)
    gradE = -q * q * (2 * pair.sum(1) + 2 * im.sum(1) + self_)
    return -gradE                                                    # in units of pi K
def heads(z, R, q=0.5):
    d = z[:, None] - z[None, :]; np.fill_diagonal(d, 1.0)
    im = 1 - z[:, None] * np.conj(z)[None, :] / R ** 2; np.fill_diagonal(im, 1.0)
    th = q * (np.angle(d) + np.angle(im)); np.fill_diagonal(th, 0)
    th_rest = th.sum(1)
    return -np.exp(2j * th_rest)
def spiral_profile(R, Re, rc):
    """v_theta(r) for a unit azimuthal body force 1/r on rc < r < Re (eta = 1), no slip at R, regular at 0 (effbc_stability.py)."""
    vp = lambda r: -0.5 * r * np.log(r); dvp = lambda r: -0.5 * (np.log(r) + 1)
    M = np.array([[rc, -rc, -1 / rc, 0], [1, -1, 1 / rc ** 2, 0], [0, Re, 1 / Re, -(Re - R ** 2 / Re)], [0, 1, -1 / Re ** 2, -(1 + R ** 2 / Re ** 2)]])
    a1, a2, b2, a3 = np.linalg.solve(M, [vp(rc), dvp(rc), -vp(Re), -dvp(Re)])
    def v(r):
        r = np.asarray(r, float)
        return np.where(r < rc, a1 * r, np.where(r < Re, a2 * r + b2 / np.maximum(r, 1e-9) + vp(np.maximum(r, 1e-9)), a3 * (r - R ** 2 / np.maximum(r, 1e-9))))
    return v
def simulate(k, R, v0, fh, K=23.0, mu=0.02, T=6000.0, dt=0.05, every=2.0, seed=1, init=None, depth_wall=4.0,
             spiral=0.0, rc=5.55, wall_supp=0.0, structured=False, barrier=0.0, a2=0.0, annular=0.0, tau_ref=63.0):
    """spiral: (zeta/eta) S of the effective-boundary body force (0 = off); the force is (spiral) sin(tau)/r on rc < r < R - 0.34 l,
    tau = mean tail tilt of the wall defects, l = 2 pi R / n_wall; only when the wall holds 2k - 2 (bulk +1). wall_supp: length a_w (cells) of the suppression of self-propulsion
    near the no-slip wall, factor tanh(depth/a_w) (0 = off).
    structured: the first 2k - 2 defects are held radially at the active layer depth (stiff radial spring, free along the wall); the last two
    start in the bulk and move freely but are kept 6 cells from the wall (tests the feedback given the observed structure). barrier: registry barrier B (kBT) of the anisotropy
    pinning (anisotropy_pinning.py), applied as a tangential force on each held defect from E_i = (B/(2k-2)) (1 - cos tau_i)/2,
    tau_i = its tail tilt, d tau_i / d(arc) = (2k - 2)/r (tilt lock). a2: second harmonic of the registry energy beyond first order
    (aniso2/relax_frank2.py: E(tau) = -A1 cos tau + a2 cos 2tau, a2 = +4..+7 kBT), shared equally by the held defects.
    annular: scale of the measured annular drive (annular_profile_k06.json, f_eff(r) = effective azimuthal body force reconstructed from the
    measured mean swirl of the k = 6, R = 31.83 circle, which includes all active forces); applied as annular * f_eff(r R0/R) (R0/R) *
    sin(tau)/sin(tau_ref), its axisymmetric no-slip flow advecting all defects (use with fh = 0 to avoid double counting)."""
    rng = np.random.default_rng(seed); N = 2 * k
    rs = R * ((N - 1) / (3 * N - 1)) ** (1 / (2 * N)) - 0.61
    z = rs * np.exp(2j * np.pi * (np.arange(N) + 0.5) / N) if init is None else init.copy()
    rl = None
    if structured:
        Nw = 2 * k - 2; rl = R * ((Nw - 1 + 4) / (3 * Nw - 1 + 4)) ** (1 / (2 * Nw)) - 0.61
        z = np.concatenate([rl * np.exp(2j * np.pi * np.arange(Nw) / Nw), np.array([0.3 * R, -0.3 * R]) * np.exp(0.5j)]).astype(complex)
    z += 0.3 * (rng.normal(size=N) + 1j * rng.normal(size=N))
    ring = np.concatenate([r * np.exp(1j * np.linspace(0, 2 * np.pi, 24, endpoint=False)) for r in R * np.array([0.3, 0.5, 0.7, 0.85])])
    out = dict(t=[], nwall=[], circ=[], tilt=[], z=[])
    nst = int(T / dt); ke = int(every / dt); piK = np.pi * K; prof_cache = {}
    ann_r = ann_v = None
    if annular > 0:
        pr = json.load(open(os.path.join(H, "annular_profile_k06.json"))); R0 = pr["R"]
        n = 400; rr = np.linspace(0, R, n + 1); hh = rr[1]
        f = annular * np.interp(rr * R0 / R, pr["r_eff"], pr["f_eff"]) * (R0 / R)
        Nn = n - 1; Am = np.zeros((Nn, Nn))
        for i in range(Nn):
            ri = rr[i + 1]; Am[i, i] = -2 / hh ** 2 - 1 / ri ** 2
            if i > 0: Am[i, i - 1] = 1 / hh ** 2 - 1 / (2 * hh * ri)
            if i < Nn - 1: Am[i, i + 1] = 1 / hh ** 2 + 1 / (2 * hh * ri)
        vv = np.zeros(n + 1); vv[1:n] = np.linalg.solve(Am, -f[1:n]); ann_r, ann_v = rr, vv / np.sin(np.radians(tau_ref))
    for s in range(nst):
        p = heads(z, R)
        F = piK * elastic(z, R)
        zi = np.repeat(z, N) / R; zj = np.tile(z, N) / R; pj = np.tile(p, N)
        mask = ~np.eye(N, dtype=bool).ravel()
        uu = np.zeros(N * N, complex); uu[mask] = green(zi[mask], zj[mask], pj[mask])
        adv = fh * uu.reshape(N, N).sum(1)
        selfp = v0 * p * (np.tanh((R - np.abs(z)) / wall_supp) if wall_supp > 0 else 1.0)
        vsp = 0.0
        if spiral > 0:
            w = (R - np.abs(z)) < depth_wall; nw = max(int(w.sum()), 2)
            tau_now = np.angle(np.mean(-p[w] * np.exp(-1j * np.angle(z[w])))) if w.any() else 0.0
            Re = R - 0.34 * 2 * np.pi * R / nw
            if Re > rc and abs(k - 0.5 * int(w.sum()) - 1.0) < 1e-9:      # spiral boundary only around a bulk +1 (inner winding 1)
                if nw not in prof_cache: prof_cache[nw] = spiral_profile(R, Re, rc)
                vsp = spiral * np.sin(tau_now) * prof_cache[nw](np.abs(z)) * 1j * z / np.maximum(np.abs(z), 1e-9)
        if annular > 0:
            w = (R - np.abs(z)) < depth_wall
            tau_now = np.angle(np.mean(-p[w] * np.exp(-1j * np.angle(z[w])))) if w.any() else 0.0
            vsp = vsp + np.sin(tau_now) * np.interp(np.abs(z), ann_r, ann_v) * 1j * z / np.maximum(np.abs(z), 1e-9)
        dz = (mu * F + selfp + adv + vsp) * dt + np.sqrt(2 * mu * dt) * (rng.normal(size=N) + 1j * rng.normal(size=N))
        if structured:
            Nw = 2 * k - 2; zw = z[:Nw]; e_r = zw / np.abs(zw); e_t = 1j * e_r
            if barrier > 0:
                tail_tilt = np.angle(-p[:Nw] * np.conj(e_r))
                dEdtau = (barrier / Nw) * 0.5 * np.sin(tail_tilt) - (a2 / Nw) * 2 * np.sin(2 * tail_tilt)   # E_i = (B/Nw)(1 - cos tau)/2 + (a2/Nw) cos 2tau
                Ft = -dEdtau * Nw / np.abs(zw)                                           # -dE/d(arc), dtau/darc = Nw/r
                dz[:Nw] += mu * Ft * e_t * dt
            rn = np.abs(zw + dz[:Nw]); dz[:Nw] += (rl - rn) * (zw + dz[:Nw]) / rn       # hold the layer radius
            zb = z[Nw:] + dz[Nw:]; far = np.abs(zb) > R - 6.0                           # keep the bulk pair out of the wall band
            zb[far] = (R - 6.0) * zb[far] / np.abs(zb[far]); dz[Nw:] = zb - z[Nw:]
        zn = z + dz; bad = np.abs(zn) > R - 0.3; zn[bad] = (R - 0.3) * zn[bad] / np.abs(zn[bad]); z = zn
        if s % ke == 0:
            depth = R - np.abs(z); w = depth < depth_wall
            # fluid circulation on the sample ring (sum of all Stokeslets)
            M = len(ring); ut = green(np.repeat(ring, N) / R, np.tile(z, M) / R, np.tile(p, M)).reshape(M, N).sum(1) * fh
            vth = np.mean((ring.real * ut.imag - ring.imag * ut.real) / np.abs(ring))
            if annular > 0:
                w2 = (R - np.abs(z)) < depth_wall
                if w2.any():
                    tau2 = np.angle(np.mean(-p[w2] * np.exp(-1j * np.angle(z[w2]))))
                    vth += np.sin(tau2) * np.mean(np.interp(np.abs(ring), ann_r, ann_v))
            if spiral > 0:
                w2 = (R - np.abs(z)) < depth_wall; nw2 = max(int(w2.sum()), 2); Re2 = R - 0.34 * 2 * np.pi * R / nw2
                if Re2 > rc and w2.any() and abs(k - 0.5 * int(w2.sum()) - 1.0) < 1e-9:
                    tau2 = np.angle(np.mean(-p[w2] * np.exp(-1j * np.angle(z[w2]))))
                    vth += spiral * np.sin(tau2) * np.mean(spiral_profile(R, Re2, rc)(np.abs(ring)))
            tail = -p; tl = np.angle(tail[w] * np.exp(-1j * np.angle(z[w]))) if w.any() else np.array([np.nan])
            out["t"].append(s * dt); out["nwall"].append(int(w.sum())); out["circ"].append(vth); out["tilt"].append(np.nanmean(tl)); out["z"].append(z.copy())
    return {k_: np.array(v) for k_, v in out.items()}
