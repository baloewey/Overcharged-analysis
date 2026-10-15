#!/usr/bin/env python3
"""Is the vortex row of the long winding channel a lattice? (user, 2026-10-13/15; cf. the dancing-disclination vortex lattice of Shendruk et al. 2017).
Long channels L 768, W 40 (winding m 24 per wall, charge 48; homeotropic m 0 as control). Flow every 10 t.u.; second half.
Per flow frame: vorticity omega(x) along the channel, averaged over the central band (|y - y_mid| < W/4), flow smoothed over 2 cells.
  * kymograph omega(x, t);
  * vortex count = number of sign domains of omega(x) (domains shorter than 4 cells merged); spacing distribution (domain lengths);
  * spatial order: structure factor |omega_q|^2 averaged over time; peak wavenumber and its width (a lattice gives a sharp peak);
  * persistence: correlation of omega(x, t) with omega(x - v dt, t + dt) maximised over the drift v -> correlation vs lag and the drift velocity;
  * defects: bulk +1/2 per vortex cell and their position relative to vortex centres (cell-relative coordinate).
-> fig_vortex_row.png (+pdf), vortex_row.txt"""
import os, sys, json, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, "/home/baloewe/Physics/Overcharged/runs/cavity"); sys.path.insert(0, os.path.join(H, "..", "energetics"))
from edge_criteria import read_with_angle  # noqa
from selfflow import load_flow  # noqa
R_ = "/home/baloewe/Physics/Overcharged/runs/channel_long"


def kymo(run):
    J = json.load(open(run + "/input.json")); L, Hd = J["domain"][:2]; y0, y1 = J["BC"][2]["Q"][1], J["BC"][3]["Q"][1]; ym = 0.5 * (y0 + y1); W = y1 - y0
    rows, ft, idx = load_flow(run); ix = rows[:, 1].astype(int); iy = rows[:, 2].astype(int); ts = ft[len(ft) // 2:]; ys = np.arange(Hd) + 0.5; band = np.abs(ys - ym) < W / 4
    K, VX = [], []
    for tt in ts:
        ii = idx[tt]; U = np.zeros((L, Hd)); V = np.zeros((L, Hd)); U[ix[ii], iy[ii]] = rows[ii, 4]; V[ix[ii], iy[ii]] = rows[ii, 5]
        Us, Vs = gaussian_filter(U, 2.0, mode="wrap"), gaussian_filter(V, 2.0, mode="wrap"); om = np.gradient(Vs, axis=0) - np.gradient(Us, axis=1)
        K.append(om[:, band].mean(1)); VX.append(Us[:, band].mean())
    return np.array(ts), np.array(K), J, np.array(VX)


def domains(w, minlen=4):
    s = np.sign(w); s[s == 0] = 1; L = len(s); start = np.flatnonzero(s != np.roll(s, 1))
    if not len(start): return np.array([L])
    lens = np.diff(np.r_[start, start[0] + L]); return lens[lens >= minlen]


if __name__ == "__main__":
    out = {}; Lines = ["VORTEX ROW in long channels (L 768, W 40; central band |y - y_mid| < W/4; second half, flow every 10 t.u.)"]
    for name in ("W40_m24_L768", "W40_m0_L768"):
        ts, K, J, VX = kymo(os.path.join(R_, name)); L = K.shape[1]; dt = np.median(np.diff(ts))
        nd = np.array([len(domains(k)) for k in K]); lens = np.concatenate([domains(k) for k in K])
        S = np.mean(np.abs(np.fft.rfft(K - K.mean(1, keepdims=True), axis=1)) ** 2, 0); qn = np.arange(len(S)); pk = np.argmax(S[1:]) + 1
        half = S[pk] / 2; lo = pk; hi = pk
        while lo > 1 and S[lo] > half: lo -= 1
        while hi < len(S) - 1 and S[hi] > half: hi += 1
        # persistence with drift: for lag l, max over shifts of corr(K[t], roll(K[t+l], -shift))
        lags = np.arange(1, 31); C = []; V = []
        Kz = (K - K.mean(1, keepdims=True)) / (K.std(1, keepdims=True) + 1e-12); Fk = np.fft.fft(Kz, axis=1)
        for l in lags:   # circular cross-correlation over ALL shifts (FFT), averaged over start frames; best shift = drift (signed, |shift| <= L/2)
            cc = np.real(np.fft.ifft(np.conj(Fk[:-l]) * Fk[l:], axis=1)).mean(0) / L; sh = int(np.argmax(cc)); sh = sh - L if sh > L // 2 else sh
            C.append(cc.max()); V.append(sh / (l * dt))
        C = np.array(C); tdec = lags[np.argmax(C < np.exp(-1))] * dt if (C < np.exp(-1)).any() else np.inf
        out[name] = (ts, K, S, C, lags * dt, lens)
        Lines.append(f"{name}: vortices (sign domains >= 4 cells) per frame {nd.mean():.1f} +- {nd.std():.1f} (cells of the row: {nd.mean()/2:.1f} pairs); domain length median {np.median(lens):.0f} cells, "
                     f"CV {lens.std()/lens.mean():.2f}; structure-factor peak at q = {pk} per L (wavelength {L/pk:.0f} cells = {L/pk/40:.2f} W), half width {hi-lo} modes; "
                     f"pattern correlation after 10/50/100/200 t.u. (best drift): {C[0]:.2f}/{C[4]:.2f}/{C[9]:.2f}/{C[19]:.2f}; decorrelation time {tdec:.0f} t.u.; "
                     f"best drift {np.median(V[:10]):+.3f} cells/t.u.; mean |<v_x>| band {np.mean(np.abs(VX)):.3f}")
    # defects per vortex cell (winding run)
    ts, K, S, C, lagt, lens = out["W40_m24_L768"]; run = os.path.join(R_, "W40_m24_L768"); J = json.load(open(run + "/input.json")); y0, y1 = J["BC"][2]["Q"][1], J["BC"][3]["Q"][1]
    t, fr = read_with_angle(run, np.array([0.0, 0.0]), 0.5); t = np.array(t); rel = []; perdom = []
    for i, tt in enumerate(ts):
        n = int(np.argmin(np.abs(t - tt))); f = fr[n]
        if not len(f): continue
        bulk = (f[:, 1] - y0 >= 4) & (y1 - f[:, 1] >= 4) & (f[:, 2] > 0); x = f[bulk, 0]
        w = K[i]; s = np.sign(w); s[s == 0] = 1; edges = np.flatnonzero(s != np.roll(s, 1))
        if len(edges) < 2: continue
        for a, b in zip(edges, np.r_[edges[1:], edges[0] + len(w)]):
            inside = ((x - a) % len(w)) < (b - a); perdom.append((inside.sum(), b - a))
            rel += list(((x[inside] - a) % len(w)) / (b - a))
    P = np.array(perdom); Lines.append(f"W40_m24: bulk +1/2 per vortex cell mean {P[:,0].mean():.2f} (per 40 cells of row {40*P[:,0].sum()/P[:,1].sum():.2f}); position within the cell (0 = edge, 0.5 = centre): "
                                     f"fraction in the middle half {np.mean((np.array(rel) > 0.25) & (np.array(rel) < 0.75)):.2f} (uniform 0.50)")
    fig = plt.figure(figsize=(15, 8.6)); gs = fig.add_gridspec(2, 3, height_ratios=[1.25, 1])
    for j, name in enumerate(("W40_m24_L768", "W40_m0_L768")):
        ts, K, S, C, lagt, lens = out[name]; ax = fig.add_subplot(gs[0, j * 1 + (0 if j == 0 else 1) - (1 if j else 0) + j])
    plt.close(fig)
    fig, axs = plt.subplots(2, 3, figsize=(16, 8.4), gridspec_kw=dict(height_ratios=[1.25, 1]))
    for j, (name, lab) in enumerate((("W40_m24_L768", "winding walls, charge 48"), ("W40_m0_L768", "homeotropic walls"))):
        ts, K, S, C, lagt, lens = out[name]; vm = np.percentile(np.abs(K), 97)
        ax = axs[0, j]; ax.imshow(K.T, aspect="auto", origin="lower", cmap="PuOr_r", vmin=-vm, vmax=vm, extent=[ts[0], ts[-1], 0, K.shape[1]])
        ax.set_xlabel("time (t.u.)"); ax.set_ylabel("position along the channel x (cells)"); ax.set_title(f"({'ab'[j]}) Mid-channel vorticity, {lab}", fontsize=10)
        axs[1, 0].plot(np.arange(1, len(S)), S[1:] / S[1:].max(), "-", lw=1.4, label=lab)
        axs[1, 1].plot(lagt, C, "o-", ms=3, label=lab)
        axs[1, 2].hist(lens, bins=np.arange(0, 200, 5), density=True, histtype="step", lw=1.6, label=lab)
    axs[0, 2].axis("off"); axs[0, 2].text(0, 0.98, "\n".join(Lines[1:]).replace("; ", ";\n").replace(": ", ":\n", 1), va="top", fontsize=6.5, wrap=True)
    axs[1, 0].set_xlim(0, 60); axs[1, 0].set_xlabel("wavenumber (modes per channel length)"); axs[1, 0].set_ylabel("$|\\omega_q|^2$ (normalised)"); axs[1, 0].set_title("(c) Spatial order", fontsize=10)
    axs[1, 1].set_xlabel("time lag (t.u.)"); axs[1, 1].set_ylabel("pattern correlation (best drift)"); axs[1, 1].axhline(np.exp(-1), color="0.7", ls=":"); axs[1, 1].set_title("(d) Persistence", fontsize=10)
    axs[1, 2].set_xlabel("vortex (sign-domain) length (cells)"); axs[1, 2].set_ylabel("probability density"); axs[1, 2].set_title("(e) Vortex sizes", fontsize=10)
    for a_ in axs[1]: a_.legend(fontsize=8, frameon=False); [a_.spines[s_].set_visible(False) for s_ in ("top", "right")]
    fig.tight_layout(); fig.savefig(os.path.join(H, "fig_vortex_row.png"), dpi=170)
    print("\n".join(Lines)); open(os.path.join(H, "vortex_row.txt"), "w").write("\n".join(Lines) + "\n")
