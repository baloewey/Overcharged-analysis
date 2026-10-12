#!/usr/bin/env python3
"""Regress measured cell velocities on the Stokes flows of the defect forces (wall +1/2, bulk +1/2, -1/2 tripoles)."""
import numpy as np, json
def load(tag):
    V = np.load(f"{tag}_vel.npz"); nf, nc = V["vx"].shape
    U = {c: np.loadtxt(f"{tag}_{c}.csv", delimiter=",").reshape(nf, nc, 2) for c in ("wall", "bulk", "minus")}
    meas = np.stack([V["vx"], V["vy"]], -1)
    return V, U, meas
def circ(x, y, v): r = np.hypot(x, y); return np.mean((x * v[..., 1] - y * v[..., 0]) / r, axis=-1)
pool = {"X": [], "y": []}
print(f"{'run':5} {'c_wall':>8} {'c_bulk':>8} {'c_minus':>8} {'R2':>6} {'R2 bulk only':>12} {'R2 wall only':>12} | circulation meas vs fit (mean |.|), corr per frame")
for tag in ("k02", "k06", "k10"):
    V, U, meas = load(tag)
    X = np.column_stack([U[c].ravel() for c in ("wall", "bulk", "minus")]); y = meas.ravel()
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    r2 = 1 - np.sum((X @ coef - y) ** 2) / np.sum(y ** 2)
    def r2_only(cols):
        Xs = X[:, cols]; c, *_ = np.linalg.lstsq(Xs, y, rcond=None); return 1 - np.sum((Xs @ c - y) ** 2) / np.sum(y ** 2)
    fit = (X @ coef).reshape(meas.shape)
    cm, cf = circ(V["x"], V["y"], meas), circ(V["x"], V["y"], fit)
    print(f"{tag:5} {coef[0]:8.2f} {coef[1]:8.2f} {coef[2]:8.2f} {r2:6.3f} {r2_only([1]):12.3f} {r2_only([0]):12.3f} | {np.abs(cm).mean():.3f} vs {np.abs(cf).mean():.3f}, corr {np.corrcoef(cm, cf)[0, 1]:.3f}")
    pool["X"].append(X); pool["y"].append(y)
X = np.vstack(pool["X"]); y = np.concatenate(pool["y"]); coef, *_ = np.linalg.lstsq(X, y, rcond=None)
print(f"pooled: c_wall {coef[0]:.2f}, c_bulk {coef[1]:.2f}, c_minus {coef[2]:.2f}, R2 {1 - np.sum((X @ coef - y) ** 2) / np.sum(y ** 2):.3f}")
