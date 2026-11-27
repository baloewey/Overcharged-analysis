#!/usr/bin/env python3
"""Cavity outlines for the talk (user, 2026-11-26): circle, ellipse (aspect 2) and Cassini oval (b/a 1.30, aspect 1.94, as runs/ovalfine), equal area,
thick black stroke, transparent background. Wall shapes from energetics.wall_polyline (the simulation geometry).
-> assets/shapes_stacked.(png|pdf) (one above another) and assets/shapes_overlaid.(png|pdf) (concentric)"""
import os, sys, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
sys.path.insert(0, "/home/baloewe/Physics/Overcharged/reports/master/theory/energetics"); import energetics as En  # noqa
A = os.path.join(os.path.dirname(os.path.abspath(__file__)), "assets"); R0 = 1.0
circle = En.wall_polyline({"R": R0, "aInv": [1, 1]})[0]
e = 2.0; ellipse = En.wall_polyline({"R": R0, "aInv": [np.sqrt(e), 1 / np.sqrt(e)]})[0]                 # semi-axes sqrt(2), 1/sqrt(2): same area as the circle
ba = 1.30; oval = En.wall_polyline({"shape": 1, "cassiniA": 1.0, "R": ba})[0]
x, y = oval[:, 0], oval[:, 1]; area = 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1))); oval = oval * np.sqrt(np.pi / area)
shapes = [circle, ellipse, oval]
close = lambda W: (np.r_[W[:, 0], W[0, 0]], np.r_[W[:, 1], W[0, 1]])
fig, axs = plt.subplots(3, 1, figsize=(3.2, 5.6))
for ax, W in zip(axs, shapes):
    ax.add_patch(Polygon(W, closed=True, fill=False, ec="k", lw=4, joinstyle="round")); ax.set_aspect("equal"); ax.axis("off"); ax.set_xlim(-1.75, 1.75); ax.set_ylim(-1.1, 1.1)
fig.subplots_adjust(hspace=0.08, left=0.02, right=0.98, top=0.98, bottom=0.02)
for ext in ("png", "pdf"): fig.savefig(os.path.join(A, f"shapes_stacked.{ext}"), dpi=300, transparent=True)
fig, ax = plt.subplots(figsize=(4.2, 3.0))
for W in shapes: ax.add_patch(Polygon(W, closed=True, fill=False, ec="k", lw=4, joinstyle="round"))
ax.set_xlim(-1.75, 1.75); ax.set_ylim(-1.1, 1.1)
ax.set_aspect("equal"); ax.axis("off"); fig.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0.02)
for ext in ("png", "pdf"): fig.savefig(os.path.join(A, f"shapes_overlaid.{ext}"), dpi=300, transparent=True)
print("wrote shapes_stacked, shapes_overlaid")
