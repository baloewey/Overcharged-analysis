#!/usr/bin/env python3
"""Crop paper-figure panels (rendered at 600 dpi without panel letters: prl/figs scripts run with PRL_NOLABEL=1 PRL_OUTDIR=assets/nolabel, then pdftoppm into assets/raw_nolabel/) for the workshop deck. Boxes are fractions (x0, y0, x1, y1)."""
import os
from PIL import Image
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "assets")
BOX = {"fig1_observation": {"f1a_all": (0.0, 0.0, 0.765, 0.485), "f1a_passive": (0.035, 0.02, 0.175, 0.385), "f1b": (0.775, 0.0, 1.0, 0.36),
                            "f1c": (0.015, 0.495, 0.625, 1.0), "f1d": (0.655, 0.5, 0.975, 1.0)},
       "fig2_mechanism": {"f2a": (0.0, 0.0, 0.5, 0.405), "f2b": (0.5, 0.04, 1.0, 0.365), "f2c": (0.0, 0.42, 1.0, 1.0)},
       "fig3_tiltlaw": {"f3a": (0.0, 0.0, 0.312, 1.0), "f3b": (0.318, 0.0, 0.685, 1.0), "f3c": (0.695, 0.0, 1.0, 1.0)},
       "fig4_shapes": {"f4a": (0.0, 0.0, 0.33, 0.475), "f4b": (0.33, 0.0, 0.645, 0.475), "f4c": (0.645, 0.0, 1.0, 0.475),
                       "f4d": (0.0, 0.5, 0.232, 1.0), "f4e": (0.232, 0.5, 0.5, 1.0), "f4f": (0.5, 0.5, 0.8, 1.0), "f4g": (0.8, 0.5, 1.0, 1.0)}}
for fig, boxes in BOX.items():
    im = Image.open(os.path.join(D, "raw_nolabel", fig + ".png")).convert("RGB"); W, H = im.size
    for name, (x0, y0, x1, y1) in boxes.items():
        im.crop((int(x0 * W), int(y0 * H), int(x1 * W), int(y1 * H))).save(os.path.join(D, name + ".png"))
print("ok")
