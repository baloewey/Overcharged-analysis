#!/usr/bin/env python3
"""Trim white margins of the crops -> assets/clean/*.png. Panel letters are already absent: the crops come from the PRL_NOLABEL renders (crop_panels.py)."""
import glob, os
from PIL import Image, ImageChops
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "assets"); O = os.path.join(D, "clean"); os.makedirs(O, exist_ok=True)
for f in sorted(glob.glob(os.path.join(D, "f[1-4]*.png"))) + [os.path.join(D, n + ".png") for n in ("anchoring", "rings", "coreflow")]:
    im = Image.open(f).convert("RGB"); W, H = im.size
    bb = ImageChops.difference(im, Image.new("RGB", im.size, "white")).convert("L").point(lambda v: 255 if v > 12 else 0).getbbox()
    pad = 20; im.crop((max(bb[0] - pad, 0), max(bb[1] - pad, 0), min(bb[2] + pad, W), min(bb[3] + pad, H))).save(os.path.join(O, os.path.basename(f)))
print("ok")
