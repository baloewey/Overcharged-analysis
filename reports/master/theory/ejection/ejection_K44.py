#!/usr/bin/env python3
"""One-constant ejection energetics at the bare constant K = 44 kBT (2026-10-22). The relaxed one-constant Frank energy with rigid cores is
proportional to K, so the K = 32.5 one-constant rows of ejection_energetics.txt and ejection_neg.txt scale exactly by 44/32.5.
-> ejection_K44.txt"""
import re
F = 44.0 / 32.5; L = [f"ONE-CONSTANT ENERGETICS AT K = 44 kBT (= K 32.5 values x {F:.4f}; anisotropic rows dropped: the response is isotropic)"]
for src in ("ejection_energetics.txt", "ejection_neg.txt"):
    L.append(f"\n-- from {src}")
    for line in open(src):
        if "one-constant" not in line and "one constant" not in line: continue
        L.append(re.sub(r"([+-]\d+(?:\.\d+)?) kBT|([+-]\d+(?:\.\d+)?)(?=[;,]|$)", lambda m: (f"{float(m.group(1))*F:+.0f} kBT" if m.group(1) else f"{float(m.group(2))*F:+.0f}"), line.rstrip()))
L = [l.replace("pi K = 102 kBT", "pi K = 138 kBT") for l in L]; open("ejection_K44.txt", "w").write("\n".join(L) + "\n"); print("\n".join(L))
