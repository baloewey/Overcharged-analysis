#!/usr/bin/env python3
"""kick_identity.py on runs/uftrack (union-find defect tracker in the solver, R_eq 1.5: circle k 6, ellipses eps 1.6 / 2.0 k 7, 2 seeds each). Question: with a cleaner tracker,
do identities live long enough to classify the bulk +1/2 near an escape (legacy tracker: median identity life 3 frames)? -> uftrack_identity.txt"""
import os, sys, numpy as np
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H)
import kick_identity as K  # noqa
K.SW = "/home/baloewe/Physics/Overcharged/runs/uftrack"; cats = ("TRANSIENT", "DYNAMICAL", "ESCAPED", "PERSISTENT", "OTHER"); L = []
groups = {"circle": [], "ellipse": []}
for r in sorted(os.listdir(K.SW)):
    if os.path.exists(f"{K.SW}/{r}/DONE"): groups["circle" if r.startswith("circle") else "ellipse"].append(r)
for g, runs in groups.items():
    res = [K.one(r) for r in runs]; S = [o for _, o in res if o]
    if not S: L.append(f"{g}: analysis failed {[rel for rel, o in res if not o]}"); continue
    Kk = sum((o["kick"] for o in S), []); C = sum((o["ctrl"] for o in S), []); E = sum(o["esc"] for o in S)
    life = sum((o.get("life", []) for o in S), [])
    L.append(f"{g} ({len(S)} runs): {E} escapes; bulk +1/2 within 3 cells: {len(Kk) / max(E, 1):.0%} of escapes vs {len(C) / max(E, 1):.0%} control")
    L.append("   kicker: " + "  ".join(f"{c} {Kk.count(c) / max(len(Kk), 1):.0%}" for c in cats) + " | control: " + "  ".join(f"{c} {C.count(c) / max(len(C), 1):.0%}" for c in cats))
    if life: L.append(f"   identity life: median {np.median(life):.0f} frames, 90th pct {np.percentile(life, 90):.0f}")
    Rc = sum((o.get("recent", []) for o in S), []); Rcc = sum((o.get("recent_ctrl", []) for o in S), [])
    if Rc: L.append(f"   escaper was in the bulk within 40 t.u. before: {np.mean(Rc):.0%} (random wall defect {np.mean(Rcc):.0%})")
print("\n".join(L)); open(os.path.join(H, "uftrack_identity.txt"), "w").write("\n".join(L) + "\n")
