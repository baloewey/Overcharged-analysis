#!/usr/bin/env python3
"""Passive -> active switch across crowding and shape, plus director-field reruns (2026-10-22; follows runs/switchon).
(1) Two-stage runs: <name>/passive (3000 t.u. passive after 1000 warm-up, checkpoint every 1000 t.u.), then <name>/active resumed from the final
checkpoint with act 0.07 (2000 t.u.; flow every 10 t.u., director every 10 t.u., defects every 2 t.u.). Circles R 31.83, k 4, 10, 16; ellipses
aspect 1.5 and 2, area 2500, k 6 (eMode 2); initial director along x (oDist 1) or random (oDist 4).
(2) Director reruns: active stage only, from runs/pinning_ic alignX and rand_s3 checkpoints (k 6), 1500 t.u. with director every 10 t.u.
Checkpoint edits (switch.py): species activity field 11 -> 0.07; total steps; output flags FLOW (token 5) = 100, AVORI (11) = 0, ORDER (12) = 100
(the reader's order; see TODO 'Checkpoint reader/writer field-order mismatch')."""
import json, os, copy, math
T = json.load(open("/home/baloewe/Physics/Overcharged/runs/pinning_ic/alignX/input.json"))
def circle(k, od, seed):
    d = copy.deepcopy(T); d["BC"][0]["customAnch"] = k; d["species"][0]["oDist"] = od; d["seed"] = seed; d["defectTracker"] = 1
    d["comment"] = f"switch-on scan, PASSIVE stage: circle R 31.83, k {k}, oDist {od} (2026-10-22)"; return d
def ellipse(eps, od, seed, k=6):
    R = math.sqrt(2500 / (math.pi * eps)); a = eps * R; Lx, Ly = int(math.ceil(2 * a + 8)), int(math.ceil(2 * R + 8))
    bc = copy.deepcopy(T["BC"][0]); bc.update({"comment": f"No-slip elliptical cavity aspect {eps}, homeotropic custom anchoring", "Q": [Lx / 2, Ly / 2, 0],
        "aInv": [eps, 1, 1], "R": R, "customAnch": k, "eMode": 2})
    d = copy.deepcopy(T); d.update({"domain": [Lx, Ly], "BC": [bc], "seed": seed, "defectTracker": 1}); d["species"][0].update({"pop": 50000, "oDist": od})
    d["comment"] = f"switch-on scan, PASSIVE stage: ellipse aspect {eps}, area 2500, k {k}, oDist {od} (2026-10-22)"; return d
runs = {}
for k in (4, 10, 16):
    for od, tag in ((1, "x"), (4, "rand")): runs[f"circle_k{k:02d}_{tag}"] = circle(k, od, 9100 + 10 * k + od)
for eps in (1.5, 2.0):
    for od, tag in ((1, "x"), (4, "rand")): runs[f"ellipse_e{eps}_{tag}"] = ellipse(eps, od, 9300 + int(10 * eps) * 10 + od)
jobs = []
for n, d in runs.items():
    d.update({"simSteps": 30000, "warmUp": 10000, "checkpointOut": 10000, "defectsOut": 20, "dirSOut": 1000})
    os.makedirs(f"{n}/passive", exist_ok=True); os.makedirs(f"{n}/active", exist_ok=True); json.dump(d, open(f"{n}/passive/input.json", "w"), indent=1)
    a = copy.deepcopy(d); a.update({"checkpointIn": "checkpoint_in.dat", "seed": -1, "checkpointOut": 0, "comment": d["comment"].replace("PASSIVE", "ACTIVE (from the passive checkpoint)")})
    json.dump(a, open(f"{n}/active/input.json", "w"), indent=1); jobs.append(f"two {n} 50000")
for ic in ("alignX", "rand_s3"):
    n = f"dir_k06_{ic}"; os.makedirs(f"{n}/active", exist_ok=True)
    a = json.load(open(f"/home/baloewe/Physics/Overcharged/runs/pinning_ic/{ic}/input.json"))
    a.update({"checkpointIn": "checkpoint_in.dat", "seed": -1, "checkpointOut": 0, "comment": f"director-field rerun of the switch from runs/pinning_ic/{ic} (2026-10-22)"})
    json.dump(a, open(f"{n}/active/input.json", "w"), indent=1)
    jobs.append(f"one {n} 45000 /home/baloewe/Overcharged/runs/switchon/{ic}/checkpoint_in.dat")
open("jobs.txt", "w").write("\n".join(jobs) + "\n"); print(len(jobs), "jobs")
