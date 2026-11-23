#!/usr/bin/env python3
"""Passive rest-pinning test at crowded spacings (user, 2026-11-23). Question: when the lock asks for registry shifts of only 1-3 cells (tau = 2 pi ds / a),
is there a pinning that holds the defects at the homeotropic sites in ABSOLUTE length (cells), i.e. do seeded tilted rings relax towards tau = 0 at small a?
Passive MPCD, no flow, circle R 31.83, eMode 2 anchoring; seeded 2k-2 wall ring at registry offset tau0/360 of the spacing (tau0 0 / 30 / 60 deg) at depth 1.6
plus the bulk pair at +-0.45 R (seed_director.py; bin/mpcd_seeddir.out, branch seeddir), as runs/rest_pinning (16 Oct) but crowded and with 3 seeds.
k 4 / 6 (spacing 33 / 20, controls) and 9 / 10 / 12 / 14 / 16 (12.5 / 11.1 / 9.1 / 7.7 / 6.7 cells). 1000 t.u., defects every 2 t.u., director every 10 t.u.
Analyse only while the wall holds 2k-2 (the passive bulk pair returns to the wall within ~100-1000 t.u.)."""
import json, os, shutil, subprocess, sys
T = json.load(open("/home/baloewe/Physics/Overcharged/runs/rest_pinning/k16_t060/input.json")); H = os.path.dirname(os.path.abspath(__file__))
SEED = "/home/baloewe/Physics/Overcharged/reports/master/theory/seeding/seed_director.py"; R = T["BC"][0]["R"]; jobs = []
cfgs = [(k, t0) for k in (16, 14, 12, 10, 9, 6, 4) for t0 in (0, 30, 60)]
def make(cfg):
    k, t0 = cfg; base = f"k{k:02d}_t{t0:03d}"; d0 = os.path.join(H, base + "_s1")
    out = []
    for s in (1, 2, 3):
        d = os.path.join(H, f"{base}_s{s}"); os.makedirs(d, exist_ok=True); J = json.loads(json.dumps(T))
        J["BC"][0]["customAnch"] = k; J.update(seed=94000000 + k * 1000 + t0 * 10 + s, dirSOut=100,
              comment=f"rest_pinning2 (2026-11-23): PASSIVE circle R 31.83, k {k}, seeded 2k-2 ring at registry offset {t0 / 360:.4f} (tau0 {t0} deg) + bulk pair; no flow; 1000 t.u.")
        J.pop("dirInitFile", None); json.dump(J, open(d + "/input.json", "w"), indent=4); out.append(d)
    bulk = f"{0.45 * R:.4f},0,0.5;{-0.45 * R:.4f},0,0.5"
    r = subprocess.run([sys.executable, SEED, d0, "--wall", str(2 * k - 2), "--bulk", bulk, "--offset", f"{t0 / 360:.6f}"], capture_output=True, text=True)
    open(d0 + "/seed.log", "w").write(r.stdout + r.stderr)
    for d in out[1:]:
        shutil.copy(d0 + "/dirinit.dat", d + "/dirinit.dat"); J = json.load(open(d0 + "/input.json")); s = int(d[-1])
        J["seed"] = 94000000 + k * 1000 + t0 * 10 + s; J["dirInitFile"] = d + "/dirinit.dat"; json.dump(J, open(d + "/input.json", "w"), indent=4)
    return [os.path.basename(d) for d in out]
if __name__ == "__main__":
    from multiprocessing import Pool
    with Pool(8) as p: res = p.map(make, cfgs)
    open(os.path.join(H, "jobs.txt"), "w").write("\n".join(n for g in res for n in g) + "\n"); print(sum(len(g) for g in res), "jobs")
