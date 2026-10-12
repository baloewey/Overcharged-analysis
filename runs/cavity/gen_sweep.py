#!/usr/bin/env python3
"""
Main circular-cavity sweep: bare linear defect density rho x winding k x seeds.

For each rho in the grid and each k = 1..kmax, the radius is R = k/(pi rho);
parameter sets with R < 5 are skipped. Directory layout:

    OUT/rho0.10/k05/s01/input.json      (--sign -1: OUT/rho0.10/km05/s01, customAnch = -5)

Output policy (per parameter set):
    seed 1       defects + director field + flow field + enstrophy
    seeds 2..n   defects (+ director field unless --dirout-first-only)
                 (+ flow field with --flow-all, for the seed-averaged circulation order parameter)

Seeds are deterministic and unique: 100000*i_rho + 1000*k + s (+ 5000000 for
--sign -1, so mirror runs never share seeds with the positive set; + --seed-offset,
e.g. 9000000 for the passive reference set), where i_rho is
the index of rho in the grid (1-based).

Also writes:
    OUT/manifest.csv   one line per run: parameters, estimated CPU hours, storage
    OUT/jobs.txt       run directories, most expensive first (for run_queue.sh)

Adding seeds later (e.g. --nseeds 10) only creates the missing seed directories.
Existing inputs are never overwritten unless --force is given.

Usage:
    gen_sweep.py OUT [--rhos 0.05 0.10 ...] [--kmax 20] [--nseeds 3] [--act 0.07]
                     [--dry-run]
Passive reference set (same geometry, activity 0, 1 seed):
    gen_sweep.py OUT --act 0 --nseeds 1 --warmup 10000 --steps 30000 --flowout 0 \
                     --enstout 0 --seed-offset 9000000
"""
import argparse
import math
import os

import numpy as np

from gen_cavity import R_MIN, geometry, make_input, write_input

# Measured single-run speeds (activeNem2D-like input, Ryzen 9 7950X), particles -> MPUPS.
# Several concurrent runs share memory bandwidth, so real throughput is lower.
_SPEED_N = [5e4, 2e5, 8e5]
_SPEED_MPUPS = [2.0, 1.3, 0.63]
_BYTES_PER_CELL_LINE = {"dirSOut": 90, "flowOut": 110}


def mpups(N):
    return float(np.interp(np.log(max(N, 1e3)), np.log(_SPEED_N), _SPEED_MPUPS)) * 1e6


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out")
    ap.add_argument("--rhos", type=float, nargs="+",
                    default=[round(0.05 * i, 2) for i in range(1, 11)])
    ap.add_argument("--kmax", type=int, default=20)
    ap.add_argument("--rmax", type=float, default=float("inf"),
                    help="skip parameter sets with R > rmax (cost control at low rho)")
    ap.add_argument("--seed-offset", type=int, default=0,
                    help="added to every seed (use a distinct offset per sweep, e.g. 9000000 for passive)")
    ap.add_argument("--sign", type=int, choices=[1, -1], default=1,
                    help="-1: mirror runs with negative winding customAnch = -k")
    ap.add_argument("--nseeds", type=int, default=3)
    ap.add_argument("--act", type=float, default=0.07)
    ap.add_argument("--warmup", type=int, default=10000)
    ap.add_argument("--steps", type=int, default=100000)
    ap.add_argument("--defout", type=int, default=50)
    ap.add_argument("--dirout", type=int, default=1000)
    ap.add_argument("--dirout-first-only", action="store_true",
                    help="director field only for seed 1")
    ap.add_argument("--flowout", type=int, default=1000, help="seed 1 only (all seeds with --flow-all)")
    ap.add_argument("--flow-all", action="store_true", help="flow field for every seed")
    ap.add_argument("--enstout", type=int, default=100, help="seed 1 only")
    ap.add_argument("--checkpoint", type=int, default=10000)
    ap.add_argument("--tag", default="",
                    help="suffix for manifest/jobs files (manifest_TAG.csv, jobs_TAG.txt) so several "
                         "batches can share one OUT tree without overwriting each other")
    ap.add_argument("--dry-run", action="store_true", help="print summary, write nothing")
    ap.add_argument("--force", action="store_true", help="overwrite existing inputs")
    a = ap.parse_args()

    runs = []
    for i_rho, rho in enumerate(a.rhos, start=1):
        for k in range(1, a.kmax + 1):
            geo = geometry(k, rho=rho)
            if geo["R"] < R_MIN or geo["R"] > a.rmax:
                continue
            for s in range(1, a.nseeds + 1):
                first = (s == 1)
                dirout = a.dirout if (first or not a.dirout_first_only) else 0
                seed = (100000 * i_rho + 1000 * k + s + (5000000 if a.sign < 0 else 0)
                        + a.seed_offset)
                inp, geo = make_input(
                    a.sign * k, rho=rho, act=a.act, seed=seed,
                    warmup=a.warmup, steps=a.steps, defout=a.defout, dirout=dirout,
                    flowout=a.flowout if (first or a.flow_all) else 0,
                    enstout=a.enstout if first else 0, checkpoint=a.checkpoint)
                d = os.path.join(a.out, f"rho{rho:.2f}", f"k{k:02d}" if a.sign > 0 else f"km{k:02d}", f"s{s:02d}")
                cells = geo["L"] ** 2
                nframes = a.steps // np.array([dirout or 10**18, (a.flowout if (first or a.flow_all) else 0) or 10**18])
                storage = cells * (nframes[0] * _BYTES_PER_CELL_LINE["dirSOut"]
                                   + nframes[1] * _BYTES_PER_CELL_LINE["flowOut"])
                cpu_h = (a.warmup + a.steps) * geo["pop"] / mpups(geo["pop"]) / 3600
                runs.append({"dir": d, "rho": rho, "k": a.sign * k, "seed": s, "R": geo["R"],
                             "L": geo["L"], "pop": geo["pop"], "cpu_h": cpu_h,
                             "storage_MB": storage / 1e6, "inp": inp})

    # ---- summary
    tot_h = sum(r["cpu_h"] for r in runs)
    tot_GB = sum(r["storage_MB"] for r in runs) / 1e3
    nsets = len({(r["rho"], r["k"]) for r in runs})
    print(f"{nsets} parameter sets x {a.nseeds} seeds = {len(runs)} runs")
    print(f"estimated CPU: {tot_h:.0f} core-hours (single-run speeds; concurrent runs are slower)")
    print(f"estimated field storage (dirSOut + flowOut): {tot_GB:.1f} GB "
          f"(defect data not included, small)")
    print(f"largest run: {max(runs, key=lambda r: r['cpu_h'])['dir']} "
          f"({max(r['cpu_h'] for r in runs):.1f} h)")
    print(f"\n{'rho':>5} {'sets':>4} {'R range':>13} {'core-h':>7} {'GB':>6}")
    for rho in a.rhos:
        rr = [r for r in runs if r["rho"] == rho]
        if not rr:
            print(f"{rho:5.2f}    0  (all R < {R_MIN})")
            continue
        Rs = [r["R"] for r in rr]
        print(f"{rho:5.2f} {len({r['k'] for r in rr}):4d} {min(Rs):6.1f}-{max(Rs):6.1f} "
              f"{sum(r['cpu_h'] for r in rr):7.0f} {sum(r['storage_MB'] for r in rr)/1e3:6.1f}")
    if a.dry_run:
        return

    # ---- write
    written = skipped = 0
    for r in runs:
        path = os.path.join(r["dir"], "input.json")
        if os.path.exists(path) and not a.force:
            skipped += 1
            continue
        write_input(r["dir"], r["inp"])
        written += 1
    os.makedirs(a.out, exist_ok=True)
    sfx = f"_{a.tag}" if a.tag else ""
    with open(os.path.join(a.out, f"manifest{sfx}.csv"), "w") as f:
        f.write("dir,rho,k,seed,R,L,pop,cpu_h,storage_MB\n")
        for r in runs:
            f.write(f"{r['dir']},{r['rho']},{r['k']},{r['seed']},{r['R']:.4f},{r['L']},"
                    f"{r['pop']},{r['cpu_h']:.2f},{r['storage_MB']:.1f}\n")
    with open(os.path.join(a.out, f"jobs{sfx}.txt"), "w") as f:
        for r in sorted(runs, key=lambda r: -r["cpu_h"]):
            f.write(os.path.abspath(r["dir"]) + "\n")
    print(f"\nwrote {written} inputs, kept {skipped} existing; "
          f"manifest{sfx}.csv and jobs{sfx}.txt in {a.out}")


if __name__ == "__main__":
    main()
