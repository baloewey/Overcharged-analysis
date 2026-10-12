#!/usr/bin/env python3
"""Mirror (negative-winding) shape sweep: the same 675 runs as gen_shape_sweep.py (same shapes, R_eq, k values, 3 seeds, steps and outputs) with k -> -k.
Directories OUTROOT/<family>/<shape>/Req*/km<|k|>/s0n; seeds 81000000 + 100 * i_point + s (distinct from the positive sweep).
Jobs file largest first (jobs_shape_mirror.txt) and manifest_shape_mirror.tsv (k negative). Usage: gen_shape_mirror.py OUTROOT [--nseeds 3] [--dry]"""
import argparse, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import gen_shape_sweep as G
def main():
    ap = argparse.ArgumentParser(); ap.add_argument("outroot"); ap.add_argument("--nseeds", type=int, default=3); ap.add_argument("--dry", action="store_true"); a = ap.parse_args()
    pts = G.points(); runs = [(p, s) for p in pts for s in range(1, a.nseeds + 1)]
    print(f"total {len(runs)} runs (mirror of the positive sweep), {sum(p['pop'] for p, _ in runs) * G.STEPS / 1e12:.2f}e12 particle-steps")
    if a.dry: return
    os.makedirs(a.outroot, exist_ok=True); jobs, man = [], [("dir", "family", "shape", "aspect", "Req_la", "k", "rho", "seed", "pop")]
    for i, p in enumerate(pts, start=1):
        q = dict(p); q["k"] = -p["k"]
        for s in range(1, a.nseeds + 1):
            name = {"circle": "circle", "oval": f"ba{p['shape']:.2f}", "ellipse": f"eps{p['shape']:.2f}"}[p["fam"]]
            d = os.path.join(a.outroot, p["fam"], name, f"Req{p['req']:.1f}", f"km{p['k']:02d}", f"s{s:02d}"); seed = 81000000 + 100 * i + s
            if not os.path.exists(os.path.join(d, "input.json")): G.write(q, d, seed)
            jobs.append((p["pop"], os.path.abspath(d)))
            man.append((os.path.relpath(d, a.outroot), p["fam"], p["shape"], f"{p['aspect']:.3f}", p["req"], -p["k"], f"{p['rho']:.4f}", seed, p["pop"]))
    with open(os.path.join(a.outroot, "jobs_shape_mirror.txt"), "w") as f:
        for _, d in sorted(jobs, key=lambda x: -x[0]): f.write(d + "\n")
    with open(os.path.join(a.outroot, "manifest_shape_mirror.tsv"), "w") as f:
        for row in man: f.write("\t".join(map(str, row)) + "\n")
    print("written jobs_shape_mirror.txt and manifest_shape_mirror.tsv in", a.outroot)
if __name__ == "__main__": main()
