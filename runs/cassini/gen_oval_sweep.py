#!/usr/bin/env python3
"""
DRAFT oval sweep for the server (not launched; needs the user's approval, decision of 2026-10-02).

Grid (matched variables, user's idea 2): shape b/a, size R_eq = sqrt(A/pi) in units of l_a, linear defect density
rho = 2k/P (k rounded to the nearest integer, rho recomputed), 3 seeds, zeta = 0.07, eMode 3, 10k + 100k steps,
defects/20, director + flow/1000 (as the circle sweep). Inputs are written with gen_cassini.py.
Usage: gen_oval_sweep.py OUTROOT [--dry]   (--dry: print the manifest and cost estimate only)
"""
import argparse, math, os, subprocess, sys
import numpy as np

LA = 22.21
BA = [3.0, 2.5, 2.25, 2.0, 1.8, 1.6, 1.45, 1.3, 1.1]
REQ = [1.0, 1.5, 2.0]
RHO = [0.04, 0.06, 0.09, 0.12]
SEEDS = [1, 2, 3]


def shape(a, b, n=20000):
    th = np.linspace(0, 2 * np.pi, n, endpoint=False)
    r2 = a * a * np.cos(2 * th) + np.sqrt(b ** 4 - a ** 4 * np.sin(2 * th) ** 2)
    x, y = np.sqrt(r2) * np.cos(th), np.sqrt(r2) * np.sin(th)
    A = 0.5 * r2.mean() * 2 * np.pi
    P = np.sum(np.hypot(np.diff(np.r_[x, x[0]]), np.diff(np.r_[y, y[0]])))
    return A, P, 2 * math.hypot(a, b) / (2 * y.max())


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("outroot"); ap.add_argument("--dry", action="store_true")
    a = ap.parse_args()
    rows = []
    for ba in BA:
        A1, P1, asp = shape(1 / math.sqrt(1 + ba * ba), ba / math.sqrt(1 + ba * ba))      # half-length 1
        for req in REQ:
            H = req * LA * math.sqrt(math.pi / A1)                                          # scale to sqrt(A/pi) = req l_a
            A, P = A1 * H * H, P1 * H
            av = H / math.sqrt(1 + ba * ba); bv = ba * av
            for rho in RHO:
                k = max(1, round(rho * P / 2)); rho_eff = 2 * k / P
                if rho_eff > 0.15:
                    continue
                for s in SEEDS:
                    rows.append(dict(ba=ba, aspect=asp, req=req, H=H, a=av, b=bv, k=k, rho=rho_eff, seed=s, pop=round(20 * A)))
    cost = sum(r["pop"] for r in rows) * 110000
    print(f"{len(rows)} runs; total work {cost / 1e12:.1f} x 1e12 particle-steps "
          f"(server: ~2.2 per day at 96-way -> ~{cost / 1e12 / 2.2:.1f} days); largest {max(r['pop'] for r in rows)} particles")
    print("b/a  aspect  R_eq/l_a  half-length   k   rho")
    seen = set()
    for r in rows:
        key = (r["ba"], r["req"], r["k"])
        if key in seen:
            continue
        seen.add(key)
        print(f"{r['ba']:4.2f}  {r['aspect']:5.2f}   {r['req']:4.2f}     {r['H']:6.1f}     {r['k']:3d}  {r['rho']:.3f}")
    if a.dry:
        return
    os.makedirs(a.outroot, exist_ok=True)
    jobs = open(os.path.join(a.outroot, "jobs_oval_sweep.txt"), "w")
    for r in sorted(rows, key=lambda r: -r["pop"]):                       # largest first
        d = os.path.join(a.outroot, f"ba{r['ba']:.2f}", f"Req{r['req']:.1f}", f"k{r['k']:02d}", f"s{r['seed']:02d}")
        seed = int(200000 + 1000 * r["seed"] + 100 * BA.index(r["ba"]) + 10 * REQ.index(r["req"]) + RHO.index(min(RHO, key=lambda x: abs(x - r["rho"]))))
        subprocess.run([sys.executable, os.path.join(os.path.dirname(os.path.abspath(__file__)), "gen_cassini.py"), d, "--k", str(r["k"]),
                        "--b", f"{r['b']:.4f}", "--a", f"{r['a']:.4f}", "--emode", "3", "--act", "0.07", "--seed", str(seed),
                        "--warmup", "10000", "--steps", "100000", "--defout", "20", "--dirout", "1000", "--flowout", "1000"],
                       check=True, stdout=subprocess.DEVNULL)
        jobs.write(os.path.abspath(d) + "\n")
    jobs.close()
    print("written", os.path.join(a.outroot, "jobs_oval_sweep.txt"))


if __name__ == "__main__":
    main()
