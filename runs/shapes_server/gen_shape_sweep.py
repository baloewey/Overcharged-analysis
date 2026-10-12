#!/usr/bin/env python3
"""
Shape sweep for the server (decided 2026-10-05: ovals before the extra circle seeds). Extends the drafted oval sweep
(gen_oval_sweep.py, kept unchanged) with a high-density column, matched ellipses and circles in the same sweep.

Grid: equivalent radius R_eq = sqrt(A/pi) = 1.0, 1.5, 2.0 l_a; target linear defect density rho = 2k/P in
{0.04, 0.06, 0.09, 0.12, 0.14} (k rounded, rho recomputed, rho_eff <= 0.15 kept); zeta = 0.07; 10k + 100k steps;
defects every 20 steps, director and flow every 1000 (as the circle sweep).
  ovals     Cassini b/a = 3.0, 2.5, 2.25, 2.0, 1.8, 1.6, 1.45, 1.3, 1.1 (eMode 3)          gen_cassini.py
  ellipses  aspect eps = 1.2, 1.4, 1.6, 2.0, 2.5 (eMode 2, tip radius >= 5 cells checked)  gen_ellipse.py
  circles   aspect 1 (ordinary circle wall, same binary)                                    gen_cassini.py --circle
Seeds: --nseeds (default 3); seed = 80000000 + 100 * i_point + s, i_point = index of (family, shape, R_eq, k) in the
generation order: unique, and adding seeds later (--nseeds 10) only writes the missing seed directories.
Jobs file: largest first. Needs the cassini-branch binary (seed honoured).
Usage: gen_shape_sweep.py OUTROOT [--nseeds 3] [--dry]
"""
import argparse, math, os, subprocess, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
LA = 22.21
BA = [3.0, 2.5, 2.25, 2.0, 1.8, 1.6, 1.45, 1.3, 1.1]
EPS = [1.2, 1.4, 1.6, 2.0, 2.5]
REQ = [1.0, 1.5, 2.0]
RHO = [0.04, 0.06, 0.09, 0.12, 0.14]
STEPS = 110000


def cassini_unit(ba, n=20000):
    """Oval with half-length 1: area, perimeter, aspect."""
    a, b = 1 / math.sqrt(1 + ba * ba), ba / math.sqrt(1 + ba * ba)
    th = np.linspace(0, 2 * np.pi, n, endpoint=False)
    r2 = a * a * np.cos(2 * th) + np.sqrt(b ** 4 - a ** 4 * np.sin(2 * th) ** 2)
    x, y = np.sqrt(r2) * np.cos(th), np.sqrt(r2) * np.sin(th)
    A = 0.5 * r2.mean() * 2 * np.pi
    P = np.sum(np.hypot(np.diff(np.r_[x, x[0]]), np.diff(np.r_[y, y[0]])))
    return A, P, 1.0 / y.max()


def ellipse_perimeter(A, B, n=20000):
    t = np.linspace(0, 2 * np.pi, n + 1)
    return np.sum(np.hypot(np.diff(A * np.cos(t)), np.diff(B * np.sin(t))))


def points():
    pts = []
    for req in REQ:
        area = math.pi * (req * LA) ** 2
        # circles
        R = req * LA; P = 2 * math.pi * R
        for rho in RHO:
            pts.append(dict(fam="circle", shape=1.0, aspect=1.0, req=req, P=P, area=area, R=R, rho_t=rho))
        # ovals
        for ba in BA:
            A1, P1, asp = cassini_unit(ba)
            H = math.sqrt(area / A1)
            for rho in RHO:
                pts.append(dict(fam="oval", shape=ba, aspect=asp, req=req, P=P1 * H, area=area,
                                a=H / math.sqrt(1 + ba * ba), b=ba * H / math.sqrt(1 + ba * ba), rho_t=rho))
        # ellipses
        for eps in EPS:
            Ax = math.sqrt(area * eps / math.pi); By = Ax / eps
            for rho in RHO:
                pts.append(dict(fam="ellipse", shape=eps, aspect=eps, req=req, P=ellipse_perimeter(Ax, By), area=area,
                                Ax=Ax, By=By, tip=By * By / Ax, rho_t=rho))
    out, seen = [], set()
    for p in pts:
        p["k"] = max(1, round(p["rho_t"] * p["P"] / 2)); p["rho"] = 2 * p["k"] / p["P"]
        key = (p["fam"], p["shape"], p["req"], p["k"])
        if p["rho"] > 0.15 or key in seen or (p["fam"] == "ellipse" and p["tip"] < 5):
            continue
        seen.add(key); p["pop"] = round(20 * p["area"]); out.append(p)
    return out


def write(p, d, seed):
    common = ["--act", "0.07", "--seed", str(seed), "--warmup", "10000", "--steps", "100000", "--defout", "20",
              "--dirout", "1000", "--flowout", "1000"]
    if p["fam"] == "ellipse":
        cmd = [sys.executable, os.path.join(HERE, "gen_ellipse.py"), d, "--area", f"{p['area']:.4f}", "--eps", str(p["shape"]),
               "--k", str(p["k"])] + common
    elif p["fam"] == "circle":
        cmd = [sys.executable, os.path.join(HERE, "gen_cassini.py"), d, "--k", str(p["k"]), "--b", f"{p['R']:.4f}", "--a", "0",
               "--circle"] + common
    else:
        cmd = [sys.executable, os.path.join(HERE, "gen_cassini.py"), d, "--k", str(p["k"]), "--b", f"{p['b']:.4f}",
               "--a", f"{p['a']:.4f}", "--emode", "3"] + common
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("outroot"); ap.add_argument("--nseeds", type=int, default=3)
    ap.add_argument("--dry", action="store_true"); a = ap.parse_args()
    pts = points()
    runs = [(p, s) for p in pts for s in range(1, a.nseeds + 1)]
    for fam in ("circle", "oval", "ellipse"):
        sub = [r for r in runs if r[0]["fam"] == fam]
        print(f"{fam:8s}: {len({(p['shape'], p['req'], p['k']) for p, _ in sub}):3d} points, {len(sub):4d} runs, "
              f"{sum(p['pop'] for p, _ in sub) * STEPS / 1e12:.2f}e12 particle-steps")
    cost = sum(p["pop"] for p, _ in runs) * STEPS
    print(f"total {len(runs)} runs, {cost / 1e12:.2f}e12 particle-steps; largest {max(p['pop'] for p in pts)} particles")
    print("rho_eff range:", f"{min(p['rho'] for p in pts):.3f}-{max(p['rho'] for p in pts):.3f};",
          "aspects:", sorted({round(p['aspect'], 2) for p in pts}))
    if a.dry:
        return
    os.makedirs(a.outroot, exist_ok=True)
    jobs, man = [], [("dir", "family", "shape", "aspect", "Req_la", "k", "rho", "seed", "pop")]
    for i, p in enumerate(pts, start=1):
        for s in range(1, a.nseeds + 1):
            name = {"circle": "circle", "oval": f"ba{p['shape']:.2f}", "ellipse": f"eps{p['shape']:.2f}"}[p["fam"]]
            d = os.path.join(a.outroot, p["fam"], name, f"Req{p['req']:.1f}", f"k{p['k']:02d}", f"s{s:02d}")
            seed = 80000000 + 100 * i + s
            if not os.path.exists(os.path.join(d, "input.json")):
                write(p, d, seed)
            jobs.append((p["pop"], os.path.abspath(d)))
            man.append((os.path.relpath(d, a.outroot), p["fam"], p["shape"], f"{p['aspect']:.3f}", p["req"], p["k"], f"{p['rho']:.4f}", seed, p["pop"]))
    with open(os.path.join(a.outroot, "jobs_shape_sweep.txt"), "w") as f:
        for _, d in sorted(jobs, key=lambda x: -x[0]): f.write(d + "\n")
    with open(os.path.join(a.outroot, "manifest_shape_sweep.tsv"), "w") as f:
        for row in man: f.write("\t".join(map(str, row)) + "\n")
    print("written jobs_shape_sweep.txt and manifest_shape_sweep.tsv in", a.outroot)


if __name__ == "__main__":
    main()
