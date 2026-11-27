#!/usr/bin/env python3
"""Talk movie run (user, 2026-11-26): active Cassini oval with a NARROW neck, b/a 1.10 (narrowest of the shape sweep) and 1.05 (peanut), k 6, same area (and so same particle number) as
runs/ovalfine/ba1.60_asp1.51_k06 (R_eq 1.5); flow + defects every 2 t.u. 2 seeds. Run on the server with the cassini binary."""
import json, os, numpy as np
base = json.load(open("../ovalfine/ba1.60_asp1.51_k06_s1/input.json")); area = base["species"][0]["pop"] / 20.0; jobs = []
for ba in (1.10, 1.05):
    th = np.linspace(0, 2 * np.pi, 20000, endpoint=False)
    def shape(a):
        b = ba * a; r2 = a * a * np.cos(2 * th) + np.sqrt(np.maximum(b ** 4 - a ** 4 * np.sin(2 * th) ** 2, 0)); r = np.sqrt(np.maximum(r2, 0)); return r * np.cos(th), r * np.sin(th)
    x, y = shape(1.0); f = 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1))); a = np.sqrt(area / f); x, y = shape(a)
    nx, ny = int(np.ceil(2 * x.max())) + 3, int(np.ceil(2 * y.max())) + 3
    print(f"a {a:.4f} b {ba * a:.4f} aspect {np.ptp(x) / np.ptp(y):.2f} neck/height {np.sqrt((ba * a) ** 2 - a ** 2) / y.max():.2f} domain {nx}x{ny}")
    for s in (1, 2):
        d = json.loads(json.dumps(base)); bc = d["BC"][0]; bc["R"] = round(ba * a, 4); bc["cassiniA"] = round(a, 4); bc["Q"] = [nx / 2, ny / 2, 0]
        d["domain"] = [nx, ny]; d["seed"] = 91300 + s; d["warmUp"] = 5000; d["simSteps"] = 23000; d["checkpointOut"] = 0
        d["comment"] = f"talk movie: Cassini b/a {ba} (narrow neck), k 6, area as ovalfine k06 (R_eq 1.5), flow + defects every 2 t.u."
        n = f"ba{ba:.2f}_k06_s{s}"; os.makedirs(n, exist_ok=True); json.dump(d, open(n + "/input.json", "w"), indent=1); jobs.append(n)
open("jobs.txt", "w").write("\n".join(jobs) + "\n")
