"""Azimuthally averaged v_theta(r) of passive imposed-swirl runs (second half) -> profile_<run>.txt (r, v_theta, n)."""
import gzip, json, os, sys
import numpy as np
ROOT = "/home/baloewe/Physics/Overcharged/runs"
for rel in sys.argv[1:]:
    run = os.path.join(ROOT, rel); J = json.load(open(run + "/input.json")); Q = J["BC"][0]["Q"]
    fn = run + "/flowfield.dat.gz" if os.path.exists(run + "/flowfield.dat.gz") else run + "/flowfield.dat"
    op = gzip.open if fn.endswith(".gz") else open
    rows = np.loadtxt((l for l in op(fn, "rt") if len(l.split()) == 7 and l.split()[0][0].isdigit()))
    tf = J.get("swirlFlipTime", 0) or 1e9; t = rows[:, 0]; T = t.max(); sel = (t > min(T / 2, tf - 500)) & (t < tf)
    x = rows[sel, 1] + 0.5 - Q[0]; y = rows[sel, 2] + 0.5 - Q[1]; vx, vy = rows[sel, 4], rows[sel, 5]; r = np.hypot(x, y)
    vt = (x * vy - y * vx) / np.maximum(r, 1e-9); edges = np.arange(0, 33, 1.0); idx = np.digitize(r, edges) - 1
    out = [(0.5 * (edges[i] + edges[i + 1]), vt[idx == i].mean(), (idx == i).sum()) for i in range(len(edges) - 1) if (idx == i).sum() > 50]
    name = rel.replace("/", "_"); np.savetxt(f"profile_{name}.txt", out, fmt="%.5f", header=f"r v_theta n; R = {Q[2] if len(Q) > 2 else 'see input'}; window {sel.sum()} rows")
    print(name, " ".join(f"{a:.1f}:{b:+.4f}" for a, b, c in out))
