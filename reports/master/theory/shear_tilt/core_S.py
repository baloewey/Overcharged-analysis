"""Order near the wall-defect cores in the passive swirl runs (steady window): cell S vs distance from the nearest wall +1/2 (depth < 4), and the Jeffery
alignment lam_eff(S) = lam (1 + S4)/(2 S2) (Maier-Saupe S4) averaged over the cells within r_c of the core -> tilt arccos(1/<lam_eff>) for r_c = 1, 1.5, 2, 3.
Also S and lam_eff vs depth away from cores (> 4 cells from any wall defect)."""
import gzip, json, os, sys
import numpy as np
from scipy.spatial import cKDTree
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from mpcd_profile import lam_eff, S4_of, ROOT, read_with_angle  # noqa
def frames_dir(run, t0, t1):
    fn = run + "/directorfield.dat"; op = open
    if not os.path.exists(fn): fn += ".gz"; op = gzip.open
    cur, tc = [], None
    with op(fn, "rt") as f:
        for line in f:
            p = line.split()
            if len(p) != 8 or not p[0][0].isdigit(): continue
            t = float(p[0])
            if t != tc:
                if cur and t0 <= tc <= t1: yield tc, np.array(cur)
                cur, tc = [], t
                if t > t1: return
            if t0 <= t <= t1 and float(p[7]) > 0: cur.append((float(p[1]), float(p[2]), float(p[7])))
def analyse(rel, t0, t1):
    run = os.path.join(ROOT, rel); J = json.load(open(run + "/input.json")); bc = J["BC"][0]; Qc = bc["Q"]; R = bc["R"]; lam = J["species"][0].get("tumble", 2)
    td, fr = read_with_angle(run, np.array(Qc[:2], float), 0.0); dmap = {round(t, 2): f for t, f in zip(td, fr)}
    bins = np.array([0, 0.75, 1.25, 1.75, 2.5, 3.5, 5, 8]); accS = [[] for _ in bins[:-1]]; core = {rc: [] for rc in (1.0, 1.5, 2.0, 3.0)}; nfr = 0
    for k, (t, c) in enumerate(frames_dir(run, t0, t1)):
        if round(t, 2) not in dmap: continue
        f = dmap[round(t, 2)]; r = np.hypot(f[:, 0], f[:, 1]); w = f[(r > R - 4) & (f[:, 2] > 0), :2]
        if not len(w): continue
        x = c[:, 0] + 0.5 - Qc[0]; y = c[:, 1] + 0.5 - Qc[1]; dd, _ = cKDTree(w).query(np.column_stack([x, y])); S = c[:, 2]; nfr += 1
        idx = np.digitize(dd, bins) - 1
        for i in range(len(bins) - 1): accS[i].extend(S[idx == i])
        for rc in core: core[rc].extend(S[dd < rc])
    L = [f"{rel} (lam {lam}, {nfr} frames):"]
    for i in range(len(bins) - 1):
        s = np.array(accS[i])
        if len(s) < 5: continue
        L.append(f"   dist {bins[i]:4.2f}-{bins[i + 1]:4.2f}: S {s.mean():.3f} (n {len(s)})  lam_eff(mean S) {lam_eff(lam, min(s.mean(), .995)):.2f}")
    for rc, s in core.items():
        s = np.clip(np.array(s), 0.02, 0.995)[:: max(1, len(s) // 3000)]; s4 = np.array([S4_of(v) for v in s]); le = lam * (1 + s4.mean()) / (2 * s.mean())   # pooled moments of the region
        L.append(f"   cells within {rc} of a wall core: <S> {s.mean():.3f}, <S4> {s4.mean():.3f}, pooled lam_eff {le:.2f} -> arccos(1/<lam_eff>) {np.degrees(np.arccos(1 / le)) if le > 1 else np.nan:.1f} deg")
    return L
if __name__ == "__main__":
    out = []
    for rel, t0, t1 in (("swirl_threshold/passive_A0.01_s2", 1500, 3000), ("swirl_threshold/passive_A0.01_s3", 1500, 3000), ("swirl_tumble/passive_lam4.0_A0.01", 1000, 1500),
                        ("swirl_tumble/passive_lam1.2_A0.01", 1000, 1500)):
        L = analyse(rel, t0, t1); print("\n".join(L), flush=True); out += L
    open("core_S.txt", "w").write("\n".join(out) + "\n")
