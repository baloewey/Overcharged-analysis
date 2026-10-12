"""Summary of kinetics_carriers.csv -> kinetics_carriers_summary.txt (numpy only)."""
import csv, numpy as np
rows = list(csv.DictReader(open("kinetics_carriers.csv")))
f = lambda k: np.array([float(r[k]) if r[k] not in ("", "nan") else np.nan for r in rows])
D = {k: f(k) for k in rows[0] if k not in ("run", "family")}
fam = np.array([r["family"] for r in rows]); asp = D["aspect"]; R = D["Req_la"]; el = asp >= 1.25
D["g_esc"] = D["r_esc"] / D["nwall"]; D["g_ref_vac"] = D["r_ref"] / np.clip(D["vac"], 0.05, None); D["net_bulk"] = 0.5 * (D["nbp"] - D["nbm"])
out = []; say = lambda s: (out.append(s), print(s))
def tab(groups, keys, title):
    say(title); say("%-24s" % "group" + "".join("%10s" % k[:9] for k in keys))
    for name, m in groups:
        if m.sum(): say("%-24s" % name + "".join("%10.3f" % np.nanmean(D[k][m]) for k in keys))
K1 = ["vac", "nwall", "r_esc", "r_ref", "g_esc", "g_ref_vac", "p_esc", "p_ref", "esc_tip", "esc_mid", "esc_flank", "tip_share", "ref_tip"]
K2 = ["vac", "nbp", "nbm", "net_bulk", "unpaired", "paired_frac", "unp_x", "unp_rho"]
G = [("%s R%.1f" % ("elong" if e else "circ-like", r), (el == e) & (R == r)) for e in (False, True) for r in (1.0, 1.5, 2.0)]
tab(G, K1, "== kinetics by class and R_eq (rates per 1000 t.u.; g_esc per wall defect, g_ref_vac per vacancy; p per 20 t.u.) ==")
tab(G, K2, "\n== bulk carriers by class and R_eq (unpaired = bulk +1/2 with no -1/2 within 4 cells) ==")
GA = [("%s %.2f" % (fa, a), (fam == fa) & (asp == a)) for fa in ("circle", "ellipse", "oval") for a in sorted(set(asp[fam == fa]))]
tab(GA, ["vac", "g_esc", "g_ref_vac", "p_ref", "esc_tip", "esc_flank", "net_bulk", "unpaired", "paired_frac", "unp_x", "unp_rho"], "\n== by family and aspect (all R_eq) ==")
say("\n== correlations across runs ==")
for name, m in (("circle-like", ~el), ("elongated", el)):
    c = lambda a, b: np.corrcoef(D[a][m], D[b][m])[0, 1]
    say("%-12s corr(vac, g_esc) %+.2f  corr(vac, p_ref) %+.2f  corr(vac, g_ref_vac) %+.2f  corr(vac, unpaired) %+.2f  corr(vac, net_bulk) %+.2f" % (
        name, c("vac", "g_esc"), c("vac", "p_ref"), c("vac", "g_ref_vac"), c("vac", "unpaired"), c("vac", "net_bulk")))
say("\n== ratio elongated / circle-like at fixed R_eq ==")
for r in (1.0, 1.5, 2.0):
    say("R_eq %.1f: " % r + "  ".join("%s %.2f" % (k, np.nanmean(D[k][el & (R == r)]) / np.nanmean(D[k][~el & (R == r)]))
                                       for k in ("vac", "g_esc", "r_esc", "p_ref", "r_ref", "g_ref_vac", "nbp", "unpaired")))
open("kinetics_carriers_summary.txt", "w").write("\n".join(out) + "\n")
