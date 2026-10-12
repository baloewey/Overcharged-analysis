#!/usr/bin/env python3
"""Copy the figures used by master.tex from runs/ into reports/master/figs/ (copies; sources untouched)."""
import os, shutil
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
FIGS = os.path.join(HERE, "..", "figs")
C, K, W, A, S = "runs/cavity", "runs/cassini", "runs/cavity/analysis_weekend", "runs/cavity/analysis_spatial", "runs/night_2026-10-01"
MAP = {
    "active_length.png": "runs/activeLength/active_length_vs_activity.png",
    "cassini_validation.png": f"{K}/tests_T1T2.png",
    "cassini_passive_pinched.png": f"{K}/P2/passive_check.png",
    "circle_maps_k.png": f"{A}/report/fig1_maps.png",
    "circle_radial.png": f"{A}/report/fig2_radial.png",
    "circle_sitedist.png": f"{A}/report/fig3_sitedist.png",
    "circle_bulkcharge_rules.png": f"{A}/report/fig4_bulkcharge.png",
    "g_same.png": f"{A}/report/fig5_g_same.png",
    "g_pm.png": f"{A}/report/fig6_g_pm.png",
    "g_act.png": f"{A}/report/fig7_g_act.png",
    "g_zones.png": f"{A}/rscan_sites_gr_zones.png",
    "edge_criteria_circles.png": f"{A}/edge_criteria_circles.png",
    "layer_diagram_circles.png": f"{W}/layer_diagram_circles.png",
    "knee_check.png": f"{W}/knee_check.png",
    "crit_bulk_maps.png": f"{W}/crit_p1_bulk_maps.png",
    "crit_rstar_fstar_knee.png": f"{W}/crit_p2_rstar_fstar_knee.png",
    "crit_transitions.png": f"{W}/crit_p3_density_transitions.png",
    "crit_bulk_vs_R.png": f"{W}/crit_p4_bulk_vs_R.png",
    "mirror_bulk_vs_k.png": f"{W}/mirror_bulk_vs_k.png",
    "mirror_Q_check.png": f"{W}/mirror_Q_check.png",
    "fourier_static_circles.png": f"{C}/analysis_fourier/fourier_p1_static.png",
    "fourier_dynamic_circles.png": f"{C}/analysis_fourier/fourier_p2_dynamic.png",
    "flow_mean.png": f"{C}/analysis_flow/flow_p1_mean_flow.png",
    "flow_enstrophy.png": f"{C}/analysis_flow/flow_p2_enstrophy.png",
    "flow_tilt.png": f"{C}/analysis_flow/flow_p3_tilt.png",
    "labels_circles_ovals.png": f"{A}/bc_labels_k06.png",
    "labels_small_ovals.png": f"{S}/bc_labels_P3.png",
    "labels_lobe_oval.png": f"{K}/P4_lobe_ref/bc_labels.png",
    "accumulated_charge.png": f"{K}/fig_accumulated_charge.png",
    "layer_depth.png": f"{K}/fig_layer_depth.png",
    "layer_area_all.png": f"{K}/fig_layer_area_all.png",
    "highk.png": f"{K}/fig_highk.png",
    "lobe_profiles_large.png": f"{K}/lobe_profiles_k06.png",
    "lobe_profiles_small.png": f"{S}/lobe_profiles_P3.png",
    "edge_criteria_small_ovals.png": f"{S}/edge_criteria_P3.png",
    "wf_maps_size.png": f"{K}/W1_size/wall_fourier_maps.png",
    "wf_S_size.png": f"{K}/W1_size/wall_fourier_S.png",
    "wf_maps_kscan.png": f"{K}/W2_kscan/wall_fourier_maps.png",
    "wf_S_kscan.png": f"{K}/W2_kscan/wall_fourier_S.png",
    "wf_maps_ellipses.png": f"{K}/W3_ellipse/wall_fourier_maps.png",
    "wf_S_ellipses.png": f"{K}/W3_ellipse/wall_fourier_S.png",
    "wall_dynamics.png": f"{K}/fig_wall_dynamics.png",
    "snap_ba1.6.png": f"{K}/movies/H48_ba1.6_k06_last.png",
    "snap_ba1.1.png": f"{K}/movies/H48_ba1.1_k06_last.png",
    "snap_ba1.02.png": f"{K}/movies/H48_ba1.02_k06_last.png",
    "snap_ba1.02_k18.png": f"{K}/movies/H48_ba1.02_k18_last.png",
    "snap_small_ba2.0.png": f"{K}/movies/P3_s01/s01_last.png",
}
os.makedirs(FIGS, exist_ok=True)
for dst, src in MAP.items():
    s = os.path.join(ROOT, src)
    if os.path.exists(s): shutil.copy2(s, os.path.join(FIGS, dst))
    else: print("MISSING", src)
print(len(MAP), "figures checked")
