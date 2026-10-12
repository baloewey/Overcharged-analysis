#!/usr/bin/env python3
"""
Build Overcharged/results/: a copy of the project's results (figures, tables, movies, reports), ordered by question.

Nothing outside results/ is modified: every file is COPIED from runs/ or reports/ (user decision 2026-10-05: don't
touch the existing layout). Rerun after new analyses to refresh the copies; the map of the results is RESULTS.md.
Each folder gets a FILES.md (what each file is, where it was copied from); results/MANIFEST.tsv lists everything.

Usage: python3 tools/build_results.py [--dry-run]
"""
import argparse, os, shutil, sys
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results")
C, K, W, F, A, S = "runs/cavity", "runs/cassini", "runs/cavity/analysis_weekend", "runs/cavity/analysis_fourier", \
    "runs/cavity/analysis_spatial", "runs/night_2026-10-01"

FOLDERS = {
    "1_validation": "Is the method right? Custom anchoring, its resolution limit, new wall shapes, tracker, seeds.",
    "2_active_length": "Length unit: active length l_a(activity) in periodic boxes, finite size, compressibility.",
    "3_circles": "Baseline: where circles form a +1 bulk with a k - 1 wall layer.",
    "4_negative_winding": "Mirror runs (winding -k): the bulk stays positive (H5).",
    "5_activity": "Does activity enter only through l_a? (H3; working activity 0.07).",
    "6_shape": "Main result of paper 1: elongated cavities (ovals, ellipses) lose the wall layer (H1).",
    "7_flow_and_wall_dynamics": "Mechanism: circulation, wall-defect tilt, wall lattice (pinned / vacant), exchange with the bulk.",
    "8_reports": "Documents: project report, collaborator brief, other write-ups, plans.",
}

# (destination in results/, source relative to the project root, description)
MAP = [
    # ---------------------------------------------------------------- 1 validation
    ("1_validation/resolution_limit.txt", f"{C}/resolution_scan/resolution_results.txt", "Anchoring resolution limit: clean to rho = 0.15, fails from 0.25 (passive and active alike)"),
    ("1_validation/resolution_limit_active_checks.txt", f"{C}/resolution_scan/act_check.txt", "Per-run audit, active resolution scan"),
    ("1_validation/resolution_limit_passive_checks.txt", f"{C}/resolution_scan/pas_check.txt", "Per-run audit, passive resolution scan"),
    ("1_validation/cassini_walls_passive_tests.png", f"{K}/tests_T1T2.png", "Cassini wall implementation: a = 0 reproduces the circle; passive ovals with a concave neck"),
    ("1_validation/cassini_passive_pinched_oval.png", f"{K}/P2/passive_check.png", "Passive pinched ovals keep all 2k defects at the wall (13-cell neck)"),
    ("1_validation/server_binary_validation_cassini.png", f"{K}/server_validation/check.png", "Cassini binary on the server = office build (2026-10-03)"),
    ("1_validation/server_benchmark.txt", "runs/server_benchmark.txt", "Server throughput benchmark"),
    # ---------------------------------------------------------------- 2 active length
    ("2_active_length/active_length_vs_activity.png", "runs/activeLength/active_length_vs_activity.png", "l_a against activity, 100^2 and 200^2 boxes, with fit"),
    ("2_active_length/active_length_vs_activity.pdf", "runs/activeLength/active_length_vs_activity.pdf", "Same, PDF"),
    ("2_active_length/finite_size_L200_vs_L100.txt", "runs/activeLength/L200/finite_size_results.txt", "Reference value l_a(0.07) = 22.2 +- 0.3; no finite-size effect"),
    ("2_active_length/active_length_L100.txt", "runs/activeLength/L100/active_length_results.txt", "l_a for 9 activities, 100^2 box"),
    ("2_active_length/active_length_L100.csv", "runs/activeLength/L100/active_length.csv", "Data, 100^2"),
    ("2_active_length/active_length_L200.csv", "runs/activeLength/L200/active_length.csv", "Data, 200^2"),
    ("2_active_length/active_length_0.06_0.11.txt", "runs/activeLength/L200/active_length_h3high.txt", "l_a at 0.06 and 0.11 (2 seeds) for the activity study"),
    ("2_active_length/active_length_0.06_0.11.csv", "runs/activeLength/L200/active_length_h3high.csv", "Data"),
    ("2_active_length/compressibility.txt", "runs/activeLength/compressibility/compressibility_results.txt", "Density fluctuations and Mach number vs activity"),
    # ---------------------------------------------------------------- 3 circles
    ("3_circles/layer_diagram_circles.png", f"{W}/layer_diagram_circles.png", "KEY: bulk charge at the knee over (rho, R/l_a); +1 layer region 0.75-2 l_a"),
    ("3_circles/layer_diagram_circles.pdf", f"{W}/layer_diagram_circles.pdf", "Same, PDF"),
    ("3_circles/main_sweep_summary.txt", f"{W}/main_sweep_summary.txt", "Main sweep (1112 runs): knee, d*, f*, site rule, wall Fourier"),
    ("3_circles/main_sweep.csv", f"{W}/main_sweep.csv", "Main sweep, one row per run"),
    ("3_circles/knee_check.png", f"{W}/knee_check.png", "Check of the knee finder (smoothing)"),
    ("3_circles/bulk_charge_criteria/bulk_vs_R_by_criterion.png", f"{W}/crit_p4_bulk_vs_R.png", "KEY: six bulk-charge criteria vs R/l_a, +k and -k"),
    ("3_circles/bulk_charge_criteria/bulk_maps_by_criterion.png", f"{W}/crit_p1_bulk_maps.png", "Same as maps over (rho, R/l_a)"),
    ("3_circles/bulk_charge_criteria/rstar_fstar_knee.png", f"{W}/crit_p2_rstar_fstar_knee.png", "Q = +1 crossing, f*, knee depth"),
    ("3_circles/bulk_charge_criteria/density_transitions.png", f"{W}/crit_p3_density_transitions.png", "Radial defect density across knee / no-knee transitions"),
    ("3_circles/bulk_charge_criteria/criteria_maps.pdf", f"{W}/criteria_maps.pdf", "All four pages, PDF"),
    ("3_circles/bulk_charge_criteria/criteria_all.csv", f"{W}/criteria_all.csv", "Per run, every criterion"),
    ("3_circles/bulk_charge_criteria/edge_criteria_circles.png", f"{A}/edge_criteria_circles.png", "Cut-free criteria (mobility, tilt) tested against the site rule"),
    ("3_circles/fixed_radius_scan/rscan_seedavg.txt", "runs/cavity/analysis_lunch/rscan_seedavg.txt", "KEY: R = 1.43 l_a, k = 2-10, 3 seeds: bulk +1.00 +- 0.02 at every k"),
    ("3_circles/fixed_radius_scan/rscan_sites_seedavg.txt", f"{A}/rscan_sites_seedavg.txt", "Same scan, pinned-site rule: 2k - 2 sites in all 27 runs"),
    ("3_circles/fixed_radius_scan/rscan_seed1_results.txt", f"{C}/Rscan_R31.83/rscan_results.txt", "Seed 1 only (first look, density cut)"),
    ("3_circles/fixed_radius_scan/rscan_circulation.csv", f"{C}/Rscan_R31.83/circulation.csv", "Circulation order parameter per run"),
    ("3_circles/fixed_radius_scan/passive_layers.txt", "runs/cavity/analysis_lunch/layer_passive.txt", "Passive references: all 2k defects at the wall, bulk 0"),
    ("3_circles/fixed_radius_scan/tilt_vs_k.txt", f"{A}/tilt_vs_k.txt", "Wall +1/2 tilt vs k"),
    ("3_circles/defect_geography/fig1_maps.png", f"{A}/report/fig1_maps.png", "Time-averaged defect maps: wall defects sit on 2k - 2 fixed sites"),
    ("3_circles/defect_geography/fig2_radial.png", f"{A}/report/fig2_radial.png", "Radial densities: layer peak, gap, bulk"),
    ("3_circles/defect_geography/fig3_sitedist.png", f"{A}/report/fig3_sitedist.png", "Distance to the nearest site (edge rule)"),
    ("3_circles/defect_geography/fig4_bulkcharge.png", f"{A}/report/fig4_bulkcharge.png", "Bulk charge by edge rule"),
    ("3_circles/defect_geography/fig5_g_same.png", f"{A}/report/fig5_g_same.png", "Pair correlations, same sign"),
    ("3_circles/defect_geography/fig6_g_pm.png", f"{A}/report/fig6_g_pm.png", "Pair correlations, +1/2 with -1/2"),
    ("3_circles/defect_geography/rscan_sites_gr_zones.png", f"{A}/rscan_sites_gr_zones.png", "Pooled g(r) by zone (layer, outer bulk, centre)"),
    # ---------------------------------------------------------------- 4 negative winding
    ("4_negative_winding/mirror_bulk_vs_k.png", f"{W}/mirror_bulk_vs_k.png", "KEY: bulk charge vs |k|, +k and -k"),
    ("4_negative_winding/mirror_bulk_vs_k.pdf", f"{W}/mirror_bulk_vs_k.pdf", "Same, PDF"),
    ("4_negative_winding/mirror_summary.txt", f"{W}/mirror_summary.txt", "570 mirror runs: bulk 0.09 (|k| 2) -> 1.0 (|k| >= 15)"),
    ("4_negative_winding/mirror.csv", f"{W}/mirror.csv", "Per run"),
    ("4_negative_winding/mirror_Q_check.png", f"{W}/mirror_Q_check.png", "Q(r) examples: positive peak, then -|k| at the wall"),
    # ---------------------------------------------------------------- 5 activity
    ("5_activity/all_activities_summary.txt", f"{A}/H3_all_activities.txt", "KEY: zeta = 0.05/0.06/0.07/0.09/0.11 at R = 1.43 l_a, k = 2-10"),
    ("5_activity/enclosed_charge_vs_activity.txt", f"{A}/enclosed_charge.txt", "Q(s) plateau flatness vs activity (flattest at 0.07)"),
    ("5_activity/pair_correlations_vs_activity.png", f"{A}/report/fig7_g_act.png", "Bulk g(r) at three activities, collapsed in r/R"),
    ("5_activity/act0.05_3seeds.txt", f"{A}/act005_seedavg.txt", "zeta = 0.05 (onset): unstable layer, vortex fails at k = 9-10"),
    ("5_activity/act0.09_3seeds_sites.txt", f"{S}/act009_3seeds_sites.txt", "zeta = 0.09, 3 seeds: clean like 0.07"),
    ("5_activity/act0.09_seed1_summary.txt", f"{A}/act009_summary.txt", "zeta = 0.09, first seed"),
    ("5_activity/density_check.txt", f"{A}/density_check.txt", "Centrifugal density stratification in the vortex (zeta 0.07/0.09/0.11)"),
    ("5_activity/empty_cells.txt", f"{A}/empty_cells.txt", "Empty / low-occupancy cells"),
    ("5_activity/pairs_by_activity.txt", "runs/cavity/analysis_lunch/pairs_by_activity.txt", "Dynamical +-1/2 pairs per frame vs activity"),
    ("5_activity/l_a_used.txt", f"{C}/H3high/l_a_used.txt", "l_a and radii used for zeta = 0.06, 0.11"),
    ("5_activity/circulation_act0.05.txt", "runs/cavity/analysis_lunch/circ_act005_3seeds.txt", "Circulation, zeta 0.05"),
    ("5_activity/circulation_act0.06.txt", "runs/cavity/analysis_lunch/circ_act006_3seeds.txt", "Circulation, zeta 0.06"),
    ("5_activity/circulation_act0.09.txt", f"{S}/circ_act009_3seeds.txt", "Circulation, zeta 0.09"),
    ("5_activity/circulation_act0.11.txt", "runs/cavity/analysis_lunch/circ_act011_3seeds.txt", "Circulation, zeta 0.11"),
    # ---------------------------------------------------------------- 6 shape
    ("6_shape/accumulated_charge_circles_vs_ovals.png", f"{K}/fig_accumulated_charge.png", "KEY: Q(c) in circles (plateau at +1, jump at the wall) vs ovals (smooth, no plateau)"),
    ("6_shape/accumulated_charge_circles_vs_ovals.pdf", f"{K}/fig_accumulated_charge.pdf", "Same, PDF"),
    ("6_shape/layer_area_fstar_vs_aspect.png", f"{K}/fig_layer_area_all.png", "KEY: area fraction f* collapses on aspect ratio for all shapes and sizes; ellipses on the same curve"),
    ("6_shape/layer_area_fstar_vs_aspect.pdf", f"{K}/fig_layer_area_all.pdf", "Same, PDF"),
    ("6_shape/layer_depth_dstar.png", f"{K}/fig_layer_depth.png", "KEY: depth of the Q = +1 level vs elongation (size-matched control)"),
    ("6_shape/layer_depth_dstar.pdf", f"{K}/fig_layer_depth.pdf", "Same, PDF"),
    ("6_shape/high_k_pinched_oval.png", f"{K}/fig_highk.png", "A crowded wall (k 12-20) partly restores the layer"),
    ("6_shape/high_k_pinched_oval.pdf", f"{K}/fig_highk.pdf", "Same, PDF"),
    ("6_shape/matching_ovals_to_circles.txt", f"{W}/matching_table.txt", "Ovals vs circles at matched density and every size measure"),
    ("6_shape/tables/layer_area_all.txt", f"{K}/layer_area_all.txt", "f* for every shape and size"),
    ("6_shape/tables/layer_depth_k06.txt", f"{K}/layer_depth_k06.txt", "d* for k = 6, circles and ovals"),
    ("6_shape/tables/layer_depth_k03_circles.txt", f"{K}/layer_depth_k03_circles.txt", "k = 3 circles (lobe-size references)"),
    ("6_shape/tables/layer_depth_highk.txt", f"{K}/layer_depth_highk.txt", "d* in the pinched oval, k 6-20"),
    ("6_shape/tables/pinch_lobe_charge.txt", f"{K}/pinch_lobe_charge.txt", "Per-lobe charge in pinched ovals: no restored layer"),
    ("6_shape/tables/edge_criteria_ovals.txt", f"{A}/edge_criteria_apply.txt", "Cut-free criteria applied to ovals: bulk 3.1-3.9 vs 1.1-1.25 in circles"),
    ("6_shape/small_ovals_1.53la/charge_summary.txt", f"{S}/cassini_charge_P3.txt", "Small ovals (half-length 1.53 l_a, 3 seeds): sharp crossover with elongation"),
    ("6_shape/small_ovals_1.53la/edge_criteria.png", f"{S}/edge_criteria_P3.png", "Edge criteria, small ovals"),
    ("6_shape/small_ovals_1.53la/edge_criteria.txt", f"{S}/edge_criteria_P3.txt", "Same, numbers"),
    ("6_shape/small_ovals_1.53la/lobe_profiles.png", f"{S}/lobe_profiles_P3.png", "Defect density and Q(c) in lobe coordinates"),
    ("6_shape/small_ovals_1.53la/lobe_profiles.txt", f"{S}/lobe_profiles_P3.txt", "Same, numbers"),
    ("6_shape/size_series_2.16la/knee.txt", f"{K}/W1_size/knee.txt", "Half-length 2.16 l_a, 3 seeds: knee per run"),
    ("6_shape/size_series_2.16la/layer_area.txt", f"{K}/W1_size/layer_area.txt", "f* with seeds: collapse confirmed"),
    ("6_shape/size_series_2.88la/layer_area.txt", f"{K}/W4_H64seeds/layer_area.txt", "Largest ovals, seeds 2-3"),
    ("6_shape/size_series_2.88la/layer_depth.txt", f"{K}/W4_H64seeds/layer_depth.txt", "Same, d*"),
    ("6_shape/k_scan/kscan_summary.txt", f"{K}/W2_kscan/kscan_summary.txt", "k = 2-12 in small ovals: bounded excess, layer returns at high k"),
    ("6_shape/ellipses/ellipse_summary.txt", f"{K}/W3_ellipse/ellipse_summary.txt", "Ellipses (aspect 1.2/1.5/2.0): lose the layer like ovals"),
    ("6_shape/lobes_and_vortices/lobe_profiles_k06.png", f"{K}/lobe_profiles_k06.png", "Lobe-coordinate profiles, larger ovals"),
    ("6_shape/lobes_and_vortices/lobe_profiles_k06.csv", f"{K}/lobe_profiles_k06.csv", "Same, data"),
    ("6_shape/lobes_and_vortices/lobe_rotation_all.txt", f"{S}/lobe_rotation_all.txt", "Two vortices per oval; co- vs counter-rotation"),
    ("6_shape/superseded/size_vs_shape_k06_first_look.txt", f"{K}/size_vs_shape_k06.txt", "SUPERSEDED (2026-10-01): first look; 'bulk' values were cut-dependent"),
    ("6_shape/superseded/layer_area_first_version.png", f"{K}/fig_layer_area.png", "SUPERSEDED by layer_area_fstar_vs_aspect.png"),
    ("6_shape/superseded/layer_area_with_seeds_intermediate.png", f"{K}/fig_layer_area_seeds.png", "SUPERSEDED by layer_area_fstar_vs_aspect.png"),
    ("6_shape/superseded/accumulated_charge_before_small_ovals.png", "reports/project_report/figs/accumulated_charge_pre.png", "SUPERSEDED: accumulated charge before the small ovals were added"),
    # ---------------------------------------------------------------- 7 flow and wall dynamics
    ("7_flow_and_wall_dynamics/wall_dynamics_ovals_vs_circles.png", f"{K}/fig_wall_dynamics.png", "KEY: pinning, upstream drift, drift vs vacancies, k scan"),
    ("7_flow_and_wall_dynamics/wall_dynamics_ovals_vs_circles.pdf", f"{K}/fig_wall_dynamics.pdf", "Same, PDF"),
    ("7_flow_and_wall_dynamics/wall_dynamics_ovals.csv", f"{K}/wall_fourier_dyn.csv", "Per oval/ellipse run"),
    ("7_flow_and_wall_dynamics/wall_dynamics_circles_k06.csv", f"{K}/wall_fourier_dyn_circles_k06.csv", "Per circle run, k = 6 (same script)"),
    ("7_flow_and_wall_dynamics/wall_fourier_circles_static.png", f"{F}/fourier_p1_static.png", "KEY: circle wall lattice: period, occupancy, order, pinning (+k and -k)"),
    ("7_flow_and_wall_dynamics/wall_fourier_circles_dynamic.png", f"{F}/fourier_p2_dynamic.png", "Lattice memory, rearrangement time, sliding speed (none)"),
    ("7_flow_and_wall_dynamics/wall_fourier_circles.pdf", f"{F}/fourier_maps.pdf", "Both pages, PDF"),
    ("7_flow_and_wall_dynamics/wall_fourier_circles.csv", f"{F}/wall_fourier_all.csv", "Per circle run"),
    ("7_flow_and_wall_dynamics/flow_mean_circulation.png", f"runs/cavity/analysis_flow/flow_p1_mean_flow.png", "KEY: mean flow and circulation collapse on R/l_a; breakdown above 2.5 l_a"),
    ("7_flow_and_wall_dynamics/flow_enstrophy.png", f"runs/cavity/analysis_flow/flow_p2_enstrophy.png", "Enstrophy vs R/l_a"),
    ("7_flow_and_wall_dynamics/wall_defect_tilt.png", f"runs/cavity/analysis_flow/flow_p3_tilt.png", "Wall +1/2 lean downstream 60-65 deg (30-40 at high rho)"),
    ("7_flow_and_wall_dynamics/flow.pdf", f"runs/cavity/analysis_flow/flow_maps.pdf", "All flow pages, PDF"),
    ("7_flow_and_wall_dynamics/flow_tilt.csv", f"runs/cavity/analysis_flow/flow_tilt.csv", "Per run"),
    ("7_flow_and_wall_dynamics/exchange_rates.txt", f"{W}/exchange_rates.txt", "Wall -> bulk escape rates: 10-20x faster in elongated ovals"),
    ("7_flow_and_wall_dynamics/wall_defect_labels_circles_ovals.png", f"{A}/bc_labels_k06.png", "Where the wall defects go: circles keep them, ovals lose them in 100-200 time units"),
    ("7_flow_and_wall_dynamics/wall_defect_labels_circles_ovals.txt", f"{A}/bc_labels_k06.txt", "Same, numbers"),
    ("7_flow_and_wall_dynamics/wall_defect_labels_small_ovals.png", f"{S}/bc_labels_P3.png", "Same for the small ovals"),
    ("7_flow_and_wall_dynamics/wall_defect_labels_lobe_oval.png", f"{K}/P4_lobe_ref/bc_labels.png", "Same for the lobe-shaped oval"),
    ("7_flow_and_wall_dynamics/wall_fourier_ovals/size_series_2.16la.txt", f"{K}/W1_size/wall_fourier.txt", "Static wall Fourier, ovals half-length 2.16 l_a"),
    ("7_flow_and_wall_dynamics/wall_fourier_ovals/size_series_2.16la_maps.png", f"{K}/W1_size/wall_fourier_maps.png", "Defect density maps"),
    ("7_flow_and_wall_dynamics/wall_fourier_ovals/k_scan.txt", f"{K}/W2_kscan/wall_fourier.txt", "k scan"),
    ("7_flow_and_wall_dynamics/wall_fourier_ovals/k_scan_maps.png", f"{K}/W2_kscan/wall_fourier_maps.png", "Defect density maps"),
    ("7_flow_and_wall_dynamics/wall_fourier_ovals/ellipses.txt", f"{K}/W3_ellipse/wall_fourier.txt", "Ellipses"),
    ("7_flow_and_wall_dynamics/wall_fourier_ovals/ellipses_maps.png", f"{K}/W3_ellipse/wall_fourier_maps.png", "Defect density maps"),
    ("7_flow_and_wall_dynamics/superseded/wall_fourier_first_test_maps.png", f"{W}/wall_fourier_test_maps.png", "SUPERSEDED: first test of the wall Fourier analysis"),
    # ---------------------------------------------------------------- 8 reports
    ("8_reports/project_report.pdf", "reports/project_report/report.pdf", "Full project report (21 pages, 5 Oct)"),
    ("8_reports/master_document.pdf", "reports/master/master.pdf", "MASTER document: paper skeleton with everything (methods in detail, circles, ovals, ellipses), grouped by shape"),
    ("8_reports/upstream_bug_report.md", "reports/upstream_bug_report.md", "Bug report for NIAMH-MPCD upstream (seed, tracker, Makefile)"),
    ("8_reports/oval_sweep_draft.md", "reports/oval_sweep_draft.md", "Planned oval sweep (not launched)"),
    ("8_reports/cassini_patch_plan.md", "reports/cassini_patch_plan.md", "Implementation plan of the Cassini walls"),
    ("8_reports/status_archive_2026-10-01.md", "reports/status_archive_2026-10-01.md", "Old status block (archive)"),
    ("8_reports/collaborator_brief/overcharged_brief.html", "reports/overcharged_brief.html", "Brief for Tyler (v7); open in a browser"),
] + [(f"8_reports/collaborator_brief/{f}", f"reports/{f}", "Figure of the brief") for f in
     ("fig_layer_area_all.png", "fig_layer_depth.png", "fig_layer_diagram.png", "fig_mirror.png", "fig_wall_dynamics.png")] + [
    ("8_reports/defect_geography/defect_geography.html", "reports/defect_geography.html", "Spatial-analysis write-up (30 Sep); open in a browser"),
] + [(f"8_reports/defect_geography/fig{i}", f"{A}/report/fig{i}", "Figure of the write-up") for i in
     ("1_maps.png", "2_radial.png", "3_sitedist.png", "4_bulkcharge.png", "5_g_same.png", "6_g_pm.png", "7_g_act.png")]

# movies: every oval movie and its last frame
for f in sorted(os.listdir(os.path.join(ROOT, K, "movies"))):
    p = os.path.join(K, "movies", f)
    if os.path.isfile(os.path.join(ROOT, p)):
        MAP.append((f"6_shape/movies/{f}", p, "Last frame" if f.endswith("_last.png") else "Movie (vorticity, director, defects)"))
MAP.append(("6_shape/movies/H34_ba2.0_k06_s01.mp4", f"{K}/movies/P3_s01/s01.mp4", "Movie, small oval"))
MAP.append(("6_shape/movies/H34_ba2.0_k06_s01_last.png", f"{K}/movies/P3_s01/s01_last.png", "Last frame"))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--dry-run", action="store_true"); a = ap.parse_args()
    dests = [d for d, _, _ in MAP]
    dup = {d for d in dests if dests.count(d) > 1}
    if dup: sys.exit(f"duplicate destinations: {dup}")
    missing, by_dir, nbytes = [], defaultdict(list), 0
    for dest, src, desc in MAP:
        s = os.path.join(ROOT, src)
        if not os.path.isfile(s): missing.append(src); continue
        by_dir[os.path.dirname(dest)].append((os.path.basename(dest), src, desc)); nbytes += os.path.getsize(s)
        if not a.dry_run:
            d = os.path.join(OUT, dest); os.makedirs(os.path.dirname(d), exist_ok=True); shutil.copy2(s, d)
    for m in missing: print("MISSING (skipped):", m)
    print(f"{len(MAP) - len(missing)} files, {nbytes / 1e6:.0f} MB{' (dry run)' if a.dry_run else ''}")
    if a.dry_run: return
    for d, rows in by_dir.items():
        top = d.split("/")[0]
        with open(os.path.join(OUT, d, "FILES.md"), "w") as f:
            f.write(f"# {d}\n\n{FOLDERS.get(top, '')}\n\nCopied by `tools/build_results.py`; regenerate rather than edit. "
                    "Files marked KEY are the main figures; SUPERSEDED files are kept for the record only.\n\n"
                    "| File | What it is | Copied from |\n|---|---|---|\n")
            for name, src, desc in rows:   # mapping order (KEY figures first in each block, PDF next to its PNG)
                f.write(f"| `{name}` | {desc} | `{src}` |\n")
    with open(os.path.join(OUT, "MANIFEST.tsv"), "w") as f:
        f.write("destination\tsource\tdescription\n")
        for dest, src, desc in MAP:
            if src not in missing: f.write(f"{dest}\t{src}\t{desc}\n")
    print("written", OUT)


if __name__ == "__main__":
    main()
