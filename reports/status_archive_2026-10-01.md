## ▶ Current status (updated 2026-09-30 ~15:15 office time) — read this first after any interruption

**Office machine (Ryzen 7950X):**

| What | State | Where / log | If the machine went down |
|---|---|---|---|
| ℓ_a sweep, 100² PBC, 9 activities | **done**; results in `runs/activeLength/L100/active_length_results.txt` | `runs/activeLength/L100/` | — |
| ℓ_a finite-size check, 200², act 0.05/0.07/0.09 | **done** 04:13; results in `runs/activeLength/L200/finite_size_results.txt` | `runs/activeLength/L200/` | — |
| Fixed-R scan R=31.83, k=2…10 (with flowOut/500) | **done** 19:07; results in `runs/cavity/Rscan_R31.83/rscan_results.txt` | `runs/cavity/Rscan_R31.83/` | — |
| Main sweep ρ ≥ 0.10 (426 runs) + passive reference set (142 runs) | **moved to the server** (see below) | — | — |
| Mirror sweep (customAnch = −k) | generator ready, not generated | `gen_sweep.py --sign -1` | — |

**Server (`baloewe@146.155.46.150`, EPYC 9454P; details in memory):**
- **Progress 2026-09-30 14:21:** batch 1 560/568 DONE (last: large ρ = 0.10, k ≈ 17–20); batches 2+3 started 00:31 as planned; 104 solvers, 370 GB free.
- **RUNNING — batch 1 (launched 2026-09-29 16:48 server time), 568 jobs, 96 in parallel**, queue `~/Overcharged/runs/cavity/jobs_batch1.txt`, log `queue_batch1.log`:
  - `main_act0.07/` — ρ = 0.10…0.50, k ≤ 20, R ≥ 5, 3 seeds, act 0.07, 10k warm-up + 100k steps; defectsOut/20, dirSOut/1000 and flowOut/1000 on all seeds, enstrophyOut/100 on seed 1.
  - `main_passive/` — same geometries, act 0, 1 seed, 10k + 30k steps, defectsOut/20, dirSOut/1000 (layer reference, H6).
  - Finished runs get a `DONE` marker and gzip-compressed outputs. Expected finish ~1–1.5 days.
  - If the server went down: `cd ~/Overcharged/runs/cavity && setsid -f nohup ./run_queue.sh jobs_batch1.txt 96 > queue_batch1.log 2>&1 < /dev/null` (skips DONE runs; unfinished ones restart from scratch).
- **QUEUED — batches 2 + 3** via `queue_chain.sh` (log `queue_chain.log`): starts one pool of 96 (`jobs_batch23.txt`, log `queue_batch23.log`) as soon as batch 1 has *started* all its jobs; batch 2 first, then batch 3.
  - **Batch 2** (544 jobs, ~35 h): positive winding, new ρ = 0.03–0.09, 0.12, 0.14, 0.16 (R ≤ 64, 3 seeds) into `main_act0.07/` + passive set into `main_passive/` (`manifest_batch2.csv`).
  - **Batch 3** (760 jobs, ~44 h): mirror, customAnch = −k, ρ = 0.03–0.10, 0.12, 0.14, 0.15, 0.16, 0.20 (R ≤ 64, 3 seeds) → `main_act0.07_neg/` (dirs `kmNN`) + passive `main_passive_neg/` (`manifest_batch3.csv`).
  - Same outputs as batch 1 (defectsOut/20, dirSOut/1000 + flowOut/1000 all seeds, enstrophy seed 1). Batch 1's ρ ≥ 0.25 part is beyond the resolution limit (methods figure only).
  - **Expected: all three batches done ≈ 2026-10-03.** Then the firmware/BIOS (CPU clock) work with Juan (a reboot kills running jobs).
  - If the server went down: relaunch batch 1 as above, then `setsid -f nohup ./queue_chain.sh > queue_chain.log 2>&1 < /dev/null` (finished runs are skipped everywhere).
- `~/Overcharged` synced (code + scripts + inputs) and built; venv `~/Overcharged/.venv` with numpy/scipy/matplotlib.
- **Benchmark done** (`runs/server_benchmark.txt`): total throughput 7.7 / 12.9 / 19.8 / **28.7 MPUPS** at 16 / 32 / 48 / **96** runs at once (office ≈ 6.1 with 9 runs). Use up to 96 at once, largest first; each run is ~4.4× slower than alone on the office machine, so the largest runs set a sweep's finish time.
- **Clock (checked 2026-09-30 14:25, 104 solvers running):** `boost` = 1 and `cpb` = 1 (Juan is right: boost is on). acpi-cpufreq only lists P-states up to 2.75 GHz (`bios_limit` 2.75 GHz); boost above that is handled by the hardware and doesn't show in `scaling_max_freq`. APERF/MPERF (`scaling_cur_freq`) reads **2.945 GHz on all 96 CPUs** under full load, i.e. slightly above base but well below the 3.8 GHz single-core max (all-core spec ~3.65). Yesterday's "2.35 GHz" reading doesn't match this and was probably measured wrongly. Governor is back to `schedutil` on all CPUs. Possible gain from amd-pstate/BIOS power settings (PPT/cTDP): at most ~2.95 → ~3.5 GHz (~15–20%), not the 1.4× estimated earlier. Low priority; nothing to do before the sweeps end.
- **Seeds:** the server binary still has the seed bug (TODO 🔴): batches 1–3 run with effectively random seeds. Results valid; exact reruns impossible. Rebuild the server binary after the sweeps.
- Storage: SSD only (387 GB free). HDD `/opt/disco2` is off-limits (unreliable). Compress finished runs.

Next steps for Claude:
1. ~~L100 analysis~~ done 18:16.
1b. ~~Server benchmark~~ done → `runs/server_benchmark.txt` (96-way best). After the driver fix: rerun `bench.sh` and compare the 96-way point.
2. ~~R-scan analysis~~ done 19:20 → H2 not supported at R/ℓ_a = 1.5; main sweep is the next step (H4 needs R/ℓ_a variation). Before analysing the sweep: data-driven edge/bulk cut + error bars in `circulation.py` (TODO).
3. ~~L200 finite-size check~~ done 04:14 → no significant finite-size effect; reference ℓ_a(0.07) = 22.2 ± 0.3.
4b. ~~Lunch tasks~~ **done** 19:17 on 09-29, analysed 2026-09-30 → `runs/cavity/analysis_lunch/` (see results log). Were: R-scan seeds 2–3 (`runs/cavity/Rscan_R31.83_seeds/`), H3 pilot at act 0.05/0.09 with R = 1.434 ℓ_a (`runs/cavity/H3pilot/`), passive references for the R-scan geometry (`runs/cavity/Rscan_R31.83_passive/`); queues `queue_office_lunch.log`, `queue_office_passive.log` (legacy tracker, fixed seeds). ETA ~3–4 h. Then: `layer_analysis.py` (with `--passive-root Rscan_R31.83_passive`) + `circulation.py` on all three.
5. ~~Reply to Juan~~ done by the user 2026-09-30 (boost is on, clock ~2.95 GHz all-core; see Clock line). **Waiting on the user:** `/home/hcorvalan` / office HDDs later.
4. ~~Defect-tracker rewrite~~ done 2026-09-29 (office build; see `TODO.md`).
4c. ~~act 0.05 seeds 2–3~~ **done** 13:46, analysed (results log 14:10). Were: act 0.05 seeds 2–3 (18 runs, 18 parallel; `runs/cavity/H3pilot_seeds/`, `queue_h3_seeds.log`). ETA ~2–3 h. Then `layer_analysis.py` + `circulation.py` + pair counts on all 3 seeds.
4d. ~~H3 test at ζ = 0.06 and 0.11~~ **done** 2026-10-01 12:52, analysed (results log 15:00). Was: chain, chain `runs/cavity/chain_h3_high.sh` (log `chain_h3_high.log`). Stage 1: ℓ_a in 200² PBC boxes, 2 seeds each (`runs/activeLength/L200/act_0.060_s{1,2}`, `act_0.110_s{1,2}`, legacy tracker, seeds 20601/20602/21101/21102), ~9–12 h. Stage 2 (automatic): ℓ_a seed mean → `runs/cavity/H3high/l_a_used.txt`, cavities R = 1.434 ℓ_a, k = 2…10 × 3 seeds (`H3high/act0.06|0.11/sNN/kNN`, seeds 60000+100s+k / 110000+100s+k, inputs as `H3pilot_seeds` otherwise). Stage 3: 54 cavities, 18 parallel (`queue_h3_high.log`), ~3–4 h. Expected done ≈ 2026-10-01 morning–noon. If the machine went down: relaunch the chain (finished stage-1 runs and DONE cavities are skipped; unfinished runs restart).
4e. **Running on the office machine since 2026-10-01 12:31: first active Cassini pilot (H8)** — `runs/cassini/P1/H{48,64}_ba{1.6,1.3,1.1}_k06` (half-length 48/64 = 2.2/2.9 ℓ_a; b/a 1.6 convex control, 1.3 and 1.1 pinched), k = 6, ζ = 0.07, eMode 3, 1 seed, 10k + 40k steps, defects/20, director+flow/500. Binary `runs/cassini/bin/mpcd_cassini_e63155b.out` (branch `cassini`), queue `run_queue_cassini.sh` → `queue_pilot1.log`. ETA ~2 h (H48) to ~5 h (H64, 178k particles).
4h. ~~Density check~~ **done** 15:41 (results log 15:45). Was: (`runs/cavity/density_check/act{0.07,0.09,0.11}_k08`, k = 8 at the H3 radii, 10k + 20k steps, `densOut` 500, `queue_density_check.log`). Question (user, 2026-10-01): at high ζ density fluctuations / near-empty cells may make the runs unusable. Empty-cell proxy from the director output (S = 0, zero director; `analysis_spatial/empty_cells.txt`): none at ζ ≤ 0.06, rare at 0.07, 1–5×10⁻⁵ of interior cells at 0.09, 6–8×10⁻⁵ at 0.11 (11–16% of snapshots have one) — all within a few cells of the cavity centre (vortex core), none near the wall. densOut will give P(n per cell) and the radial density profile.
- **Working conclusion (user, 2026-10-01): ζ = 0.07 is the sweet spot** (clean +1 plateau, deepest gap, no onset effects, negligible empty cells). High-ζ end possibly unusable if occupancy is poor; then test the high side at 0.09–0.10 instead of 0.11.
4g. **Running on the office machine since 2026-10-01 14:49: ζ = 0.09 seeds 2–3** (`runs/cavity/H3pilot_seeds/act0.09/s0{2,3}/`, 18 runs, 10 parallel, `queue_h3_009_seeds.log`; seeds 91500+k / 92500+k, inputs as `H3pilot/act0.09` otherwise, production binary). ETA ~2–2.5 h. Then: `edge_sites.py`, `circulation.py`, enclosed charge, gap profile on all 3 seeds.
4i. **Running since 2026-10-01 16:25: Cassini P3** (`runs/cassini/P3`, `queue_pilot3.log`, 12 parallel): (a) high k on H48 b/a 1.02 (perimeter 240): passive checks k = 18, 20 first, then active k = 20, 18, 15, 12 (10k + 40k steps, 1 seed). Limit: imposed density 2k/P ≤ 0.15 ⇒ k ≤ 18 clean, k = 20 (0.167) marginal, k ≥ 25 fails (needs a larger oval, P ≥ 400 for k = 30). (b) **smaller ovals for the size control of H1**: half-length 34 (1.53 ℓ_a, where circles are clean), b/a = 2.5, 2.0, 1.6, 1.3, k = 6, 3 seeds, 10k + 100k steps, defects/20, director + flow/1000 (as the server circles). ETA: passive minutes, high k ~1.5 h, small ovals ~evening.
4f. **Running on the office machine since 2026-10-01 14:19: Cassini pilot P2** (`runs/cassini/P2`, `queue_pilot2.log`): active k = 6, ζ = 0.07, H64 b/a = 1.05/1.02/1.01 and H48 b/a = 1.02 (waists 28/18/13/13.5 cells). Passive neck checks at H64 1.01 and H48 1.02 passed (alignment 0.93, audit 100%, 12 +½). P1 H64 runs still going.
6a. **Edge/bulk:** pinned-site rule (`edge_sites.py`) now preferred over d_ρ (no radial cut, shape-independent); use it in `circulation.py` / main-sweep analysis.
6. **Next:** data-driven edge cut (d_ρ, as in `layer_analysis.py`) in `circulation.py` — its fixed `--wall 2.5` cut gives ⟨q_b⟩ ≈ 1.1–1.5, P(q_b=1) ≈ 0.2–0.8 on the same runs where the density cut gives 1.00 / 0.9. Then an H3 follow-up (act 0.05, 3 seeds, maybe act 0.04/0.06).

---

