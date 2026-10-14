# TODO — Overcharged / NIAMH-MPCD

Urgency: 🔴 **High** (can bias current results) · 🟠 **Medium** (needed before a planned step) · 🟢 **Low** (nice to have)

Paths are relative to `mpcd/`. Research questions, planned runs and results: see `RESEARCH.md`.

---

## Physics / measurement correctness

### ✅ Input seed was silently ignored → no run was reproducible (FIXED: office build 2026-09-29; server production binary 2026-10-05)
- `readJson()` calls `accessibleVolume()` (Monte Carlo, `read.c` ~l.1619/1721) before `initvar()` seeds the RNG (`init.c:1183`). The first draw made xoshiro self-seed from time+pid (`X_genrand_int32`), and `X_RandomSeedSRD()` then skipped the input seed because the generator was "already seeded", while returning the input seed → `synopsis.dat` reported a seed that was never used. Same input + same seed gave different runs from the first frame.
- Fix (office): seed right after `seed` is read in `readJson()` (`in->seed = RandomSeedSRD(in->seed)`, so the synopsis reports the seed actually used, also the time-based one for seed 0), and `X_RandomSeedSRD()` always (re)initialises. Verified: identical runs for identical seeds; the tracker choice doesn't change the dynamics.
- **Consequences:** all runs before the fix (ℓ_a sweeps, validations, R-scan, server batches 1–3) have unknown, effectively random seeds. They're statistically independent, so results stand, but they can't be regenerated exactly, and the seeds in their inputs/synopses are meaningless. State this in the methods if reproducibility is claimed; new runs should use the fixed binary.
- **Server done 2026-10-05:** `~/Overcharged/mpcd` source replaced by repo `main` (fa68326; md5 of every .c/.h = main), rebuilt; old source + binary kept in `~/Overcharged/mpcd_legacy_2026-09-28/`. Server diff was only additive (seed fix, union-find tracker + `defectTracker` switch, eMode 2 branch): circle physics unchanged. Test `~/Overcharged/runs/seedfix_test_2026-10-05/` (k = 6, ρ = 0.06, 2.2k steps): same seed → flowfield, directorfield, defects, enstrophy byte-identical; other seed → different; `defectTracker` 0 vs 1 → byte-identical outputs including defects. The cassini binary (`~/Overcharged/niamh-overcharged-cassini`) already had the fix.
- Remaining: report upstream (affects every NIAMH-MPCD user with a fixed seed). Check whether checkpoint recovery restores the RNG state or re-seeds (`readchckpnt`).

### ✅ Defect clustering can over-count defects — DONE 2026-09-29 (office build; server binary still legacy)
- **Done:** union-find tracker is the default (`"defectTracker": 1`), legacy kept as `topoChargeAndDefectsOut_legacy()` (`"defectTracker": 0`), choice + tolD in `synopsis.dat`, documented in `InputGuide.md`. Validated: == independent Python re-implementation in 101/101 frames (100² PBC, act 0.07); PBC net charge 0 in 98.0% of frames (legacy 64.4%); cavity defect lists identical to legacy in 251/251 frames. Upstream patch `reports/upstream_patches/fixes/0002-*`. All data so far (ℓ_a, R-scan, lunch runs, server batches 1–3) is legacy-tracker; fine for cavities (identical lists), ℓ_a uses the (L−2)² effective area. Remaining: server rebuild after batches 1–3.
- `mpcd/subroutines/pout.c`, `topoChargeAndDefectsOut()` (~l.1860–1910).
- Charged cells are grouped into defects in a single raster pass that only looks at already-visited neighbours. A cluster with a U-shape (or similar) gets split into two defects → N too high → active length too small.
- "First column" branch compares `topoC[i][j]*topoC[i][j]` (always > 0) instead of the neighbour's charge, so opposite-sign cells there can be merged.
- Options: proper connected-component labelling (union-find) in C, or post-merge same-sign defects closer than ~1.5 cells in `runs/activeLength/active_length.py`.
- Affects the active-length measurement → check before quoting ℓ_a. (First synthetic test: ℓ_a changes by ~2% only.)
- **Planned after the L100 active-length sweep finishes**, together with the PBC item below:
  - New tracker: union-find connected components (8 neighbours), periodic wrapping when `XYZPBC`, seam-aware centre of mass.
  - Keep the old tracker as a separate function `topoChargeAndDefectsOut_legacy()`, selected by a JSON key (e.g. `"defectTracker": 1` new [default] / `0` legacy). Record which tracker was used in `synopsis.dat`. Document it in `InputGuide.md`.
  - Validate against an independent Python charge calculation on `dirSOut` output, frame by frame (including defects crossing the seam, and net charge = 0 with PBC).

### 🟠 Ellipse and channel support in the analysis and anchoring (before those runs)
- ✅ **Arclength mode `eMode=2` implemented (office build, 2026-09-29)** for circles/ellipses: `anchArclengthFraction()` + branch in `oriBC()`. Validated (`runs/emode2Validation/`): circle identical to eMode 0 (alignment 0.984/0.985, audit 100%); ellipse a=20, b=10, k=4: tip winding fraction 0.384 vs arclength fraction 0.396 (eMode 1: 0.478 vs 0.500; eMode 0: 0.694 vs 0.709, and eMode 0 fails the audit in 14% of frames at its tip density ≈ 0.25). Remaining: flat walls (channel, 🟢), rotated walls (O[2] item below applies to eMode 1 and 2 as well).
- (original plan) - **Arclength-based anchoring mode** (e.g. `eMode=2`) in `oriBC()` (`lc.c`). **Effectively required**, not optional: the local resolution limit ρ(s) ≲ 0.15 makes eMode 0 unusable beyond a/b ≈ 2 and tightens eMode 1 with eccentricity; arclength-uniform keeps ρ(s) uniform (see `RESEARCH.md`, ellipse notes). anchoring angle ψ = ψ₀ + 2π·k·s/P, with s the arclength along the wall at the collision point and P the period (perimeter). One implementation for all shapes:
  - circle: s = Rφ, P = 2πR (same as the current hack);
  - ellipse: s(t) from an inverse-arclength lookup table per wall (elliptic integral), built at init. With `eMode=1` (parameter-uniform) the imposed density is a/b higher at the tips (see `RESEARCH.md`, ellipse notes);
  - 🟢 *(channel, only if needed: paper 1 is circle + ellipse; channel possibly for the referee reply)* **plane (channel walls)**: s = coordinate along the wall (y for vertical walls), P = L_y (domain period). The current normal-angle hack gives **no winding on a flat wall** (constant normal). Design the arclength mode so this is a small addition later.
  - Note: the anchoring must then use the particle position `pp->Q` (in the wall frame after `rotateBC`), not only the normal `n`.
  - Prototype on the circle first: it must reproduce the current `customAnch` results exactly.
- 🟢 *(channel, only if needed)* **Channel inputs:** two planar walls (P=[1,1,1], aInv with only the x component) at x = 0 and x = W, PBC in y. `domainWalls` / BC setup to give PBC in y only. Winding per wall per period ∈ ℤ/2; opposite sense first (see `RESEARCH.md`).
- **Analysis for non-circular walls:** `cavity_check.py` / `circulation.py` assume a circle (radius from `C`, circles for Q(r), v_φ about the centre). Generalise: wall-distance contours (or level sets of W) for Q and the edge/bulk split, bulk = distance from wall > d. The charge audit (`runs/cavity/defect_audit.py`) already works for any shape.
- **Generator:** `gen_cavity.py` for ellipses (aInv semi-axes, `eMode`, box size from a and b, population = dens·π·a·b, check b²/a ≥ 5).
- Settle the `O[2]` item below before any rotated ellipse.

### 🟠 Rotated ellipses: `customAnch` + `eMode=1` ignores `O[2]` consistently
- `mpcd/subroutines/lc.c`, `oriBC()`.
- Circle branch subtracts `O[2]` from the normal angle and adds it back; ellipse branch subtracts it but never rotates `nO` back. Also, `rotateBC()` may already put `n` in the wall frame, so the subtraction itself may be double-counting.
- Harmless for `O[2]=0`. Must be resolved before using rotated ellipses.

### ✅ Defect/charge detection ignores the outer cell layer and doesn't wrap PBCs — DONE 2026-09-29 (new tracker, see above)
- `pout.c`, `topoChargeAndDefectsOut()`: charge computed only for cells `1..L-2` (angle only for `2..L-3`); clustering doesn't wrap across the periodic seam.
- The mean count is fine if you use the effective area `(L-2)²` (the analysis does). But any defect in the strip is lost, so frames aren't charge-neutral. With L=100 the strip is ~4% of the area, so with ~40 defects ~80% of frames miss at least one. That makes charge-conservation filtering of frames useless, and it corrupts defect tracking and pair statistics near the seam.
- Fix: when the domain walls are periodic, compute `topoChargeLocal()`/`topoAngleLocal()` with wrapped neighbour indices on all cells, and make the clustering periodic. Keep the current behaviour for solid walls.
- After this fix, net charge ≠ 0 in a PBC frame flags only genuine clustering errors → becomes a clean diagnostic.

---

### 🟢 Cassini-oval walls + eMode 3 (for H8, after paper 1)
- **Implemented 2026-10-01** on branch `cassini` (e63155b), passive validation passed; active pilot running. Not merged into `mpcd/`.
- 🟠 Likely cause found 2026-10-05 (not yet tested): the checkpoint writer (`pout.c` ~l.2500) prints **31 ints + 1 float** on the output-flags line (order …FLOWOUT, SWFLOWOUT, VELOUT, DENSITYOUT…), the reader (`read.c` l.727) scans **32 ints + 1 float** in a different order (…FLOWOUT, VELOUT, DENSITYOUT, SWFLOWOUT, DENSITYOUT…) → every later field is shifted on recovery. Also `outFlag.CHCKPNTrcvr` is never initialised (mpcd.c uses a local variable), so checkpoints store garbage there (the only byte difference between two same-seed checkpoints).
- 🟠 Found along the way: **checkpoint restarts segfault** (input `{"checkpointIn": ...}`), in the production build too. Not used by our workflow (run_queue restarts from scratch), but needed if long runs ever have to be resumed.
- Plan: `reports/cassini_patch_plan.md`. Do it on a `cassini` branch of `niamh-overcharged` (private GitHub repo), not in `mpcd/`; merge only after validation and never mid-sweep.

---

## Performance

### 🟠 Sort particles by cell for cache locality
- Cost per particle update grows with system size because cell linked lists jump randomly through the particle array (128 B/particle; exceeds L3 above ~200k particles).
- Measured (sample `activeNem2D` input, single run, Ryzen 9 7950X): 50² → 2.0 MPUPS, 100² → 1.3 MPUPS, 200² → 0.63 MPUPS.
- Fix: every N steps, counting-sort `SRDparticles` by cell index and rebuild the lists. Make it opt-in (e.g. `"sortEvery": 50`). Caveat: reorders per-particle outputs.
- Validate: statistically identical results on the 50² sample; benchmark 200².

### 🟢 Redundant end-of-step `localPROP()`
- `mpcd/subroutines/mpc.c:~4987`. Full per-cell pass (moment of inertia + order tensor + eigensolve) re-run after re-binning; the next step recomputes it before any collision. ~15–18% of runtime.
- Guard with "output due this step || MD/swimmers active". Verify nothing else reads these cell values first.

### 🟢 OpenMP parallelisation of cell loops
- Collision and `localPROP` loops are embarrassingly parallel except for the global Mersenne Twister RNG → needs per-thread RNG streams. Breaks bit-reproducibility with the serial code.

---

## Bugs (reporting / tooling)

### 🟢 PUPS integer overflow
- `mpcd/mpcd.c:341–342`: `int mpcUpdates = (simSteps+warmupSteps) * GPOP` overflows for large runs (e.g. 25000 × 800k). Use `double`/`long long`. Cosmetic only.

---

## Workflow / tooling

### 🟠 Cavity analysis: data-driven edge/bulk cut and error bars (before analysing the main sweep)
- `runs/cavity/circulation.py` uses a fixed bulk cut r < R − 2.5. The R-scan showed this is unsafe: the wall layer is diffuse and ~4 cells thick at low k, sharp at high k → a fixed d fakes a ρ-dependence (see `runs/cavity/Rscan_R31.83/rscan_results.txt`).
- Add: per-run edge/bulk division from the minimum of the radial +½ density (and/or the Q(r) plateau near 1); report bulk charge vs d as a sensitivity check; keep a fixed-d option for comparison.
- Add statistical errors: block averages over time per run (P(q_b = 1), ⟨q_b⟩, Φ), then seed-to-seed errors per parameter set.
- Also add radial ±½ density profiles and Q(r) (from `dirSOut`) to the per-set outputs.

### 🟠 Cloud backup of the server data (once the server runs finish; added 2026-10-11)
- Office machine already backed up nightly: `backup/backup_overcharged.sh` (rclone copy → `onedrive:Backups/Overcharged/office`, cron 03:00; see RESEARCH.md 2026-10-11). Server `~/Overcharged` (~63 GB, ~43k files on 2026-10-11; growing with seeds s04–s10 and the mirror sweep) has no off-site copy yet.
- When: after the circle seeds (`run_queue.sh jobs_seeds10.txt`) and the chained mirror shape sweep finish, so the first upload is complete and doesn't compete with production I/O.
- Steps: (1) copy `~/bin/rclone` (v1.75.2) to the server's `~/bin`; (2) own token for the server (don't share the office token): user runs `ssh -N -L 53682:localhost:53682 baloewe@146.155.46.150` on the laptop, then `~/bin/rclone config` on the server (remote `onedrive`, OneDrive Business, drive "OneDrive (business)", auto-config y); (3) adapt the script → `onedrive:Backups/Overcharged/server`, same filters (skip `deleteme*`, live run dirs); SSD only, never `~/disco2`; (4) dry run, first upload, then a nightly cron offset from the office one (e.g. 04:00).
- Check afterwards: `rclone check` file counts for `runs/cavity` and `runs/shapes`.

### 🟠 `run_queue.sh` stop instructions are wrong (found 2026-10-05)
- The header says `pkill -f run_queue.sh; then kill the running solvers`. Killing the script leaves its `xargs -P` pool alive, which starts the next job each time a solver is killed (192 runs started and killed on 5 Oct). Correct: kill the xargs PID first, then `pkill -u baloewe -x mpcd.out`. Fix: trap SIGTERM in `run_queue.sh` to kill its xargs child, or write the xargs PID to a file; update the header in `runs/cavity/run_queue.sh`, `runs/cassini/run_queue_cassini.sh` and the server copies.

### 🟢 Resume interrupted runs from checkpoints in `run_queue.sh`
- `runs/cavity/run_queue.sh` currently restarts an unfinished run (no "Simulation successful" in `log.txt`) from scratch. After a power loss, the big runs (up to ~7 h for ρ ≥ 0.10, ~50 h for ρ = 0.05) would lose all progress.
- All runs write `checkpoint.dat` every `checkpointOut` steps (10k in the sweep). Per `InputGuide.md`, setting `seed: -1` loads the checkpoint.
- To do: check how the solver finds the checkpoint (path, input vs output dir), whether it keeps the warm-up/production step count, and whether output files are appended or overwritten on resume (otherwise `defects.dat` data before the checkpoint could be lost or duplicated). Then make `run_one` detect a partial run with a checkpoint and resume it instead of restarting.
- Test: kill a small run mid-way, resume, and check the defect time series is continuous with no gaps or duplicate times.
- Also useful to extend runs that haven't reached stationarity (the ℓ_a sweeps already write checkpoints for this).

### 🟢 Scalar circulation output (angular momentum about the cavity centre)
- H4 (see `RESEARCH.md`) needs ⟨v_φ⟩ over space as a time series. It currently comes from full `flowfield.dat` frames: ~17 GB for the ρ ≥ 0.10 sweep, and only at the `flowOut` interval (each frame is an interval average).
- Add an output flag (e.g. `circOut`) that writes per step or per N steps: t, ⟨v_φ⟩, ⟨|v|⟩, L_z about a given centre (first BC's `Q`, or a key). Tiny files, full time resolution → can resolve vortex reversals and correlate them with bulk-charge switches.
- Check that `avVelOut` isn't a substitute: it's the net velocity (≈ 0 for a vortex).

### 🟠 Makefile doesn't track header dependencies
- `mpcd/Makefile` rebuilds `%.o` only when its `.c` changes. Changing a header (e.g. adding a field to a struct in `SRDclss.h`) leaves stale objects compiled against the old layout → silent memory corruption. Hit on 2026-09-29 (tracker rewrite); fixed by `make clean && make`.
- Fix: add automatic dependency files (`-MMD -MP` and `-include $(objects:.o=.d)`), or make every object depend on `$(headers)`. Until then: **always `make clean && make` after touching a header** (office and server).

### 🟢 `topoSmallestAngle()` modifies its argument
- `lc.c`: flips `v` in place; the legacy `topoChargeLocal()` passes `CL[..].DIR` directly, so computing charges flips stored cell directors (n → −n). Physically harmless for a nematic, and `localPROP()` recomputes `DIR` every step before use. The new tracker passes copies. Consider making it take `const` inputs and work on a local copy.

### 🟢 Unknown code change shifted the activity scale (user decision 2026-09-29: don't track now)
- Since Kozhukhov, Loewe & Shendruk, Commun. Phys. (2024) ("Mitigating density fluctuations in particle-based active nematic simulations"), an unknown edit shifted the activity needed for active turbulence to higher values (paper: onset α ≈ 0.02; now even 0.05 is only mildly turbulent). Cause and date unknown.
- Consequence: bare activity values can't be compared with the 2024 paper (or older runs). Mitigation in this project: normalise by the **measured** ℓ_a, and set compressibility cutoffs from measured density fluctuations / Mach number, not from the paper's α†.
- If it's ever tracked: bisect upstream git history (server clone `~/NIAMH-MPCD`) with a fixed input (e.g. `activeNem2D`), comparing defect counts.

---

## Documentation

### 🟢 `sigWidth`/`sigPos` docs say collOp #20 only
- `InputGuide.md` (species table): "specifically for collision operator #20", but `lc.c` (~l.1877) applies the sigmoid falloff to both `DIPOLE_DIR_SIG` (20) and `DIPOLE_DIR_SIG_SUM` (21, the one this project uses). Note there too that at the mean density the falloff factor is ½(1 + tanh 1) ≈ 0.88 with the defaults.

- Mention: winding = `customAnch`, φ is the wall-normal angle, `eMode=1` uses the ellipse parameter (see `Hack_correction.nb`).

---

## Done
- ✅ Documented `customAnch` and `eMode` (0, 1, 2) in `InputGuide.md` (2026-09-29).
- ✅ `oriBC()`: custom-anchoring block now only runs when `customAnch ≠ 0` (previously overrode normal anchoring on every wall and gave NaN on `aInv=[0,1,0]` planes).
- ✅ `oriBC()`: ellipse factor corrected from `b²/a²` to `b/a` (matches `Hack_correction.nb`).

- **(medium) Stale objects when switching branches in niamh-overcharged/** (found 2026-10-07): the Makefile has no header dependencies, so after `git checkout` of another branch `make` relinks old .o files (a rebuild of unmodified `main` differed from production). Always `make clean && make` after switching branches. With clean builds, production = main = anchmod branch (anchAmp 0) bit for bit.
- **(low) Queue scripts waiting on `pgrep -f <pattern>`** (found 2026-10-08): the pattern also matches any shell whose command line contains it (e.g. my own monitoring commands), so the halfk queue waited ~1 h after the anchmod runs had ended. Use `ps -eo args | grep -q "^<full binary path>"` (anchors on the executable) instead.
- **(medium, theory) Anisotropic Hessian of the passive wall ring — is "K3 fits the fluctuations" earned?** (user, 2026-10-08) The one-constant Hessian needed K_eff = 17.8 (factor 1.35); exact MC with grid rounding matches the along-wall rms with K = K3 = 23 (2–4%, N = 8–20; N = 4 overshoots 2.65 vs 2.35). The master text ("tangential displacements cost bend", §29.6) is a post-hoc reading. Test: anisotropic Hessian of the passive 2k ring with K1 = 42, K3 = 23.1, no fit (reuse reports/master/theory/aniso2/relax_frank2.py: rigid cores, full splay/bend energy; finite-difference second derivatives of the relaxed energy w.r.t. defect positions), plus grid rounding → compare with fluct_measured.json. If it fails, reword §29.6 and treat the K3 match as coincidental.


- [REVIEW with Benjamin, 2026-10-14] Shielding/vacancy maps: missing wells near some layer defects, circle asymmetry, very well shielded ellipse tips. Details and candidate checks: RESEARCH.md "TO REVIEW WITH BENJAMIN". PRL text using them is flagged.
