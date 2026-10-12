# NIAMH-MPCD: issues found while validating cavity simulations

*Benjamin Loewe, September 2026. Checked against `master` at **a2cf7b2** (2026-08-03) and at 2a38960 (2025-10-01); every issue below is present in both.*

While setting up active-nematic simulations in cavities with custom wall anchoring, we validated the solver more closely than usual (fixed-seed reproducibility, charge conservation of the defect tracker, performance scaling). This turned up several issues in the core code. Three have tested fixes ready as patches (one commit each, on top of current `master`); the rest are suggestions.

| # | Issue | Severity | Status |
|---|---|---|---|
| 1 | Input `seed` is silently ignored; no run is reproducible | **High** | patch 0001 |
| 2 | Defect tracker: loses the outer cell layer, no periodic wrapping, can split defects | **Medium** | patch 0002 (legacy tracker kept as an option) |
| 3 | Makefile doesn't track header dependencies | **Medium** | patch 0003 |
| 4 | Performance drops ~3× per particle for large systems (cache locality) | Low | suggestion |
| 5 | PUPS report: integer overflow, warm-up not counted | Low | suggestion |
| 6 | Legacy `.inp` input leaves `tolD` uninitialised | Low | fixed in patch 0002 |
| 7 | Docs: `sigWidth`/`sigPos` said to apply to collision operator 20 only | Low | suggestion |
| 8 | Activity scale has shifted since the 2024 density paper | Observation | not investigated |

---

## 1. Input `seed` is silently ignored (high)

**What happens.** Two runs with the same input and the same fixed `seed` differ from the first output frame. `synopsis.dat` reports the input seed ("Random seed: 4242"), but that seed was never used.

**Why.** `readJson()` calls `accessibleVolume()` (Monte Carlo; `read.c:1612` and `read.c:1714`) before `initvar()` applies the seed (`init.c:1183`). The first random draw finds xoshiro unseeded and seeds it from time + pid (`X_genrand_int32()`, `rand.c:258`). When `initvar()` later calls `X_RandomSeedSRD(seed)`, the `if (X_seeded == 0)` guard (`rand.c:239`) skips the input seed, and the function still returns it.

**Impact.** Every run since this ordering was introduced has an effectively random seed. Runs are statistically independent, so results stay valid, but they can't be regenerated exactly, and the seeds recorded in inputs and synopses are not the ones used.

**Fix (patch 0001).** Seed immediately after `seed` is read in `readJson()` (so `in->seed` holds the seed actually used, including the time-based one for `seed: 0`), and let `X_RandomSeedSRD()` always initialise, as `MT_RandomSeedSRD()` already does. Tested: identical input and seed now give byte-identical output; `seed: 0` reports the seed it drew, so such runs can be reproduced afterwards.

**Open point.** Checkpoint recovery (`readchckpnt`) restores the seed but, as far as we can see, not the generator state; worth checking.

## 2. Defect tracker (medium)

**What happens.** In periodic boxes, many frames aren't charge-neutral: 36% of frames in a 100×100 active box at activity 0.07 (about 25 defects per frame) have net charge ≠ 0, and the fraction grows with the number of defects (roughly 1 − 0.96ᴺ for N defects).

**Why** (`topoChargeAndDefectsOut()`, `pout.c:1829`):
- charge is computed only on cells `1..L−2` (`pout.c:1840`) and the orientation on `2..L−3`, with no periodic wrapping, so defects in the outer strip (~4% of a 100×100 box) are lost;
- clustering is a single raster pass that only looks at already-visited neighbours, so a U-shaped cluster can be split into two defects;
- the "first column" branch compares `topoC[i][j]*topoC[i][j]` (`pout.c:1882`, `1884`), always positive, instead of the neighbour's charge;
- `topoSmallestAngle()` flips its second argument in place (`lc.c:2401`) and is called with `CL[..].DIR`, so computing charges flips stored cell directors. Harmless for a nematic, since `localPROP()` recomputes `DIR` before use, but a side effect.

**Fix (patch 0002).** A new default tracker (`"defectTracker": 1`): charge and orientation wrap periodic directions (`XYZPBC`) and work on copies of the directors; defects are connected components (union-find over 8 neighbours, same sign, periodic) of cells with ||q| − ½| < `tolD`; the centre of mass is correct across the periodic seam. Output format and ordering are unchanged. The old tracker is kept, unmodified, as `topoChargeAndDefectsOut_legacy()` (`"defectTracker": 0`). The choice and `tolD` are written to `synopsis.dat` and documented in `InputGuide.md`.

**Evidence** (same seed with both trackers, after fix 1):

| check | 100×100 periodic box, activity 0.07 | circular cavity, k = 5 |
|---|---|---|
| dynamics with either tracker | byte-identical | identical (printed directors differ only by the legacy sign flip) |
| charges on interior cells, new vs legacy | identical | identical |
| new tracker vs an independent Python re-implementation (from `directorfield.dat`) | 101/101 frames agree | — |
| frames with the correct net charge | **98.0%** (legacy 64.4%) | 99.6% (both) |
| defects at the periodic boundary | 166 found (legacy: none) | — |
| defect lists, new vs legacy | — | identical in 251/251 frames |

So with walls the new tracker reproduces the legacy one exactly; with periodic boundaries it recovers the missing strip.

## 3. Makefile doesn't track header dependencies (medium)

**What happens.** `%.o: %.c` rebuilds an object only when its `.c` changes (`Makefile:33`). Changing a header, for example adding a field to a struct in `SRDclss.h`, leaves objects compiled against the old struct layout: silent memory corruption unless `make clean` is run. We hit this during the tracker change.

**Fix (patch 0003).** Compile with `-MMD -MP`, include the generated `.d` files, and remove them in `make clean`. Tested: touching `SRDclss.h` recompiles all 23 dependent files; a second `make` recompiles none.

## 4. Performance at large system sizes (suggestion)

Cost per particle update grows with system size (`activeNem2D` sample input, single run, Ryzen 9 7950X):

| box | particles | per-particle memory | MPUPS |
|---|---|---|---|
| 50² | 50k | 6 MB | 2.0 |
| 100² | 200k | 26 MB | 1.3 |
| 200² | 800k | 102 MB | 0.63 |

Cells reach their particles through linked lists that jump randomly through the particle array (128 bytes per particle), so once the array outgrows the L3 cache most accesses go to main memory. **Suggestion:** every N steps (opt-in, e.g. `"sortEvery": 50`), sort the particle array by cell index and rebuild the lists; `bin()` already rebuilds them every step. Separately, the `localPROP()` call after re-binning at the end of each step (`mpc.c:4953`) is repeated before any collision in the next step (`mpc.c:4641`); guarding it with "an output is due this step" saves ~15–18% of runtime (gprof).

## 5. PUPS report (suggestion)

`int mpcUpdates = (simSteps + warmupSteps) * GPOP` (`mpcd.c:342`) overflows a 32-bit int for large runs (e.g. 25,000 steps × 800k particles), and `warmupSteps` has already been reset to 0 at that point (`mpcd.c:276`), so warm-up steps are not counted although their CPU time is. Use `double` and save the warm-up count before resetting it.

## 6. Legacy `.inp` input: `tolD` uninitialised (fixed in patch 0002)

`readin()` (`read.c:53`) never sets `tolD` (and `inputList` is a local in `main()`), so the defect tolerance is undefined with legacy input. Patch 0002 sets the same defaults as `readJson()`.

## 7. Docs: `sigWidth` / `sigPos` (suggestion)

`InputGuide.md:177–178` say these apply "specifically for collision operator #20", but `lc.c` applies the sigmoid fall-off to both `DIPOLE_DIR_SIG` (20) and `DIPOLE_DIR_SIG_SUM` (21). Worth adding that with the defaults the fall-off factor at the mean density is ½(1 + tanh 1) ≈ 0.88, i.e. the effective activity is 12% below the nominal value everywhere.

## 8. Observation: activity scale shift (not investigated)

With the parameters of Kozhukhov, Loewe & Shendruk, *Commun. Phys.* (2024) (200², ρ̄ = 20, `collOp` 21), turbulence set in near α ≈ 0.02 there; with the current code even 0.05 is only mildly turbulent. We haven't tracked down which change caused this. A bisection over the git history with a fixed input (e.g. `activeNem2D`, comparing defect counts) would find it.

---

## Separately: custom wall anchoring (`customAnch`, `eMode`)

Our cavity work uses a custom wall anchoring in which the anchoring axis winds `customAnch` = k times around the wall (it originates in my old `BenjaminCustomBC` branch, 2022). Patch `customanch/0001` adds it cleanly on top of the three fixes, with `eMode` 0/1/2 (normal-angle, ellipse-parameter and arclength-uniform winding), documentation and validation. It is independent of the fixes above and can be reviewed separately, or not at all.

## Proposed next step

We're happy to follow whatever workflow the group uses now: open issues, or one pull request per patch (seed, tracker, Makefile) from a branch for review. The patches apply cleanly to current `master`.
