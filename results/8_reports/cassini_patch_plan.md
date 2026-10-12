# Cassini-oval cavity: patch plan

Written 2026-09-30. For H8 (figure-8 cavity, `RESEARCH.md`): a cavity with one wall and fixed winding k, pinched so that each lobe can host its own vortex.

**Status 2026-10-01: steps 1–3 implemented on branch `cassini` (commits 362ca8b, e63155b), validation steps 1–3 passed (`RESEARCH.md` results log), step 4 (active pilot) running in `runs/cassini/P1`. Differences from the plan: the JSON key is `"shape": 1` (int) rather than a string; `P` stays [2,2,2,2] in the input.**

## Where to do it: a separate copy of the code

This is a major change, so it must **not** be made in the production tree `mpcd/`, which the office runs use and the server build is synced from. Work in a separate copy under version control:

1. Use the git repository `~/Physics/Overcharged/niamh-overcharged/` (private GitHub `git@github.com:baloewey/NIAMH-MPCD-overcharged.git`, created 2026-09-30). Its `main` is the production code: upstream 2a38960 (PR #160) + seed fix + union-find tracker + customAnch/eMode 0/1/2, source identical to `mpcd/`.
2. Create the branch from `main`: `git switch -c cassini main`.
3. Build in the repository (not in `mpcd/`) and check a short passive and a short active cavity give the same result as the production binary with the same seed.
4. One commit per step below; push the branch to `origin` (never to `upstream`, whose push URL is disabled).
5. Only after validation (step 4 below) does it go into `mpcd/`, and never while a sweep is using that binary.


## Why a patch is needed

NIAMH-MPCD walls are not a free implicit function. `calcW()`/`surf_func()`/`normal()` only handle W = Σ|A_i(x_i − Q_i)|^{P_i} − R^{P_4}, plus a wavy term B₀cos(B₁θ) and a non-4-fold branch (`ROTSYMM`). The Cassini oval ((x−a)² + y²)((x+a)² + y²) = b⁴ has an x²y² cross term, so it fits none of them.

No-patch alternative for a quick pilot: the existing wavy wall with P = 2, B = [B₀, 2, 0] gives r² = R² − B₀cos2θ (circle at B₀ = 0, two touching circles at |B₀| = R², concave neck for |B₀| > R²/3). The Cassini family is preferred for the study: concave neck for b/a < √2, lemniscate at b = a, two separate ovals for b < a.

## 1. Input and data

- `bc` struct (`SRDclss.h`): `int SHAPE` (0 = current power-sum form, 1 = Cassini), `double CASSA` (focal half-distance a). b reuses `R`; `aInv` still scales x, y (1 = plain Cassini).
- JSON (`read.c`): `"shape": "cassini"`, `"cassiniA": a`; defaults keep every existing input unchanged.
- Dispatch on `SHAPE` first everywhere, so a Cassini wall never enters the P = 2 quadratic shortcut in `crosstime()`.
- Checkpoint write/read (`pout.c` ~l. 2525, `read.c` ~l. 789) and synopsis: include the new fields. Otherwise a restart silently turns the wall into a circle.
- Guards: exit with an error for a mobile Cassini wall, or with MD/swimmers (`calcW_BC`, `calcW_MD`, `calcW_swimmer` have their own copies of the formula).
- Document in `docs/InputGuide.md`.

## 2. Surface function and normal

Use W = d₁d₂ − b² (d₁, d₂ distances to the foci) instead of the quartic d₁²d₂² − b⁴: same surface, but W scales like the circle's r² − R², so the secant solver's fixed tolerances behave as for circles, and a = 0 gives exactly the circle's W.

```c
/// Cassini oval d1*d2 = b^2 (foci at ±a along x, in the wall frame); b = WALL.R, a = WALL.CASSA.
double cassiniW( bc WALL, double POS[] ) {
	double X = WALL.A[0]*(POS[0]-WALL.Q[0]), Y = WALL.A[1]*(POS[1]-WALL.Q[1]);
	double a = WALL.CASSA;
	double d1 = sqrt( (X-a)*(X-a) + Y*Y ), d2 = sqrt( (X+a)*(X+a) + Y*Y );
	double W = d1*d2 - WALL.R*WALL.R;
	if( WALL.INV ) W *= -1.0;
	return W;
}

/// Gradient of d1*d2 (not inverted, same convention as the ellipse branch of normal()).
double *cassiniNormal( double *n, bc WALL, double *point ) {
	double X = WALL.A[0]*(point[0]-WALL.Q[0]), Y = WALL.A[1]*(point[1]-WALL.Q[1]);
	double a = WALL.CASSA;
	double d1 = sqrt( (X-a)*(X-a) + Y*Y ), d2 = sqrt( (X+a)*(X+a) + Y*Y );
	n[0] = WALL.A[0] * ( d2/d1*(X-a) + d1/d2*(X+a) );
	n[1] = WALL.A[1] * ( d2/d1 + d1/d2 ) * Y;
	n[2] = 0.0;
	return n;
}
```

Singular only at the foci (inside the cavity, away from the wall for b > a).

Hook at the top of `calcW()` (`bc.c`), `surf_func()` (`mtools.c`), `normal()` (`bc.c`), and send Cassini walls to `secant_time()` in `crosstime()` / `crosstimeReverse()`:

```c
if( WALL.SHAPE==SHAPE_CASSINI ) return cassiniW( WALL,P.Q );
```

## 3. Anchoring: eMode 3, arclength by polar angle

eMode 2 recovers the wall position from the normal angle, which runs backwards on a concave neck. A Cassini oval with b > a is star-shaped about its centre, so the polar angle θ of the particle position labels the wall uniquely.

```c
else if( feq(WALL->EMODE,3.0) && WALL->SHAPE==SHAPE_CASSINI ) {
	double th = fmod( atan2( pp->Q[1]-WALL->Q[1], pp->Q[0]-WALL->Q[0] ) + 2.0*pi, 2.0*pi );
	double tau = cassiniArclengthFraction( WALL->CASSA, WALL->R, th );   // s(th)/P, cached table
	phi_bc = WALL->CA*2.0*pi*tau;
	nO[0] = cos(phi_bc + WALL->O[2]);
	nO[1] = sin(phi_bc + WALL->O[2]);
}
```

- `cassiniArclengthFraction()`: same structure as `anchArclengthFraction()` (cached table, linear interpolation); integrand √(r² + r′²) with r² = a²cos2θ + √(b⁴ − a⁴sin²2θ), r′ by numerical differentiation.
- `oriBC()` runs after `posBC()`, so `pp->Q` has already moved by the bounce; the θ error is tiny over one step. Pass the collision point explicitly if exactness matters.
- Assumes `O[2]` = 0 (the open rotation item in `TODO.md` applies).

## 4. Validation, in order

1. **a = 0 vs a circle of radius b**, same seed: W is identical, so statistics must match (passive custom anchoring: alignment, charge audit). Small differences expected from the secant vs quadratic crossing time.
2. **Accessible volume** vs the numerically computed Cassini area; no particle outside the wall after streaming (the "Darn" debug check in `MPC_BCcollision`).
3. **Passive** b/a = 1.6 (convex) and 1.2 (concave neck), k = 2…5: wall director follows the imposed arclength winding; charge audit net = k in every frame, including at the neck.
4. **Active pilot** with R_lobe ≳ 1.4 ℓ_a per lobe: flow field, one vortex or two, and the rotation sense of each.

## 5. Around the solver

- `gen_cavity.py --shape cassini --a … --b …`: box from the half-length √(a² + b²) and waist √(b² − a²); population from the numerical area.
- `circulation.py`: circulation per lobe (vortex centres from the vorticity).
- `edge_sites.py` and the charge audit (`defect_audit.py`) already work for any shape.

## Effort

About half a day for the C and the input/checkpoint plumbing, plus a day of validation runs.
