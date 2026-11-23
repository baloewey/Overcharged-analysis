# Curvature (1/R) corrections to the tilt law — first-order theory (2026-11-23, user)

Circle of radius R, no-slip wall, single vortex with azimuthal velocity v(r) (signed with the circulation), depth below the wall y = R − r.
Director angle measured from the local tangent (flow direction), α, positive with the circulation. Flow-aligning parameter λ (Leslie angle
θ_L: cos 2θ_L = 1/λ). Small parameters: ε = y/R (curvature) and ℓ_E/R; the flat-wall law is the ε → 0 limit.

## Link A — kinematics of a curved no-slip wall (exact, then first order)

For any axisymmetric azimuthal flow, vorticity and strain are

    ω = v' + v/r,     e_rθ = ½ r d(v/r)/dr = ½ (v' − v/r),     e_rr = e_θθ = 0,

so the extension axis is always at 45° from the local tangent (curvature does not rotate it). Only the ratio of rotation to strain changes:

    Λ(y) ≡ |ω| / (2|e_rθ|) = (|v'| − v/r) / (|v'| + v/r)        (v' < 0 near the wall: speed grows inward)

Near the wall, v = s y + ½ s₂ y² + … (s = wall shear rate). Then v/r = s y/R + O(y²/R²) and |v'| = s + O(y), so

    Λ(y) = 1 − 2y/R + O(y²/R², y² s₂/(s R))                                  [A1]

At first order the correction is purely geometric and independent of the flow profile: a curved no-slip wall is more extensional than simple
shear by 2y/R. Deeper in, the profile matters (solid-body rotation near the vortex centre, Λ → ∞: tumbling); e.g. a vortex driven by a force
∝ r gives v ∝ r(1 − r²/R²) and Λ = |2 − 1/x²| (x = r/R), aligned (Λ < λ = 2) only for x > ½.

## Link B — Jeffery alignment with elasticity

2D Jeffery (Larson) rotation of the director: θ̇ = ω/2 + λ |e| sin 2(χ − θ), χ = extension axis. Steady, flow-dominated (no elasticity):

    cos 2α_J(y) = Λ(y)/λ   →   α_J(y) = θ_L + (y/R) / √(λ² − 1) + O(ε²)          [B1]

(λ = 2: θ_L = 30°, +33° per unit y/R.) With one-constant elasticity, the director angle relative to the tangent, θ = φ + π/2 + α(r), has
∇²θ = α'' + α'/r (φ is harmonic), so the torque balance reads

    K (α_yy − α_y / R) = γ₁ λ |e| sin 2(α − α_J(y)) ≈ (K/ℓ_E²) (α − α_J(y)),     ℓ_E² = K / (γ₁ s √(λ² − 1))      [B2]

α relaxes from its wall value to α_J(y) over the alignment length ℓ_E; curvature enters (i) through α_J(y) [B1] and (ii) through the −α_y/R
term, which changes the relaxation length to ℓ_E (1 + ℓ_E/2R + …) — a second-order effect on the tilt when α_y is small far from the wall.

## Link C — the zero mode: layer tilt = 2 × rotation of the field above it

Flat wall: sliding the screening +½ by D rotates the field above it by κD and the comet by 2κD (θ = φ/2 + θ₀ around a +½: rotating the
surroundings by β turns the tail by 2β). Circle: the one-constant circle maps conformally onto the flat wall, the lock τ_j = (2k − 2)φ_j + c
is exact, and the 2k − 2 ring slides rigidly: **no curvature correction in link C**. Hence

    τ = 2 α(y*) = arccos(1/λ) + 2 (y*/R) / √(λ² − 1) + O(ε², ℓ_E/R)              [C1]

where y* is the depth at which the field "above the layer" is read (between the core depth d and d + ℓ_E). λ = 2: +66° per unit y*/R;
for y* ≈ 2 cells, +6° at R = ℓ_a and +3° at R = 2ℓ_a.

## Other size effects that are NOT curvature (must be separated)

- Wall shear s grows with R/ℓ_a (stronger vortex): enters only through ℓ_E (B2), i.e. through y* and the depth correction.
- Wall-defect spacing a = πR/(k − 1): crowding lowers the tilt (measured; mechanism open). Test at fixed ρ = 2k/P.
- Active self-flow of each wall +½ (extensional, |E|/W 0.87–0.91 at the core in active circles vs 1.06–1.12 passive; defect_flow).

## Predictions to test (in this order)

| Link | Prediction | Result (979 circles, R 0.8–2.2 ℓ_a, curvature_profiles.txt; 23 Nov) |
|---|---|---|
| A | 1 − Λ = 2ℓ_v/R, ℓ_v = v/\|v'\| (= y only for a linear profile) | **Confirmed, 1/R scaling clean**: sparse, y 2.25: 1 − Λ = 0.105 / 0.069 / 0.051 vs 2ℓ_v/R = 0.098 / 0.069 / 0.052 (R 22 / 32 / 42). ℓ_v ≈ y/2 near the wall (velocity rises faster than linearly across the defect layer) |
| B | α → α_J(y) = ½ arccos(Λ/λ) beyond ℓ_E | **Fails, with a 1/R correction of the opposite sign, ~10× the Jeffery one**: α_J changes by 0.5° (30.9 → 30.4°) near the wall, while the measured α − α_J at y 1.25 is −8.8 / −4.5 / −3.8° (sparse) and −10.8 / −6.2 / −2.9° (crowded) for R 22 / 32 / 42; deeper, α falls towards 0 and Λ > λ (tumbling) beyond 6–8 cells. Hypothesis: elastic coupling of the aligned annulus to the tumbling interior, whose distance from the wall scales with R → pull ∝ 1/R. Next: solve [B2] in the annulus with the measured α_J(y) and a matching condition at the tumbling radius |
| C | τ = 2α(y*) | **Fails for the ring-averaged α**: sparse tilts 64–65° (flat in R) while 2α ≤ 53° at every depth. The tilt is not set by the ring average; October's defect-centred analysis (defect_flow_main10) points to the local flow at each defect → reformulate link C locally |

The curvature correction of the kinematics is real and small near the wall; the director's size dependence is dominated by the elastic connection to the
tumbling interior; and the defect tilt follows neither (it is set locally at each defect).
