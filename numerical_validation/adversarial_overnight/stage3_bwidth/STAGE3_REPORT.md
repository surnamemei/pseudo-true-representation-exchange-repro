# Stage 3 — parameter width in b

## 3A. Numerical scan (continuum tangent problem; diagnostic only, never 'certified')

Independent Gauss–Legendre (256-node) realified least-squares implementation, cross-checked at b = 10 against the frozen closed form: max |ΔR| = 4.8e-17, |Δλ∞| = 1.8e-16 (pass).

| b | λ∞(b) | v_A | v_B | R_vv,A | R_vv,B | β_A | β_B | slope ∂λ(R_A−R_B) | coalescent gap | nearest regular competitor gap | tail margin | pass |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 9.5 | 0.06618436 | -4.4685 | 12.2206 | 1.237e-04 | 8.169e-04 | 0.1016 | -0.0593 | 0.1015 | 8.908e-04 | 2.640e-03 | 3.153e-03 | True |
| 9.6 | 0.06621017 | -4.3889 | 12.2608 | 1.257e-04 | 8.444e-04 | 0.1061 | -0.0603 | 0.1028 | 8.856e-04 | 2.738e-03 | 3.266e-03 | True |
| 9.7 | 0.06625613 | -4.3083 | 12.3012 | 1.277e-04 | 8.726e-04 | 0.1109 | -0.0613 | 0.1041 | 8.790e-04 | 2.839e-03 | 3.382e-03 | True |
| 9.8 | 0.06632145 | -4.2266 | 12.3418 | 1.297e-04 | 9.013e-04 | 0.1160 | -0.0623 | 0.1054 | 8.711e-04 | 2.942e-03 | 3.502e-03 | True |
| 9.9 | 0.06640533 | -4.1439 | 12.3827 | 1.317e-04 | 9.306e-04 | 0.1216 | -0.0634 | 0.1066 | 8.618e-04 | 3.049e-03 | 3.624e-03 | True |
| 10.0 | 0.06650697 | -4.0600 | 12.4239 | 1.337e-04 | 9.604e-04 | 0.1275 | -0.0644 | 0.1078 | 8.512e-04 | 3.159e-03 | 3.749e-03 | True |
| 10.1 | 0.06662558 | -3.9751 | 12.4653 | 1.357e-04 | 9.907e-04 | 0.1339 | -0.0655 | 0.1090 | 8.391e-04 | 3.272e-03 | 3.878e-03 | True |
| 10.2 | 0.06676037 | -3.8891 | 12.5070 | 1.377e-04 | 1.021e-03 | 0.1408 | -0.0665 | 0.1101 | 8.257e-04 | 3.388e-03 | 4.009e-03 | True |
| 10.3 | 0.06691054 | -3.8021 | 12.5489 | 1.397e-04 | 1.053e-03 | 0.1483 | -0.0676 | 0.1112 | 8.109e-04 | 3.507e-03 | 4.144e-03 | True |
| 10.4 | 0.06707528 | -3.7139 | 12.5910 | 1.417e-04 | 1.084e-03 | 0.1565 | -0.0686 | 0.1123 | 7.946e-04 | 3.629e-03 | 4.281e-03 | True |
| 10.5 | 0.06725379 | -3.6247 | 12.6334 | 1.437e-04 | 1.116e-03 | 0.1654 | -0.0697 | 0.1134 | 7.770e-04 | 3.754e-03 | 4.421e-03 | True |

**Largest symmetric interval supported by the frozen numerical criteria: b ∈ [9.5, 10.5] (radius 0.5, the largest radius tested).** Evidence level GLOBAL_NUMERICAL (tangent problem); diagnostic, not certified.

Secondary: finite N = 21, z = 2 record (full-torus crossing analysis at each b):

| b | verdict | ε_c | J* | third margin (rel.) | min Hessian eig |
|---|---|---|---|---|---|
| 9.5 | present | 0.24756731 | 0.719220 | 0.7410 (1.030) | 0.0478 |
| 9.6 | present | 0.24755953 | 0.737415 | 0.7679 (1.041) | 0.0491 |
| 9.7 | present | 0.24762039 | 0.755782 | 0.7955 (1.053) | 0.0505 |
| 9.8 | present | 0.24774747 | 0.774311 | 0.8239 (1.064) | 0.0518 |
| 9.9 | present | 0.24793836 | 0.792990 | 0.8530 (1.076) | 0.0532 |
| 10.0 | present | 0.24819063 | 0.811808 | 0.8830 (1.088) | 0.0546 |
| 10.1 | present | 0.24850187 | 0.830752 | 0.9137 (1.100) | 0.0561 |
| 10.2 | present | 0.24886966 | 0.849809 | 0.9452 (1.112) | 0.0576 |
| 10.3 | present | 0.24929155 | 0.868964 | 0.9775 (1.125) | 0.0591 |
| 10.4 | present | 0.24976508 | 0.888201 | 1.0106 (1.138) | 0.0606 |
| 10.5 | present | 0.25028774 | 0.907503 | 1.0445 (1.151) | 0.0622 |

b = 10 is not a tuned point: across b ∈ [9.5, 10.5] the continuum crossing coefficient changes by < 1.6% and every margin stays positive; the certified finite z = 2 exchange persists at every b with ε_c ∈ [0.2476, 0.2503] and third-competitor margins ≥ 1.03·J*.

## 3B. Attempted explicit interval certificate (continuum problem, uniform in b)

Verdict: **NO_EXPLICIT_B_INTERVAL_CERTIFIED** (runtime 108.0 min of the 150-min pre-registered cap; 0/24 attempted 0.01-wide blocks passed).

| candidate b interval | blocks | primary (mpmath.iv) status | Arb replay |
|---|---|---|---|
| [9.75, 10.25] | 50 | failed | — |
| [9.90, 10.10] | 20 | failed | — |
| [9.95, 10.05] | 10 | failed | — |
| [9.98, 10.02] | 4 | failed | — |
| [9.99, 10.01] | 2 | failed | — |

**Why no interval was certified (diagnosis).** Every one of the 24 completed blocks (b ∈ [9.88, 10.12]) passed all *local* obligations at all 100 sub-tiles: parametric Krawczyk inclusion (max contraction ≤ 0.67), positive reduced curvatures (min R_vv,A ≥ 4.7e-5), nonzero satellite coefficients (β_A ≥ 0.120, |β_B| ≥ 0.063), crossing slope ≥ 0.106, strict convexity of both ±0.1 windows, coalescent margin ≥ 8.3e-4 and tail margin ≥ 3.59e-3 (minima over all 24 blocks). Each block then exhausted the 2.5-million-visit far-field budget while sweeping v upward from −100: the sweep stopped at v ≈ −1.32. The flat region between well A and the central chart consumed the budget (≈ 5.7e5 leaves in v ∈ [−2, −1.32], 2.5e5 in [−3, −2]). The near-zero region and the whole B neighbourhood were never reached. The cause is dependency overestimation of naive mpmath.iv enclosures in b over flat regions: many cells had to be split down to single 1e-4 sub-tiles. It is not a certified counterexample: at every sampled b the numerical margins are positive (3A). The second-wave blocks (b beyond ±0.12) were abandoned at the 108-min primary deadline, and their workers were terminated after the verdict was written. A feasible future attempt would use Arb as the primary evaluator (measured ≈ 40× faster per cell here) and centred (mean-value) forms in b for the cost and gradient predicates.

Method: b tiled into 0.01 blocks × 100 sub-tiles of width 1e-4; per sub-tile a parametric Krawczyk inclusion for (v_A, v_B, λ) with a mean-value form in b; strict convexity of R on ±0.1 windows around the root boxes (unique stationary point per window); full-line exclusion on [−100, 100] by cost/gradient/concavity/monotone predicates with adaptive (v, b)-bisection; per sub-tile coalescent and |v| ≥ 100 tail margins; exact rational coverage of every sub-tile's partition; independent python-flint/Arb (256-bit) replay recomputing every enclosure, preconditioner, λ-enclosure and incumbent from the stored geometry only. The earlier Stage-23 single-box attempt failed because its boxes were ~50× too wide for the Hessian enclosure.