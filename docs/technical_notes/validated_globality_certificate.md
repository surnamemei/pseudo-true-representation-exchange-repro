# Validated globality certificate — 21-sample baseline

**Model.** For `t = −10,…,10`, let `s=t/21` and `x_t=2 cos(2s)−ε exp(i10s)`. We minimize `J(u1,u2)=min_{c1,c2∈C} ||x−c1 exp(iu1s)−c2 exp(iu2s)||²`, where `u_j=21ν_j`, over all unordered frequencies modulo `2π`. The weak-tone phase is π. This file indexes the reproducible [machine-readable certificate](validated_globality_certificate.json); all outer terminal cells and their interval lower bounds are in [interval_boxes.csv](interval_boxes.csv), and inner terminal cells are in [stage2_inner_cells.csv](stage2_inner_cells.csv).

## Result

| ε case | ε | global branch | A/B interval objective ordering | outer excluded cells | smallest excluded J margin | coalescent-strip margin |
|---|---:|---|---|---:|---:|---:|
| below | 0.2131906301722773982 | A | `J_A < J_B`, gap ≈0.272124 | 12,288 | 2.7084×10⁻⁵ | 0.0178572 |
| crossing reference | 0.2481906301722773982 | A/B exchange | both ≈0.81180787191753444 | 14,287 | 5.8979×10⁻⁶ | 0.0053768 |
| above | 0.2831906301722773982 | B | `J_B < J_A`, gap ≈0.306964 | 10,961 | 8.6576×10⁻⁶ | 0.0028048 |

The exact crossing lies in

`[0.248190630172276398203609851355…, 0.248190630172278398203609851355…]`.

The interval slope `d(J_A−J_B)/dε` over the bracket is `[8.25905103823729955…, 8.25905105329243140…]`, so the crossing is transverse. At its endpoints, interval `J_A−J_B` is strictly negative and positive respectively. Both local Krawczyk inclusions, both full inner searches, and all **14,287** outer interval exclusions were repeated across the **entire** ε bracket; see [stage2_crossing_global.json](stage2_crossing_global.json). Thus the globally selected branch changes inside that bracket.

At the reference crossing, the fitted frequency pairs are `A=(-4.7724473539186143, 0.6051703202768620)` and `B=(-0.0576535619266798, 12.463419356801660)` in `u=Nν` units. Their internal separations are ≈5.37762 and ≈12.52107, so neither solution is coalescent. The absolute pair displacement is much larger than the certified local boxes.

## How the interval exclusion works

1. **Separated boxes.** All inner frequency boxes have radius `2.00000001` in each `u_j`; their Gram determinants stay positive. `three_to_two_tone_stage2_interval.py` evaluates the exact finite exponential sums and analytically differentiated reduced 2×2 Gram formula with interval second-order automatic differentiation. A floating Hessian inverse is used only as a fixed Krawczyk preconditioner; the inclusion and exclusion inequalities use intervals. Taylor objective lower bounds, mean-value gradient exclusion, and Krawczyk exclusion partition each inner box. The six searches visit 105–1,571 cells each, leave **zero** unresolved cells, and admit cells only inside larger certified root boxes of radius `0.0003` (A) or `0.01` (B). Those larger boxes have Krawczyk inclusion, contraction below one, and positive-definite Hessians at all three ε values. The final root coordinates are enclosed in radius-`10⁻¹⁰` boxes.
2. **Full torus.** Write `m=(ν1+ν2)/2` and `h=(ν2−ν1)/2`, choosing the shorter circular separation. `m∈[−π,π]` and `h∈[0,π/2]` cover every unordered pair. The outer bound uses the continuous orthonormal basis `q1∝e^{imt} cos(ht)` and `q2∝e^{imt} t sinc(ht)`. Its second column tends to `t e^{imt}` as `h→0`, so no ill-conditioned direct two-column inverse appears at coalescence. Taylor bounds for the basis change bound the projector distance from each cell center. The triangle inequality gives `J_cell ≥ max(0, sqrt(J_center)−||x||·||P−P_center||)²`. Float64 only proposes an adaptive partition; every excluded cell is re-evaluated with `mpmath.iv` at 65 decimal digits. Each accepted cell, inflated by `10⁻¹²` radians per `m,h` side, was independently checked to lie inside an inner box. The smallest acceptance containment margin is at least `0.004268` in `u` units. The terminal partition has no unresolved leaves.
3. **Coalescence.** The explicit strip `0≤h≤0.03` includes the confluent limit and corresponds to `0≤u2−u1≤1.26`. It intersects 199, 209, and 203 excluded cells below, at, and above the crossing; **none** is accepted. The table's positive margins bound the entire strip. Therefore a near-coalescent or confluent approximation cannot win in the decisive cases.
4. **Local minima.** At all six larger root boxes, interval Hessian eigenvalue lower bounds are positive. The smallest of these conservative bounds is `0.030761` (below/A). Local existence and uniqueness follow from strict Krawczyk inclusion and contraction; strict local minimality follows from Hessian positive definiteness. The global partition reduces every competitive point to these roots.

## Arithmetic and reproducibility

Reference roots and objectives use 110-digit `mpmath`; local and crossing intervals use 90 digits; outer cells use 65 digits; inner cells use 45 digits. Library: `mpmath 1.3.0` (`mpmath.iv` with directed endpoint rounding). The [certificate JSON](validated_globality_certificate.json) records versions, every root box, Hessian bounds, cell counts, margins, and the crossing bracket. `interval_boxes.csv` is the full outer terminal-cell log. Scripts `three_to_two_tone_stage2_{reference,interval,krawczyk,outer_global,inner_global,crossing,crossing_global,cell_audit,certificate_build}.py` reproduce the computations.

The independent float64 continuation gives `εc` within `1.39×10⁻¹⁶` of the 110-digit reference, but its fitted frequencies differ by as much as `4.86×10⁻⁸` and objectives by about `10⁻¹⁵`; float64 alone is insufficient for the narrow crossing sign test.

**Scope of “validated.”** The inequalities close using the installed interval implementation and the stated exact finite-sum formulas. The interval library itself was not independently audited or cross-checked against a second outward-rounded package. The result is a reproducible computer-assisted numerical certificate for this baseline geometry and three ε cases (plus the narrow crossing bracket), not a symbolic theorem for all nearby geometries.

