# N=31 independent finite-spacing global certificate

## Exact instance and conclusion

For `t=-15,...,15`, `s=t/31`, fit two undamped complex exponentials with free complex amplitudes and free real frequencies to

`x_t = exp(-2 i s) + exp(2 i s) - ε exp(10 i s)`.

Thus `z=NΔ=2`, the strong frequencies are `ω=±2/31`, the weak frequency is `ω=10/31` radians/sample, the strong amplitudes are one, and the weak tone has phase π and positive magnitude ε. The N=21 comparator uses the same normalized frequencies `u=Nω=(-2,2,10)`, so the physical frequencies scale as `1/N`. The fitted objective is the exact sample-domain variable-projection loss over **all unordered two-frequency pairs modulo 2π**, including the confluent limit.

Within the bracket below, branches A and B are the only globally competitive fits. They exchange global optimality once, transversely, without fitted-frequency coalescence. This is an independent N=31 point certificate; it does not establish a connected `z` interval or a theorem for arbitrary N.

## Crossing and local conditions

| Quantity | N=31 result |
|---|---:|
| Multiprecision stationary equal-cost root `εc` | `0.2491430898525017772068990221341981098904119007017468553511825943724723` |
| Validated `ε` bracket | `[0.2491430898524917772068990221341981098904119007017468553511825943724723, 0.2491430898525117772068990221341981098904119007017468553511825943724723]` |
| `J_A=J_B` at numerical root | `1.21785722117124279980487926228355670033549596953796050332254` |
| Branch A fitted `u=Nω` | `(-4.7323081555003052296, 0.6135741125336810595)` |
| Branch B fitted `u=Nω` | `(-0.0573137795963794694, 12.440325828105324947)` |
| A / B fitted separation in `u` | `5.3458822680` / `12.4976396077` |
| Interval slope `d(J_A-J_B)/dε` | `[12.2910730810731, 12.2910731034833]` |
| Ordering | A wins below the bracket root; B wins above it. |

The lower endpoint has `J_A-J_B∈[-1.2291085264,-1.2291060921]×10⁻¹³`; the upper endpoint has `J_A-J_B∈[1.2291060921,1.2291085264]×10⁻¹³`. Krawczyk inclusions and positive Hessians hold at both endpoints and **for the full bracket**. The slope interval is strictly positive, giving a unique A/B crossing within the bracket. See [N31_crossing_certificate.csv](N31_crossing_certificate.csv) and [N31_global_certificate.json](N31_global_certificate.json).

Each radius-`10⁻¹⁰` root box is strictly included by its Krawczyk image. The larger boxes used by global exclusion, radius `0.0003` in each A frequency and `0.01` in each B frequency, also have strict Krawczyk inclusion, contraction below one, and positive-definite interval Hessians throughout the epsilon bracket. Their Gram determinant lower bounds are `933.45127` (A) and `960.92616` (B). All fitted amplitude intervals are finite; the A amplitudes lie within `[0.339865,0.340577]` and `[1.656277,1.656813]`, and B within `[1.724691,1.733339]` and `[-0.247117,-0.228489]`. These wide bounds refer to the larger root boxes, not numerical root uncertainty.

## Global exclusion, including coalescence

Parameterize every unordered pair by a center frequency `m∈[-π,π]` and shortest half separation `h∈[0,π/2]`. The full-domain float search supplies only a partition. Every proposed exclusion is checked again with directed interval arithmetic over the **entire epsilon bracket**. Each accepted cell is independently checked, after `10⁻¹²` radian side inflation, to lie inside one of the two radius-`2.00000001` boxes in `u` coordinates. The inflated root domain includes the tiny float rounding gap at the exact `π` boundaries.

The outer bound uses the continuous orthogonal basis `exp(i m t) cos(h t)` and `exp(i m t) t sinc(h t)`. At `h=0` the second vector tends to `t exp(i m t)`; thus the coalescent limit is included without dividing by a nearly singular two-tone Gram matrix. Directed Taylor bounds on basis variation give a projector-distance bound and then a residual lower bound for each cell.

| Full-domain audit | Count / bound |
|---|---:|
| Terminal outer cells | 27,584 |
| Directed-interval excluded cells | 27,481 |
| Cells independently contained in A/B neighborhoods | 103 (69 A, 34 B) |
| Failed / unresolved outer cells | 0 / 0 |
| Minimum excluded objective margin over incumbent | `>1.08091×10⁻⁴` |
| Cells intersecting `0≤h≤0.03` | 281, all excluded |
| Minimum directed objective bound in that strip | `>1.2218118780` (incumbent near `1.2178572212`) |

Inside the accepted A and B neighborhoods, interval Taylor objective bounds, interval-gradient exclusion, and Krawczyk exclusion visit 1,371 and 277 cells respectively; neither search leaves an unresolved cell. Only cells fully contained in the larger, independently validated root boxes survive. The complete outer log is [N31_interval_cells.csv](N31_interval_cells.csv); the inner terminal-cell log is [N31_inner_interval_cells.csv](N31_inner_interval_cells.csv). The executable implementation is [stage7_n31_certificate.py](stage7_n31_certificate.py), with the inner-log replay in [stage7_n31_inner_log.py](stage7_n31_inner_log.py).

## Analytical comparison

The N=31 tangent calculation gives `λ₃₁=0.0662657040848183134471...` and `c₃₁=-0.00103037885628279811372...` for the **local** small-`z` branch law. At `z=2`:

| Approximation | Value | Absolute error from numerical crossing | Relative error |
|---|---:|---:|---:|
| `λ₃₁ z²` | `0.26506281633927325379` | `0.01591972648677147658` | `6.38979%` |
| `λ₃₁ z²+c₃₁ z⁴` | `0.24857675463874848397` | `0.00056633521375329324` | `0.227313%` |
| Interval-certified crossing | the bracket above | width `2×10⁻¹⁴` | — |

This is a pointwise accuracy comparison at `z=2`, not a validated asymptotic remainder bound through `z=2`.

## Arithmetic and scope

The root was initialized from Stage 5 and refined with 100-digit `mpmath`; local intervals use 90 decimal digits, outer intervals 70, and inner intervals 45. Float64 proposes cells and a Hessian inverse; all deciding inclusion, positivity, containment, and exclusion inequalities are interval checked. The directed-rounding claim relies on `mpmath.iv`; no independent interval-library audit was performed. The result is for this exact symmetric geometry and narrow epsilon bracket.

**Verdict: SECOND-N CERTIFIED.**
