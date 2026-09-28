# Finite-amplitude certificate audit and attempted extension

## Historical status

The N=21 cases at εc−0.035, εc, εc+0.035 are **full interval-certified point cases**, not merely high-precision illustrations. The outer logs contain 12,288/14,287/10,961 excluded cells and 76/90/51 accepted cells. The respective inner searches complete with zero unresolved cells. The original clean replay regenerated those records. Figure 1 may call these three amplitudes separately certified; it may not claim the whole interval between them.

## Finite-width attempt

The requested [εc−0.03,εc+0.03] was split into 12 slabs of width 0.005. Each slab used a multiprecision root proposal at its midpoint, parameter-interval gradients/Hessians and Krawczyk boxes of radii 0.005,0.01,0.025,0.05,0.1 in normalized frequency. A uniform incumbent came from a fixed feasible frequency pair evaluated with interval ε. Forty weakest historical outer cells were tested in each slab using interval ε.

None of these slab attempts closed: no tested root-box radius passed both inclusion and positive-Hessian checks, and none of the 40 weakest reused spatial cells passed in any slab. These failures are recorded rather than silently discarded. The original natural interval expression overestimates frequency/amplitude dependency; shrinking a scalar bracket alone is not a proof of a finite-width regime. No complete new outer/inner partition was generated after this failed prerequisite diagnostic. The finite-width result is **CERTIFICATION-LIMITED**, not certified and not mathematically refuted. A predictor-centered parameterized Krawczyk/Taylor model and a new spatial partition would be needed to continue this attempt.

## Location of the minimum historical margin

The attaining physical-frequency cell is

    m∈[0.0337475773334841,0.03681553890925539]
    h∈[0.07363107781851078,0.07669903939428206].

Its center is u≈(−0.83755351019,2.31937895128), and its coordinate enclosure is

    u1∈[−0.90198070328,−0.77312631709],
    u2∈[2.25495175819,2.38380614437].

The historical lower bound is 0.8118137698122954, only about 5.8978946e−6 above the incumbent. The actual center cost is 1.27373900320810, about 0.46193 higher. Its center is at L∞ distance 3.93489 from A and 10.14404 from B, or 1.93489 and 8.14404 outside the corresponding radius-two boxes. A multiprecision unconstrained stationary solve from this center converges to A outside the cell. More decisively, a 31-node directed interval-gradient subdivision excludes stationary points throughout the cell. It is **enclosure slack, not a third basin**.

All attempted slab boxes, failed cells, bounds, midpoint roots, and the gradient-exclusion result are in `finite_epsilon_certificate.json`. The pre-existing narrow crossing bracket remains valid; no finite-width claim has been substituted for it.
