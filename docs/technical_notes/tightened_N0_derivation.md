# Margin-driven uniform large-\(N\) certificate

**Result.** For \(b=10,\phi=\pi\), every odd \(N\ge 10{,}001\) satisfies C1–C8 of the small-spacing theorem. This replaces the former sufficient floor \(35{,}377\) by a factor \(3.54\). It is a proof threshold, not a claimed physical or minimal length. Five smaller odd lengths remain covered by their separate direct certificates.

## Kernel remainder without a prescribed \(10^{-30}\) target

Write \(t=N^{-2}\), \(d_t(q)=d(q)H(q^2t/4)\), and \(H(w)=\sqrt w/\sin\sqrt w=\sum_{m\ge0}h_mw^m\). The checker uses exact rational \(h_0,\ldots,h_8\). On \(|w|=4\), the sine series gives
\[
\left|\frac{\sin\sqrt w}{\sqrt w}-1\right|
\le \frac{\sinh 2}{2}-1<\frac56,
\qquad |H(w)|<6.
\]
Cauchy's bound is therefore \(|h_m|<6/4^m\). The integral formula \(d(q)=\int_{-1/2}^{1/2}e^{iqs}\,ds\) gives \(|d^{(j)}(q)|\le2^{-j}\) for real \(q\). Applying Leibniz's rule to each omitted term and using \((2m)_j\le(2m)^j\), for every \(0\le k\le46\) and \(|q|\le110\),
\[
\left|\partial_q^k\left[d_t(q)-d(q)\sum_{m=0}^8h_m(q^2t/4)^m\right]\right|
\le 6\sum_{m=9}^\infty R^m
\max\left\{1,\left(\frac12+\frac{2m}{110}\right)^{46}\right\},
\quad R=\frac{110^2}{16N^2}.
\]
The machine checker sums the first 92 terms as exact rationals and bounds the remaining positive tail geometrically. At \(N=10{,}001\) this upper bound is below \(4.846\times10^{-46}\), so the directed evaluator uses \(\pm10^{-45}\). The bound decreases with \(N\). No arbitrary \(10^{-30}\) target enters the threshold choice.

The near-central chart evaluates the analytic deflations \(A_t(v)/v^4\) and \(B_t(v)/v^2\) by Taylor polynomials through degree 40. On a complex \(v\)-circle of radius 10, the finite-kernel representation at \(N\ge10{,}001\) and \(|\lambda|<0.067\) bounds the deflated functions by \(10^6\). Cauchy's coefficient bound, including two \(v\)-derivatives on \(|v|\le1\), puts the omitted tail below \(10^{-25}\); the code widens by that amount. This enclosure is far below the tightest final predicate margin and does not set the floor.

## Direct margin inequalities

The proof does not use an unpropagated comparison between the kernel remainder and a loss margin. It evaluates each *complete interval predicate* for the whole parameter interval \(0\le t\le10{,}001^{-2}\). The joint Krawczyk map is strictly included, its contraction is below \(0.142794\), the tangent Gram factors and A/B curvatures are positive, the satellite coefficients have fixed nonzero signs, and the equal-cost slope exceeds \(0.1078127\). The two root boxes remain inside their designated compact-cover cells.

The 27,708-cell full-line continuum partition is replayed as a uniform finite-\(N\) cover. Four frozen predicates fail at this floor. Two of these cost cells pass after one exact bisection each; the other two pass by a stronger gradient or monotonicity predicate without splitting. Thus four archived parents are replaced by six positive-margin terminal leaves. Every other original cell passes. The smallest resulting compact-cover margins are:

| Predicate | Uniform lower margin |
|---|---:|
| cost exclusion | \(4.65369\times10^{-9}\) |
| gradient exclusion | \(5.43958\times10^{-9}\) |
| monotone no-root cell | \(4.12411\times10^{-10}\) |
| A/B root-cell curvature | \(2.50644\times10^{-5}\) / \(1.01184\times10^{-4}\) |
| coalescent cost gap | \(8.51154\times10^{-4}\) |
| finite-period tail cost gap | \(3.13564\times10^{-3}\) |

The normalized finite period expands with \(N\); the separate summation-by-parts tail bound is replayed for all allowed \(t\). Positive margins for the other conditions, their continuum counterparts, and the resulting *joint sufficient* \(N\) floor are in [largeN_condition_thresholds.csv](largeN_condition_thresholds.csv). The column \(N=10{,}001\) is a certified common floor, not a condition-specific minimum.

For orientation, predicate-wise tested sufficient floors are also tabulated. Local Gram, curvature, satellite, slope, coalescent, and finite-period-tail signs already pass on the wider trial interval ending at \(N=1{,}501\). The joint root Krawczyk inclusion fails there but passes at 5,001. The compact regular-competitor cover still has 19 failed cells at 5,001 and passes after targeted refinement at 10,001; it sets the tested joint floor. These are tested sufficient floors, not smallest possible odd lengths. A lower failed interval predicate does not demonstrate failure of the theorem there.

## Active bottleneck and limits

The previous floor was triggered by requiring a specific kernel error of \(10^{-30}\). At the new floor, kernel truncation is much smaller than the proof margins; the active difficulty is preservation of a few cost/derivative predicates in the frozen competitor partition. At \(N_0=5{,}001\), 19 frozen cells fail and the original A-root cell does not contain the parameterized root enclosure. At \(N_0=1{,}501\), the initial local box and at least 30 cells fail. These are failures of the reused enclosure, not evidence that the branch exchange disappears. Direct finite-\(N\) interval blocks below \(10{,}001\) and a recentered root partition could lower the floor further.

Increasing the series order beyond \(m=8\) would shrink a remainder that is already inactive at \(N_0=10{,}001\). A different compact partition or parameterized direct finite sum is the relevant route to a much smaller threshold. No minimality claim is made. [tightened_N0_certificate.json](tightened_N0_certificate.json) records the composite directed proof and input hashes.
