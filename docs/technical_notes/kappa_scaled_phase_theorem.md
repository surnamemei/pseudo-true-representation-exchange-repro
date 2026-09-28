# Phase-offset scaling \(\delta=\kappa z\): local theorem and numerical boundary

Let \(p=z^2\), \(\epsilon=\lambda p\), and \(\delta=\kappa z\) with fixed real \(\kappa\). On centered samples,

\[
2\cos(zs+\delta)=2\cos[z(s+\kappa)]
=2-p(s+\kappa)^2+\frac{p^2}{12}(s+\kappa)^4+O(p^3).
\]

The exact order-\(p\) target for the phase-opposed weak tone is
\(q_\kappa(s)=-(s+\kappa)^2-\lambda e^{ibs}\). The constant \(-\kappa^2\) is absorbed by the central complex amplitude. The remaining real-linear term \(-2\kappa s\) cannot be absorbed by a real frequency shift of one central exponential, but it is part of the tangent target. The regular one-central/one-satellite tangent loss is

\[
R_{N,\kappa}(v,\lambda)=
\min_{\alpha,\beta\in\mathbb C,\,\xi\in\mathbb R}
\sum_t|q_\kappa(s_t)-\alpha-2i\xi s_t-\beta e^{ivs_t}|^2.
\]

## Exact separation and symmetry

With \(D_N,S_2,G_0\) as in `fixed_delta_scaling.md`, real-even/imaginary-odd and real-odd/imaginary-even parts of the centered least-squares design decouple. The exact identity is

\[
R_{N,\kappa}(v,\lambda)=R_{N,0}(v,\lambda)
+4\kappa^2P_N(v),\qquad
P_N(v)=S_2-\frac{D_N'(v)^2}{N-D_N(v)^2/N}\ge0.
\]

The identity was checked independently against direct real least squares. The confluent chart contains complex \(1,s\) and real \(s^2\), so its tangent loss is independent of \(\kappa\). The transformation \(q_\kappa(s)\mapsto\overline{q_\kappa(-s)}=q_{-\kappa}(s)\) leaves the fitted frequency set and optimized cost unchanged. Therefore every isolated branch cost, its stationary frequency, and its local crossing \(\lambda_c(\kappa)\) are even functions of \(\kappa\). No linear \(\kappa\) term is allowed.

## Certified-local structural-stability corollary

At \(\kappa=0\), the existing finite-\(N\) certificates verify two distinct positive-curvature roots, a nonzero cost slope, and strict margins against every other regular or confluent competitor. The tangent limiting-chart exhaustion depends on the fitting class, which is unchanged by the shifted target. On the compact finite frequency quotient, choose disjoint neighborhoods of A and B; the complement has a positive cost gap at the baseline. The perturbation \(4\kappa^2P_N\) is continuous and bounded, so this gap and the coalescent gap persist for some **existential** \(\kappa_0(N)>0\). The implicit-function theorem continues the roots and unique transverse equality. Hence, for each certified finite \(N\) and \(|\kappa|<\kappa_0(N)\), the global noncoalescent small-\(z\) exchange persists with

\[
\epsilon_c(z,\kappa)=\lambda_N(\kappa)z^2+O(z^4).
\]

The base large-\(N\) theorem supplies the same pointwise conclusion for each covered odd \(N\ge10001\); no \(N\)-uniform phase radius is established. The continuum full-line certificate also has a strict compact-region, coalescent, and tail margin; \(P_\infty\le1/12\) keeps the tail stable for sufficiently small \(|\kappa|\). Thus an existential continuum \(\kappa_0(\infty)>0\) exists. These are deductions from the frozen certified base cases, **not** interval certification of any displayed nonzero \(\kappa\) value.

## Crossing curvature at the aligned case

Let \(\Gamma_N=\partial_\lambda(R_A-R_B)>0\) at \(\kappa=0\). Envelope differentiation gives

\[
\lambda_N(\kappa)=\lambda_N(0)
-\frac{4\{P_N(v_A)-P_N(v_B)\}}{\Gamma_N}\kappa^2
+O(\kappa^4),\qquad \lambda_N'(0)=0.
\]

For \(N=21,b=10\), the quadratic coefficient is \(1.92838866512018\ldots\). In the continuum, use \(D_\infty(v)=\int_{-1/2}^{1/2}e^{ivs}\,ds=\operatorname{sinc}(v/2)\), \(S_2=1/12\), and
\(P_\infty(v)=1/12-D_\infty'(v)^2/(1-D_\infty(v)^2)\). The independently certified aligned continuum crossing is \(\lambda_\infty(0)=0.0665069703923955\ldots\); its quadratic coefficient is \(1.94417110591517\ldots\). These displayed coefficients are high-precision evaluations of an analytic envelope formula, not new interval enclosures.

## Requested N21 map and its proof boundary

`kappa_crossing_map.csv` records the stationary A/B solve, curvatures, slope, separation, coalescent gap, and a full-period **numerical** third-minimum check for \(\kappa=0,\pm0.05,\pm0.10,\pm0.20,\pm0.30,\pm0.50\). The map is exactly even up to rounding. All requested sampled values have positive local curvatures, positive transverse slopes, positive coalescent gaps, and A/B as the sampled best pair. Selected positive-\(\kappa\) values are:

| \(\kappa\) | \(\lambda_c\) | \(v_A\) | \(v_B\) | Coalescent gap | Sampled third-regular gap |
|---:|---:|---:|---:|---:|---:|
| 0 | 0.0659804 | -4.15344 | 12.46809 | 0.0180476 | 0.0644497 |
| 0.05 | 0.0709173 | -3.49100 | 12.34843 | 0.0153244 | 0.0746620 |
| 0.10 | 0.0855110 | -2.39370 | 12.07072 | 0.0111358 | 0.109630 |
| 0.20 | 0.132226 | -1.04357 | 11.56045 | 0.00506306 | 0.272019 |
| 0.30 | 0.188224 | -0.365517 | 11.24974 | 0.00120220 | 0.570991 |
| 0.50 | 0.308064 | 0.307574 | 10.92823 | 0.00212978 | 1.59954 |

The A branch approaches the confluent coordinate \(v=0\) between \(\kappa=0.38\) and \(0.39\). Solving the **numerical** conditions \(\lambda=V_2/Q\), B-stationarity, and equality of B with the confluent cost gives \(|\kappa|\approx0.3882278594\) at \(N=21\); the continuum analogue is \(0.4000361843\). There the strict C7 margin vanishes, so the certified noncoalescent A/B theorem cannot be continued across that point without a different chart. The positive sampled margin at \(\kappa=0.50\) does not repair this connected-theorem break. No explicit interval-validated lower bound on \(\kappa_0\) is claimed.
