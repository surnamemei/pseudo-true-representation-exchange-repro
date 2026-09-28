# Fixed nonzero strong-pair phase: first-order regime

**Model.** On centered samples (s_t=t/N), (N=21), take

\[
x_t=2\cos(zs_t+\delta)-\epsilon e^{ibs_t},\qquad b=10,
\]

with fixed (0<|\delta|<\pi/2) as (z\to0). The range excludes the separate zero-leading-constant case \(\cos\delta=0\). Expanding the pair gives

\[
2\cos(zs+\delta)=2c-2zds-z^2cs^2+\tfrac13z^3ds^3+O(z^4),
\quad c=\cos\delta,\ d=\sin\delta.
\]

## Why one central exponential cannot absorb the linear term

If a single central fitted exponential has amplitude \(2c+z\alpha\) and real frequency \(z\xi\), its first-order perturbation is \(\alpha+2ic\xi s\). Its frequency direction is **imaginary** linear in \(s\), whereas the strong pair supplies **real** \(-2ds\). Since \(c\ne0\) and \(d\ne0\), no choice of complex amplitude perturbation and real \(\xi\) represents that term. A noncentral satellite must help, or both fitted frequencies must approach the center. The exact two-strong-tone fit at \(\epsilon=0\) belongs to the latter class.

## Correct competing tangent charts and scale

Set \(\epsilon=\lambda z\). For one central fitted tone plus a noncentral satellite at limiting normalized frequency \(v\ne0\), the first-order fit is

\[
y=2c+z\{\alpha+2ic\xi s+\beta e^{ivs}\}+O(z^2).
\]

Its residual target is \(q=-2ds-\lambda e^{ibs}\). In the two-central closure, both fitted frequencies tend to zero and the first-order tangent space contains complex \(1,s\) and a real \(s^2\) coefficient. The latter is realizable by splits of order \(\sqrt z\): leading real amplitudes summing to \(2c\), with vanishing first moment, yield a quadratic coefficient \(c\xi_1\xi_2\), of either sign. Complex amplitude corrections of order \(\sqrt z\) produce the required real-linear term. Thus the coalescent loss is \(\lambda^2 K_N\), where \(K_N\) is the residual weak-tone energy after projection onto that class; the strong-pair linear term is absorbed exactly.

If \(\epsilon=o(z)\), the exact strong pair fitted at order two has loss \(O(\epsilon^2)=o(z^2)\), while a fit with a satellite bounded away from the center cannot represent the real-linear strong-pair term at order \(z\). If \(z=o(\epsilon)\), a central-plus-weak-tone fit has loss \(O(z^2)=o(\epsilon^2)\), while any two-central fit leaves nonzero weak-tone projection error of order \(\epsilon^2\). The competition between these representation classes therefore occurs at **\(\epsilon=\Theta(z)\)**. This establishes the exponent \(\alpha=1\) as the natural fixed-\(\delta\) scale; it does not alone certify a unique crossing coefficient.

## Exact scalar tangent reduction

Let \(D_N(v)=\sum_t e^{ivs_t}\), \(S_k=\sum_t s_t^k\),
\(G_0(v)=N-D_N(v)^2/N\), and

\[
P_N(v)=S_2-\frac{D_N'(v)^2}{G_0(v)}\ge0.
\]

Let \(H_N(v)\) be the optimized regular-chart squared loss for target \(-e^{ibs}\) over the real tangent span \(\{1,i,is,e^{ivs},ie^{ivs}\}\). Real-parity blocks are orthogonal, so the fixed-phase tangent loss is **exactly**

\[
R_{N,\delta}^{(1)}(v,\lambda)
=\lambda^2H_N(v)+4\sin^2\delta\,P_N(v).
\]

Explicitly, with \(A=N-D(v)^2/N-D'(v)^2/S_2\) and \(B_w=-D(b-v)+D(v)D(b)/N+D'(v)D'(b)/S_2\),
\[
H_N(v)=N-D(b)^2/N-D'(b)^2/S_2-B_w^2/A.
\]
Writing \(V_2=S_4-S_2^2/N\), the confluent coefficient is
\[
K_N=N-D(b)^2/N-D'(b)^2/S_2-
\{D''(b)+S_2D(b)/N\}^2/V_2.
\]

For \(N=21,b=10\), the coalescent cost is \(\lambda^2K_{21}\), with \(K_{21}=17.64367452454865\ldots\). For \(\sin\delta\ne0\), divide by \(4\sin^2\delta\) and put \(\mu=(\lambda/(2|\sin\delta|))^2\). The universal normalized scalar problem is \(P_{21}(v)+\mu H_{21}(v)\). This formula was cross-checked against direct real least squares.

## N21 numerical diagnosis at five fixed phases

A full-period numerical scan and high-precision stationary/equal-cost solve find a candidate two-well crossing at

\[
\mu_c=0.09842255880858805\ldots,\quad
v_A=1.618255087376875\ldots,\quad v_B=10.304770817765262\ldots .
\]

Both scalar curvatures are positive; the normalized coalescent margin is \(0.0526083799168584\ldots\), and the sampled third-regular-minimum margin is \(1.8349580\ldots\). The candidate critical ratio is

\[
\frac{\epsilon_c}{z}=2\sqrt{\mu_c}|\sin\delta|
=0.6274473963882169\ldots |\sin\delta|.
\]

| \(\delta/\pi\) | Candidate \(\epsilon_c/z\) | Coalescent gap | Sampled third-regular gap |
|---:|---:|---:|---:|
| 0.05 | 0.0981544 | 0.00514967 | 0.179618 |
| 0.10 | 0.193892 | 0.0200946 | 0.700892 |
| 0.20 | 0.368804 | 0.0727030 | 2.53585 |
| 0.30 | 0.507616 | 0.137731 | 4.80398 |
| 0.40 | 0.596738 | 0.190339 | 6.63894 |

The same stationary coordinates appear at all five phases because the normalized tangent problem is independent of \(\delta\). `fixed_delta_phase_map.csv` records costs and curvatures. **The crossing and full-domain ordering in this table are numerical, not interval-certified.** No finite-width or all-\(z\) statement follows. As \(\delta\to0\) with \(z\) fixed, this first-order chart becomes nonuniform; that is why the \(\delta=\kappa z\) regime must be analyzed separately.
