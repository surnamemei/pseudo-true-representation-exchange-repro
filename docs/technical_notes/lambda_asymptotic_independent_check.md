# Independent check of the large-record crossing coefficient

**Result: PASS after a precision repair in the validation generator.** An independent SymPy coefficient algebra and implicit-root calculation agrees with the corrected multiprecision calculation, without fitting finite-(N) data.

For centered odd (N), $s_n=n/N$ and $t=N^{-2}$. Exact midpoint summation gives

\[
d_t(q)=\frac{\sin(q/2)}{N\sin(q/(2N))}
=d(q)+t\frac{q^2d(q)}{24}+t^2\frac{7q^4d(q)}{5760}+O$t^3$,
\quad d(q)=\frac{2\sin(q/2)}q.
\]

The independent symbolic calculation expands the reduced tangent loss as $R_0(v,\lambda)+tR_1(v,\lambda)+t^2R_2(v,\lambda)$, including $\mu_{2,t}=(1-t)/12$ and $\mu_{4,t}=1/80-t/24+7t^2/240$. It then perturbs the **joint stationary/equal-cost system**, not just the crossing equation at fixed frequencies. If $S=\partial_\lambda(R_{0,A}-R_{0,B})$, $H_j=R_{0,vv,j}$, and $G_j=R_{0,v\lambda,j}a_2+R_{1,v,j}$, then

\[
a_2=-\frac{R_{1,A}-R_{1,B}}{S},\qquad
a_4=-\frac{\tfrac12 F_{\lambda\lambda}a_2^2+F_{\lambda t}a_2+$R_{2,A}-R_{2,B}$-G_A^2/$2H_A$+G_B^2/$2H_B$}{S}.
\]

The branch-frequency drifts $-G_j/H_j$ are included. The derived values are

\[
\begin{aligned}
\lambda_\infty&=0.0665069703923955354559\ldots,\\
a_2&=-0.2315601086498065983416\ldots,\\
a_4&=-0.2831432749096774046231\ldots,\\
v_{A,2}&=-40.495636462185502392\ldots,\\
v_{B,2}&=19.158387625772463283\ldots .
\end{aligned}
\]

Here $S=0.10781751364959401460\ldots>0$, and both continuum curvatures are positive. The negative signs of $a_2,a_4$ are therefore consistent with the finite lengths approaching $\lambda_\infty$ from below. The expansion convention is $N^{-2}$ for the **actual odd record length (N)**; no $N\pm1$ substitution is used. Centered symmetry removes odd inverse powers.

The earlier script `stage14_asymptotics.py` contained two floating-point ingress paths: `(1-t)/12` evaluated with integer `t=0` produced binary `1/12`, and the finite kernel evaluated `10/(2*N)` with integer arguments before conversion to multiprecision. The first changed archived coefficient decimals only beyond about 15 digits. The second spoiled the $N^6$-scaled validation residual at large (N). Both are repaired by entering multiprecision before division; the independent SymPy calculation matches the repaired coefficients. The manuscript's 11-digit displayed coefficients and all five pre-existing finite-(N) certificates were unaffected. The revised numerical validation rows and figures are regenerated from the corrected helper.

`stage16_residual_scaling.csv` reports the requested nine lengths, high-precision $\lambda_N$, the continuum-only prediction, and predictions through $N^{-2}$ and $N^{-4}$. The first five input $\lambda_N$ values are from existing interval-certified anchors; $N=101,201,501,1001$ are numerical root solves only. In the corrected table, $N^4(\lambda_N-\lambda_\infty-a_2/N^2)$ moves from (-0.3011078) at $N=11$ to (-0.28314526) at $N=1001$, approaching $a_4=-0.28314327\ldots$. After the quartic correction, $N^6$ times the residual moves from (-2.17371) to (-1.99101), with the $N=501$ row at (-1.99107). This is a numerical scaling check, not an additional interval theorem. The coefficient derivation is recorded in `stage16_symbolic_asymptotics.py` and its JSON output.
