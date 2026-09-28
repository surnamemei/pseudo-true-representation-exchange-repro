# Checkable conditions for global finite-N branch exchange

Fix odd N>=3, normalized weak location b modulo 2*pi*N, and phase phi. Put s=t/N, p=z^2, epsilon=lambda*p, q=-s^2+lambda*exp(i*phi)*exp(i*b*s). Norms are unnormalized sample Euclidean norms. Define R_N by eliminating alpha,beta in C and kappa in R from ||q-alpha-2i*kappa*s-beta*exp(i*v*s)||^2.

All following conditions concern a fixed lambda_N and two roots on the frequency circle, distinct modulo 2*pi*N. They are not asserted for every N,b,phi.

|Condition|Required statement|Verification|
|---|---|---|
|C1|Two regular stationary satellites v_A,v_B, each nonzero modulo the period; real five-column tangent Gram positive definite.|Joint stationary/equal-cost interval root; finite-sum Gram bounds.|
|C2|R_vv(v_A),R_vv(v_B)>0.|Directed interval curvature bounds.|
|C3|Optimized complex beta_A,beta_B are nonzero.|An interval component excludes zero or a positive lower norm bound.|
|C4|R_A=R_B at lambda_N>0.|A three-variable Krawczyk inclusion for (R_v(A),R_v(B),R_A-R_B); numerical equality alone is insufficient.|
|C5|S=partial_lambda(R_A-R_B)!=0.|Envelope derivatives on the root box. Positive S means A below, B above.|
|C6|Every other **regular** tangent minimum has cost at least R_*+delta_reg, delta_reg>0.|Complete periodic stationary isolation and strict cost comparisons; no unexamined interval. This is not a uniform gap arbitrarily close to A/B.|
|C7|C_coal(lambda_N)>R_*+delta_coal, delta_coal>0.|Projection onto complex{1,s}+real{s^2}. This includes the strong-pair limiting fits and the v->0 endpoint, which is NOT counted in C6.|
|C8|All O(p)-residual competitors reduce to the one-central-node chart or the both-central-node class; closure minimizers exist.|The three formal lemmas prove this for the present undamped, unrestricted-complex-amplitude finite-record family, N>=3.|

For fixed odd N>=3 the regular Gram has full real rank off v=0. The scalar loss extends analytically through v=0 after deflation; its endpoint limit is the coalescent projection. Thus complete stationary isolation and C7 provide a strictly positive cost gap outside fixed small A/B neighborhoods by compactness. C6 is used with this analytic compactification, not as a claim that all nonstationary points have the stationary margin.

## General finite-sum formula

Let D(v)=sum exp(iv s), S_k=sum s^k, c=cos(phi), d=sin(phi). Define

\[q_0=-S_2+\lambda c D(b),\quad q_1=-2\lambda c D'(b),\quad
q_3=D''(v)+\lambda c D(b-v),\]
\[A=N-D(v)^2/N-D'(v)^2/S_2,\quad
B=q_3-D(v)q_0/N+D'(v)q_1/(2S_2),\]
\[C=N-D(v)^2/N,\quad H=\lambda d[D(b-v)-D(v)D(b)/N].\]
Then, for regular v,

\[R_N=S_4+2\lambda cD''(b)+N\lambda^2-q_0^2/N-q_1^2/(4S_2)
-\lambda^2d^2D(b)^2/N-B^2/A-H^2/C.\]

The fitted satellite coefficient is beta=B/A+i*H/C. Its central coefficient is
alpha=(q_0-D(v)*Re(beta))/N+i*(lambda*d*D(b)-D(v)*Im(beta))/N,
and kappa=(q_1+2D'(v)*Re(beta))/(4S_2).

An explicit general-phase coalescent loss is

\[C_{\rm coal}=\|q\|^2-\frac{|\sum q|^2}{N}
-\frac{|\sum s q|^2}{S_2}
-\frac{[\Re\sum(s^2-S_2/N)q]^2}{S_4-S_2^2/N}.\]

Here sum q=-S_2+lambda*exp(i phi)*D(b), sum s q=-i lambda*exp(i phi)*D'(b). At phi=pi this reduces to the expression in the manuscript. Taking v->0 in the regular five-column projection spans precisely this real five-dimensional class; deflate A/v^4, B/v^2, C/v^2 and H/v before evaluation.
