# Theorem: checkable finite-N global branch exchange

Fix odd N>=3 and (b,phi). Suppose C1–C8 in `generalN_conditions.md` hold at lambda_N>0. Then there exist p_0,eta>0 and two analytic strict local-minimum branches A/B of the divided loss such that, for 0<p<p_0 and |lambda-lambda_N|<eta:

1. The global compactified two-tone minimum is attained by an ordinary separated fit. Up to permutation and periodic identification, every global fit is A or B.
2. Their central frequency is O(p); their satellite frequencies tend to the distinct regular roots v_A,v_B; the satellite amplitudes are p*beta_j+O(p^2) and nonzero. Their full ordinary-coordinate Hessians are positive definite for each p>0.
3. Exactly one analytic transverse equal-cost curve lambda_c(p) passes through lambda_N. At equality precisely A and B are global. On either side the lower-cost branch is the unique global fit.
4. With p=z^2,

\[\epsilon_c(z)=\lambda_Nz^2+c_Nz^4+O(z^6),\quad
c_N=-\frac{T_A-T_B}{\partial_\lambda(R_A-R_B)}.\]

If the denominator is positive, A wins below and B above. The existence radii and remainder constant are not numerical bounds.

## Proof

**Local branches.** Write y=(2+alpha*p)exp(i*kappa*p*s)+p*beta*exp(i*v*s). The divided residual (x-y)/p has a removable analytic extension to p=0. At each tangent minimum the five linear-coordinate Hessian is twice its positive Gram matrix. Its Schur complement is R_vv>0. The six-real-coordinate Hessian is therefore positive definite. The analytic implicit-function theorem produces unique stationary minima near the two tangent solutions. C3 and regular nonzero satellite frequencies make the map to ordinary two-tone coordinates invertible for p>0; Hessian definiteness is preserved by congruence. Applying the implicit-function theorem to their cost difference using C5 gives lambda_c(p).

**Existence and exhaustion.** The compactified-subspace lemma gives a minimizing plane for every p. A/B trial fits have loss O(p^2); its projected fitted vector y obeys y=2+O(p). Any sequence of minimizing planes has limiting unit-circle nodes and bounded h=(y-2)/p. The recurrence forces at least one central limiting node. Exactly one central node gives, by the uniform scaling lemma, a central-plus-satellite tangent fit; two central nodes give the containing quadratic class from the recurrence lemma.

Suppose global minimizers escaped both tangent wells for a sequence p->0 and lambda near lambda_N. For both central nodes the limiting divided cost is at least C_coal, contradicting C7. For a regular satellite, optimizing the five linear coordinates bounds the limiting cost below by R_N. Complete scalar tangent classification and C6–C7 exclude all limits outside A/B neighborhoods. Shrink eta: continuity preserves this strict exclusion, including the central endpoint. The uniform scaling estimate makes the remaining scaled fit coordinates bounded; their limits must be the unique optimal tangent coefficients. They therefore enter the six-dimensional neighborhoods where the local implicit-function theorem gives the unique minimum branches. Thus the minimizing plane is nonconfluent and the ordinary minimum is attained. This proves globality and uniqueness away from equality.

**Coefficient.** At fixed tangent coordinates,

\[r_0=q-\alpha-2i\kappa s-\beta e^{ivs},\qquad
r_1=s^4/12-i\kappa\alpha s+\kappa^2s^2,\]
\[J_j=p^2[R_j+pT_j+O(p^2)],\qquad T_j=2\Re\sum\overline{r_{0,j}}r_{1,j}.\]

Stationarity cancels coordinate-derivative terms in the first optimized p derivative. Implicit differentiation of the equal-cost equation gives c_N as displayed. Analyticity supplies the uniform local remainder and
J_A-J_B=p*(epsilon-epsilon_c)*[S+O(p+|lambda-lambda_N|)].

## Corollary and evidence

**Corollary.** For b=10 and phi=pi, C1–C7 are interval-certified for every N in {11,15,21,31,41}; C8 follows from the lemmas. Thus the theorem applies to all five instances.

The verified subset and interval bounds are recorded in `generalN_certified_instances.csv`; full finite-sum root, partition and margin certificates are in `stage12_generalN`. C1–C7 are computational finite-dimensional checks. C8 is proved by the three lemmas, rather than being assumed from a numerical sweep. No conclusion for untested N or a connected bridge to z=2 follows from this theorem.
