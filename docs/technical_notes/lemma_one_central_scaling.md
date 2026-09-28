# Lemma: uniform scaling with one central node

Let N be odd, N>=3, s=t/N on N consecutive centered integer samples, and e_v=(exp(iv s_t))_t. Frequencies are identified modulo 2*pi*N. Fix v_* not congruent to zero and choose a closed frequency neighborhood V of v_* avoiding zero. Shrink V if necessary.

Define the five-real-dimensional parameter vector

\[\theta=(\Re(a-2),\Im(a-2),u,\Re b,\Im b),\quad F(\theta,v)=a e_u+b e_v.\]

Thus F(0,v)=2 and

\[L_v=D_\theta F(0,v)=[1,i1,2is,e_v,ie_v]\]

is a real linear map into C^N with its Euclidean real norm. These columns are independent. Indeed, if A+B e_v+C s=0, taking second differences gives B*(exp(iv/N)-1)^2*e_v=0; hence B=0, and A=C=0. This proves independence even with complex C, and therefore for the restricted purely imaginary coefficient in L_v. The argument needs N>=3.

Continuity on V gives sigma=min_{v in V} sigma_min(L_v)>0. Smoothness supplies M<infinity and r>0, uniform in v, such that

\[\|F(\theta,v)-2-L_v\theta\|\le M\|\theta\|^2,\qquad \|\theta\|\le r.\]

Taking r<=sigma/(2M),

\[\|\theta\|\le (2/\sigma)\|F(\theta,v)-2\|
\le (2/\sigma)(\|F(\theta,v)-x_p\|+\|x_p-2\|).\]

This is the required local left-inverse estimate. A residual-only bound relative to x_p would be false without accounting for the O(p) perturbation of x_p.

Suppose x_p=2+p q+O(p^2), u_p->0, v_p->v_*, and ||F(theta_p,v_p)-x_p||=O(p). The two-column Gram matrix stays uniformly invertible; convergence of the fitted vectors to 2 first implies a_p->2 and b_p->0. Hence theta_p eventually enters the uniform neighborhood above and theta_p=O(p). Along a subsequence,

\[(F-2)/p\longrightarrow \alpha+2i\kappa s+\beta e_{v_*},\qquad
(x_p-F)/p\longrightarrow q-\alpha-2i\kappa s-\beta e_{v_*}.\]

No rate for v_p-v_* is required: its coefficient b_p is already O(p). This distinction is essential to the limiting-chart argument.
