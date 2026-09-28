# Lemma: both nodes central, without amplitude bounds

Let p_j->0+, and let y_j lie in an ordinary or confluent two-tone plane with physical nodes xi_{1j},xi_{2j} on the unit circle, both tending to 1. Assume

\[y_j=2+p_j h_j,\qquad \sup_j\|h_j\|<\infty.\]

Choose normalized frequency lifts u_{kj}->0 with xi_{kj}=exp(i u_{kj}/N). Every such fitted vector obeys the exact recurrence

\[y_{t+2}-(\xi_1+\xi_2)y_{t+1}+\xi_1\xi_2 y_t=0.\]

The same recurrence holds for a repeated node and its confluent vector t*xi^t. Substitution yields

\[2(1-\xi_1)(1-\xi_2)/p+
 h_{t+2}-(\xi_1+\xi_2)h_{t+1}+\xi_1\xi_2h_t=0.\]

Since N>=3, at least one recurrence equation is available. Bounded h implies (1-xi_1)(1-xi_2)/p=O(1). The identity

\[(1-e^{iu_1/N})(1-e^{iu_2/N})
=-\frac{u_1u_2}{N^2}(1+o(1))\]

holds uniformly for local lifts, including zero factors. Consequently u_1u_2/p is bounded and real. Extract a subsequence with u_1u_2/p->gamma in R and h_j->h. Passing to the limit gives

\[\Delta_t^2 h_t=2\gamma/N^2,\qquad
h_t=\alpha+\beta s_t+\gamma s_t^2,\quad \alpha,\beta\in\mathbb C.\]

This characterizes a containing class of normalized **fitted perturbations**, not the residual itself. If x_p=2+p q+O(p^2), the limiting residual is q-h and its loss is at least

\[C_{\rm coal}=\min_{\alpha,\beta\in\mathbb C,\gamma\in\mathbb R}
\|q-\alpha-\beta s-\gamma s^2\|^2.\]

The argument imposes no relative approach rate, no bounded-amplitude assumption, and no exclusion of divergent cancellation between the two fitted coefficients. Not every element of the containing class must be attainable: a lower bound over a larger class is sufficient for exclusion. For the reported phase-pi instances its optimizer has negative gamma and is attainable by frequencies ±sqrt(-gamma)*sqrt(p) and suitable amplitude corrections.

Finally, for any O(p)-residual sequence, the limiting recurrence of y->2 implies (1-xi_1)(1-xi_2)=0. Thus at least one limiting node is central; there is no fourth chart with both limiting nodes noncentral.
