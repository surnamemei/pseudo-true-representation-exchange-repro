# All-sufficiently-large-odd-$N$ persistence theorem

**Stage 15 effective update.** The existential $N_0$ below has now been made explicit by a uniform interval-assisted certificate: every odd $N\ge35{,}377$ is covered for the specified $b,\phi$ family. See `explicit_N0_derivation.md` and `explicit_N0_certificate.json`. The text below records the original Stage 14 existence argument.

**Theorem (all sufficiently large odd records, specified family).** Fix $b=10$ and $\phi=\pi$. There exists an integer $N_0$ such that every odd $N\ge N_0$ satisfies the finite-$N$ tangent conditions C1–C7 at a unique transverse A/B crossing $\lambda_N$ near $\lambda_\infty=0.0665069703923955\ldots$. The existing limiting-chart exhaustion C8 applies to each such $N$. Therefore each odd $N\ge N_0$ obeys the established small-$z$ global noncoalescent branch-exchange theorem, with a possibly $N$-dependent positive small-$z$ radius. The theorem does not assign a numerical value to $N_0$.

**Certified continuum premise.** The local interval Krawczyk certificate verifies two distinct regular roots, positive Gram factors and curvatures, nonzero satellite coefficients, and a transverse equal-cost slope. The full-line interval partition/replay, coalescent bound, and analytic $|v|\ge100$ tail bound exclude every other global competitor; see continuum_global_certificate.md and its JSON/replay objects. The third-party numerical scan is not used to establish this premise.

**Proof.** Choose disjoint closed frequency wells around the two continuum minima, avoiding zero. Positive Gram and curvature bounds persist on these wells. The continuum global gap gives (i) a strictly higher cost on the compact complement of the wells and an endpoint neighborhood and (ii) a strictly higher tail cost beyond some fixed $V$. On the compact regular part, centered midpoint quadrature gives $C^2$ convergence $R_N/N\to R_\infty$ at order $N^{-2}$, uniformly in $\lambda$ near $\lambda_\infty$. Analytic deflation of $A/v^4$ and $B/v^2$ gives uniform convergence through the endpoint, preserving its positive gap.

The expanding finite frequency circle requires a separate argument. On the principal period $|v|\le\pi N$, discrete summation by parts bounds $|d_N^{(k)}(v)|$ by $2\pi 2^{-k}/|v|$ for $k=0,1,2$. The corresponding $b-v$ bound uses $|v|-b$. Thus the satellite projection gain is uniformly small for $|v|\ge V$, independently of sufficiently large $N$. No fitted satellite can escape to a high-frequency alias while competing globally.

The implicit-function theorem now continues the two continuum stationary wells for all large odd $N$. Their positive curvature and nonzero coefficients persist. The nonzero cost-difference slope gives one transverse zero $\lambda_N$. Compact, endpoint, and tail gaps exclude all other global competitors at that zero and nearby. Conditions C1–C7 follow, and the finite-$N$ compactification/scaling/recurrence lemmas already establish C8. Applying the existing finite-$N$ theorem completes the proof.

The strict continuum predicates are verified by continuum_global_certificate.json (verified_global=true), stage14_continuum_global_replay.json (replay_verified=true), and stage14_continuum_local_interval.json (verified_local=true). No numerical $N_0$ is claimed until a uniform quantitative replay supplies one.

## Analytic crossing expansion

With $t=N^{-2}$, the exact kernel $d_N(q)=d(q)\{1+q^2t/24+7q^4t^2/5760+O(t^3)\}$ is analytic near $t=0$ on the A/B wells. Analytic implicit continuation therefore gives

$$
\lambda_N=\lambda_\infty+aN^{-2}+bN^{-4}+O(N^{-6}).
$$

The coefficients follow by differentiating the stationary and equal-cost equations, as detailed in lambdaN_asymptotic_fit.md. Their high-precision values are $a=-0.2315601086498066\ldots$ and $b=-0.2831432749096778\ldots$. Coefficient interval enclosures and an effective $N_0$ are separate verification tasks; the existence theorem does not require them.
