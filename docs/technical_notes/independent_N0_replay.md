# Stage 16 independent replay of the effective all-large-odd-N certificate

**Result: PASS.** A new python-flint/Arb implementation verified the uniform parameter box $0\le t=N^{-2}\le 35{,}377^{-2}$, including $t=0$. The frozen compact partition contains 27,708 cells; all pass. Ten original cells required Arb-only subdivision (26 refinement nodes); none failed. This is an independent arithmetic and evaluator replay of the frozen obligations, not a new parameter search.

## Independent implementation and scope

Run `stage16_independent_arb.py` with python-flint 0.9.0. It imports none of the Stage 14/15 objective or interval functions. It constructs its own Arb kernel, second-order dual derivatives, reduced loss, Krawczyk map, zero-frequency Taylor chart, coalescent bound, and tail bound. The only reused inputs are the candidate center and frozen cell coordinates/predicate labels. `independent_N0_replay.csv` gives a decision for every cell; `independent_N0_replay.json` gives exact counts and directed margins. The script uses 85 decimal digits of Arb midpoint precision.

The 27,708 frozen cells exactly cover ([-100,100]) with matching rational-decimal endpoints. Their labels comprise 25,274 cost exclusions, 1,420 gradient exclusions, 1,012 monotone no-root exclusions, and one root cell for each branch. The Arb replay passed 27,708, refined 10 original cells using 26 subdivision nodes, and failed 0. The root cells contain the parameterized Krawczyk roots. A subdivision of the wider B root cell certifies positive curvature throughout while opposite endpoint-gradient signs give a unique root. Other refined cells are excluded by cost, derivative sign, or monotonicity on their subcells. The least original-width passing Arb cost-cell margin is $4.0740\times10^{-9}$; gradient and curvature-sign margins exceed $1.03\times10^{-8}$ and $1.02\times10^{-7}$, respectively. The refined subcells are checked separately, so these three minima are not claims about the initial enclosure of every frozen cell.

## Root persistence and signs

The parameterized three-variable Krawczyk map lies strictly inside both branch-frequency boxes and the crossing-coefficient box. Its infinity-norm contraction bound is $0.054366<1$. Over the whole parameter interval, the Arb lower bounds are:

| Quantity | Directed lower bound |
|---|---:|
| A reduced curvature | $1.3030081160\times10^{-4}$ |
| B reduced curvature | $9.6031586246\times10^{-4}$ |
| A satellite coefficient | (0.1274886219) |
| magnitude of B satellite coefficient | (0.0644223989) |
| A-minus-B equal-cost slope | (0.1078144621) |

The two frequency boxes remain disjoint. The signed satellite coefficients stay away from zero; positive reduced curvature and the one-dimensional root-cell tests exclude coalescent or flat competing representations in the regular chart. Since $t=0$ is included, these checks independently enclose the continuum root and its curvatures as well as every covered finite odd length. The interval inclusion proves one joint root in the stated box for each (t); it does not assert an (N)-independent small-(z) radius.

## Finite kernel and monotone remainder

Independently expand $H(w)=\sqrt w/\sin\sqrt w$ by multiplying it with $\sin\sqrt w/\sqrt w$; the rational recurrence gives coefficients through $w^8$. For $|q|\le110$, use a complex (q)-circle of radius 55. The Cauchy bound for every (q) derivative through order 46 is controlled by

\[
 B(N)=\frac{e^{55/2}}{1-\eta}\frac{z_N^9}{1-z_N},\qquad
 z_N=(165/N)^2,\quad \eta=\frac{\sinh(1/2)}{1/2}-1.
\]

Indeed $46!/55^{46}<1$; the (q)-circle reaches modulus at most 165. Arb gives $B(35{,}377)=9.98991003324\times10^{-31}<10^{-30}$ and $B(35{,}375)=1.00000813551\times10^{-30}>10^{-30}$. Since (B(N)) decreases with $N>165$, the chosen $10^{-30}$ envelope applies uniformly to all odd $N\ge35{,}377$. The previous odd length only misses this *chosen envelope*; it is not a failure of the branch exchange.

For the zero chart, the sinc Maclaurin tail after index 37 is below $10^{-50}$ on $|q|\le1$ through derivative order 46 and on $|q|\le2$ for the derivatives through order four used there. The deflated $A(v)/v^4$ and $B(v)/v^2$ Taylor polynomials retain degrees through 40. On the complex (v)-circle of radius 10, the exact parameterized kernel and its first two derivatives bound the deflated functions by $10^7$, using $|d_t(q)|\le e^{|\operatorname{Im}q|/2}/(1-\eta_{20})$, Cauchy derivative bounds on unit circles, $|\lambda|<0.067$, and $\mu_{2,t}>0.083$. The resulting value/first/second derivative tail on $|v|\le1$ is below $3\times10^{-31}$; the propagated kernel truncation is below $10^{-25}$. This validates the conservative $10^{-25}$ zero-chart enclosure used by the Arb checker.

## Coalescent class and expanding tail

The analytic coalescent extension is evaluated with $\mu_{2,t}=(1-t)/12$, $\mu_{4,t}=1/80-t/24+7t^2/240$, and the deflated quadratic coefficient. Its Arb cost gap above the feasible A incumbent is at least (0.0008511669596), uniformly in (t). At $t=0$, the independently enclosed limiting gap remains positive.

For odd finite (N), represent the frequency on the principal period $|v|\le\pi N$. Abel summation for $N^{-1}\sum s_n^k e^{iv s_n}$ yields the conservative envelopes $|d_N^{(k)}(v)|\le 2\pi 2^{-k}/|v|$ for $k=0,1,2$. For the shifted kernel $d_N(10-v)$, its wrapped distance on the principal period is at least (|v|-10); hence $|d_N(10-v)|\le2\pi/(|v|-10)$. These bounds hold for the full $|v|\ge100$ part of the principal period, including the endpoint strip where $|10-v|>\pi N$. Substituting into the positive Gram denominator and reduced loss gives an Arb tail gap of at least (0.0031356600897) above the same incumbent for every covered (N). This explicitly closes the compact-partition/tail join at $|v|=100$.

## Interpretation

The independent replay supports the specified $b=10,\phi=\pi$ all-large-odd-(N) small-spacing theorem. **35,377 is a conservative sufficient certification threshold from uniform bounds.** Direct separate certificates already cover $N=11,15,21,31,41$. Neither implementation identifies the smallest true length, gives a uniform positive small-spacing radius, or connects the local theorem continuously to the isolated $z=2$ certificates.
