"""Assemble the major-revision reports from the recorded calculations."""
import csv,json,hashlib
from pathlib import Path
from decimal import Decimal,getcontext
ROOT=Path(__file__).resolve().parent;getcontext().prec=65
def read(n):return json.loads((ROOT/n).read_text())
def write(n,s):(ROOT/n).write_text(s.strip()+'\n',encoding='utf-8')
def table(rows,keys):
    return '| '+' | '.join(keys)+' |\n| '+' | '.join(['---']*len(keys))+' |\n'+'\n'.join('| '+' | '.join(str(r[k]) for k in keys)+' |' for r in rows)
def main():
    c=read('smallz_global_certificate.json');p=read('stage11_phase_validation.json');f=read('finite_epsilon_certificate.json')
    tref=read('stage3_tangent_reference.json');base=read('stage2_reference.json')
    pred=4*Decimal(tref['lambda_N'])+16*Decimal(tref['lambda_quartic']);err=Decimal(base['epsilon_cross'])-pred
    rt=[dict(v=f"{r['v']:.10f}",type=r['kind'],cost=f"{float(r['R'][0]):.12f}") for r in c['stationary_points']]
    write('smallz_global_tangent_classification.md',r'''
# Small-spacing global classification — N=21, b=10, phase π

## Result and proof scope

The central-plus-satellite tangent circle has **18 minima and 18 maxima**, all isolated by directed intervals. A and B are its only global minima at the jointly certified λ21. Every other stationary minimum is higher by more than 0.06444970. The best central-coalescent tangent coefficient exceeds their cost by more than 0.01804760. These are margins to other stationary minima and the coalescent stratum, **not** a uniform margin to points arbitrarily close to A/B.

The compactness argument below upgrades the local theorem to an **existential small-z global theorem**. It gives no numerical positive radius and no connected certificate reaching z=2.

## Reproducible interval enumeration

Run `python stage11_tangent_global.py`. The λ box has radius 10^-25 about the recorded 80-digit center; a new 3-variable Krawczyk check verifies the crossing within it. The cover [-66,66] contains the full period [-21π,21π]. Its small redundant end pieces contain no root. Derivative exclusion and fixed-sign second derivative plus opposite endpoint derivative signs exhaust the cover: 32,406 visited cells, zero unresolved cells. Root brackets, curvature bounds and costs are in `stage11_tangent_roots.csv`; the terminal cover is `stage11_tangent_partition.csv`.

Away from zero, all Dirichlet derivatives use finite trigonometric sums. On |v|≤1, the exact zeros are deflated before evaluating: A=v^4 Ã, B=v² B̃, R=C−B̃²/Ã. Degree-40 Taylor polynomials for Ã,B̃ and their first two derivatives have independently outward-rounded coefficients. Each omitted tail is bounded by 10^-35. A conservative exponential coefficient majorant is 10^4/k! because 2 max|s|<1; after deflation, the derivative tails are bounded by 10^4 sum_{k≥41} k²/(k+2)! <10^-43. Thus the chosen remainder is deliberately loose. No division by an interval containing zero is used in this chart.

The stationary points (display rounding only) are:

'''+table(rt,['v','type','cost'])+r'''

## The concrete strong-pair objection

Put p=z² and q=−s²−λ exp(ibs). With the fitted frequencies fixed at ±z and amplitudes optimized, the exact leading loss is

    Jstrong/z⁴ → λ² [N−D(b)²/N−D′(b)²/S2].

Indeed the strong pair itself is exactly in the fitted span, and that span converges to complex span{1,s}. This coefficient is **0.08550038682621012**, compared with A/B **0.05876262329534976**. Keeping unit amplitudes would instead give Nλ²=0.09142170770222293; that weaker comparison is not used.

Allowing both fitted frequencies to change on the z scale yields the larger tangent class

    h(s)=α+βs+γs²,  α,β complex, γ real.

Its exact minimum is

    Ccoal = λ² [N−D(b)²/N−D′(b)²/S2
                 −(D″(b)+S2 D(b)/N)²/(S4−S2²/N)]
          = 0.07681023119887829.

The optimum quadratic coefficient is γ=−0.725519235… . It is attained to leading order by frequencies ±sqrt(−γ) z, with suitable amplitude corrections supplying the constant and linear terms. Thus this is an attainable coalescent competitor, not merely a lower bound obtained by an impossible relaxation. If frequencies are constrained to ±z+o(z), γ=−1 and the fixed-strong-pair coefficient applies. Neither beats A/B.

## Exhaustion of limiting charts, including divergent amplitudes

Take any p_n→0 and two-tone fits y_n with ||x_n−y_n||=O(p_n). Then y_n=2+p_n h_n, where h_n is bounded. Extract convergent subsequences of h_n and of the two frequency nodes ξ_j=exp(iu_j/N) on the unit circle. Every y_n satisfies the exact second-order recurrence

    y[t+2]−(ξ1+ξ2)y[t+1]+ξ1ξ2 y[t]=0.

At the limit y=2 this forces (1−ξ1)(1−ξ2)=0. At least one limiting frequency is therefore central, irrespective of amplitude divergence.

**One central node, one distinct node v.** The ordinary two-column Gram is nonsingular. The central amplitude tends to 2 and the other amplitude to 0. The real derivative columns {1,i,2is,exp(ivs),i exp(ivs)} are independent for v≠0 modulo 2πN. Taylor expansion and a uniformly bounded local left inverse force the central frequency and coefficient departures to be O(p). The limit of (y−2)/p is exactly the central-plus-satellite tangent class already enumerated. Taking v=b is included; its coefficient is 0.1036226553321117. A limiting satellite amplitude equal to zero is also included.

**Both nodes central, arbitrary rates and arbitrary amplitudes.** Substitution of y=2+ph into the recurrence gives

    2(1−ξ1)(1−ξ2)/p + h[t+2]−(ξ1+ξ2)h[t+1]+ξ1ξ2h[t]=0.

The second term is bounded, so u1 u2/p is bounded, using local frequency representatives tending to zero. Its subsequential limit γ is real. Passing to the limit yields Δ_t² h=2γ/N², hence h=α+βs+γs². This proof imposes no bound on the individual fitted amplitudes, no common rate for u1,u2, and no separation lower bound. It covers near-coincident poles, the strong-pair scale, one pole approaching much faster than the other, and divergent-coefficient confluent fits. Exact confluent closures follow by continuity of the same recurrence. Both nodes converging to a noncentral frequency cannot fit 2 with error O(p).

These are all cases by compactness of the frequency torus. Their leading losses are bounded below by the two explicitly examined tangent problems.

## Global theorem and passage from leading loss

For N=21,b=10,φ=π, there exist p0,η>0 such that for 0<p<p0 and |λ−λ21|<η, every globally optimal ordinary two-tone fit is A or B (modulo permutation). Both attain the minimum at the unique analytic crossing εc(z)=λ21 z²+c21 z⁴+O(z⁶); A alone wins below it and B alone above it within this scaled neighborhood.

Proof: compactify the fitted subspaces by adjoining the confluent subspaces. Their projectors are continuous and the compactified loss attains a minimum. The A/B trials bound that loss by O(p²). If global minimizers escaped fixed disjoint neighborhoods of A/B along a sequence p→0, the recurrence classification would supply either a coalescent limit of loss at least Ccoal, or a separated tangent limit outside the A/B wells. The former contradicts the strict 0.01804760 margin. For the latter, the interval classification and compactness outside the two wells supply a strictly positive margin (possibly smaller than the displayed stationary margin), also a contradiction. Uniformity for λ near λ21 follows by continuity and by shrinking η. Inside each well the positive desingularized Hessian and the analytic implicit-function theorem give the unique local minimizer already proved. The same exclusion shows that the minimizing closure point is an ordinary separated fit. The positive crossing slope then gives the unique global switch.

The numerical interval predicates are computer-assisted ingredients; the recurrence/compactness argument is analytical. No bound on the finite-p remainder is needed for this existential theorem. Such bounds would be needed to publish an explicit p0 or a bridge to z=2.
''')
    c.update(global_smallz_theorem=True,global_theorem_scope='existential p0 and lambda neighborhood; N21,b10,phi pi',analytic_proof='smallz_global_tangent_classification.md',certified_connected_bridge_to_z2=False)
    write('smallz_global_certificate.json',json.dumps(c,indent=2))
    write('general_phase_tangent_analysis.md',r'''
# General weak-tone phase: exact tangent law and tested mechanism

For q=−s²+λ exp(iφ) exp(ibs), use the real inner product Re<a,b>. The columns [1,exp(ivs),2is] form one Gram block and [i,i exp(ivs)] the other. Centered-record parity block-diagonalizes the Gram **for every φ**; the right-hand sides, not the Gram, change with φ. Write c=cosφ, d=sinφ,

    q0=−S2+λcD(b), q1=−2λcD′(b), q3=D″(v)+λcD(b−v),
    A=N−D(v)²/N−D′(v)²/S2,
    B=q3−D(v)q0/N+D′(v)q1/(2S2),
    C=N−D(v)²/N,
    H=λd[D(b−v)−D(v)D(b)/N].

Then the exact scalar loss is

    R(v,λ,φ)=S4+2λcD″(b)+Nλ²−q0²/N−q1²/(4S2)
               −λ²d²D(b)²/N−B²/A−H²/C.

This agrees with the full five-column real least-squares elimination and reduces to the original formula at φ=π. The quantities A/v⁴, B/v², C/v² and H/v have analytic limits at v=0; numerical evaluation uses those deflated series there.

## Mechanism

The nonconstant part of the fit gain consists of two projection terms. B mixes the strong pair's real even curvature with λ cosφ through D(b−v), and can produce competing satellite wells on opposite sides of the central frequency. The quadrature component supplies H²/C, weighted by sin²φ, and reshapes those wells. At π the quadrature term vanishes and the weak component opposes the central curvature. Rotating phase changes a specific interference cross term in the finite-record projection; there is no phase-invariant orthogonality approximation that retains this mechanism.

## Results and their evidence levels

The 65-point map is `phase_crossing_map.csv`. For every successful crossing sample, the full-period derivative scan found 18 separated tangent minima, with A/B the two lowest observed minima. Failed continuation rows explicitly mean failure to establish a crossing, not proof of its absence. At failed rows, λ and minimum count are left blank rather than assigned fabricated values.

Independent 65-digit roots plus directed Krawczyk checks at φ/π=1,0.75,0.5,0.4 verify positive A/B curvatures and positive equal-cost slopes. The corresponding λc values are approximately 0.06598041113, 0.07207353166, 0.06928102462, 0.06216540978. They certify local crossings at those phase points, not full-period globality for each phase. Each regular point implies an open phase neighborhood, so exchange is not a measure-zero consequence of the π symmetry. The strict global tangent margins at π further imply an existential open global phase neighborhood of π, without a numerical width.

The tracked equal-cost curve numerically ends at a cusp:

    v*=5.11554281869418, λ*=0.0581287218835350,
    φ*/π=0.339281485775696,
    Rv=Rvv=Rvvv=0, Rvvvv≈0.00101860386027>0.

The two minima and intervening maximum merge; the two-branch cost slope also tends to zero as the minima become identical. This is a cusp/annihilation mechanism, not merely a solver's first failed step. The cusp solve is numerical, not an interval proof of global nonexistence below φ*. The branch topology suggests exchange for φ/π in (0.3392815,1], with one exceptional central collision at φ/π≈0.46495168622: v_A=0, v_B≈8.33021724, λ≈0.06681386255. At that point the separated satellite chart fails, although the deflated scalar cost continues. Thus the numerically observed **noncoalescent** phase family excludes this isolated value. For the connected family containing π, the candidate phase interval is (0.46495168622π,π], not the entire cusp-to-π interval.

The absence of tracked crossings at φ=0,π/4 in the former finite-z search is consistent with this tangent mechanism, but the cusp does not by itself prove a finite-z no-crossing theorem.

Run `python stage11_phase_resolution.py phase` and `python stage11_phase_validation.py`. The latter records all interval enclosures and the high-precision cusp/central-collision solves in `stage11_phase_validation.json`.
''')
    rs=list(csv.DictReader(open(ROOT/'resolution_sweep.csv')))
    tr=[dict(z=f"{float(r['z']):.7g}",bins=f"{float(r['separation_DFT_bins']):.5f}",epsilon=f"{float(r['epsilon_cross']):.10f}",third_margin=f"{float(r['other_minimum_margin']):.7g}") for r in rs]
    write('resolution_regime_analysis.md',r'''
# Resolution-regime audit

The true strong-pair separation is 2z/(2π)=z/π DFT bins. At z=2 it is 0.63661977 bins. “Sub-bin” describes that spacing; it is not an estimator-independent statistical resolution verdict.

The requested seven points all have separated positive-curvature A/B equal-cost roots:

'''+table(tr,['z','bins','epsilon','third_margin'])+r'''

Each requested point used a 512×512 full-frequency-torus search, local refinement from up to 200 detected grid minima and 32 random starts, explicit refinement initialized at the strong generating pair, and multiprecision stationary checks. This is high-accuracy numerical global ordering; only z=2 has a full-domain interval certificate. At z=5 and the decimal approximation 6.283185307179586, local A/B stationary roots and positive Hessians were additionally checked with directed intervals. The latter is a point extremely close to 2π, not an interval certificate for exact 2π.

The initial direct jump from z=4 to z=5 failed to converge. Focused continuation in increments 0.025 recovered the branch, whose B curvature becomes small near z≈4.7 and then increases again. Treating the first solver failure as a resolution boundary would have been wrong. `stage11_resolution_boundary.json` records this diagnosis and the local validations.

In every tested crossing case, refinement initialized at (−z,z) lands on A. A fixed-frequency strong-pair fit has larger cost. There is no separate winning “direct strong-pair” branch among these tested crossings. Below the crossing A remains the candidate preferred representation; above it B does. No z_res was identified: exchange survives numerically to two bins, and no universal Rayleigh limit is claimed.

Two scope qualifications matter. At z=5 the crossing amplitude is 1.05236, so the designated omitted component is no longer weaker than either unit-amplitude generating tone. At z≈2π its amplitude is 0.98343, nearly equal. Also, keeping b=10 fixed makes the positive strong tone move closer to the omitted tone as z grows. The experiment therefore follows this particular three-tone family; it does not isolate strong-pair resolution independently of all three pairwise spacings.

No connected global continuation or above-one-bin full-domain proof is inferred from these searches. `resolution_sweep.csv` preserves frequencies, losses, Hessian eigenvalues, direct-fit controls and the observed third-minimum margins.
''')
    mc=f['minimum_outer_cell']
    write('finite_epsilon_certificate.md',r'''
# Finite-amplitude certificate audit and attempted extension

## Historical status

The N=21 cases at εc−0.035, εc, εc+0.035 are **full interval-certified point cases**, not merely high-precision illustrations. The outer logs contain 12,288/14,287/10,961 excluded cells and 76/90/51 accepted cells. The respective inner searches complete with zero unresolved cells. The original clean replay regenerated those records. Figure 1 may call these three amplitudes separately certified; it may not claim the whole interval between them.

## Finite-width attempt

The requested [εc−0.03,εc+0.03] was split into 12 slabs of width 0.005. Each slab used a multiprecision root proposal at its midpoint, parameter-interval gradients/Hessians and Krawczyk boxes of radii 0.005,0.01,0.025,0.05,0.1 in normalized frequency. A uniform incumbent came from a fixed feasible frequency pair evaluated with interval ε. Forty weakest historical outer cells were tested in each slab using interval ε.

None of these slab attempts closed: no tested root-box radius passed both inclusion and positive-Hessian checks, and none of the 40 weakest reused spatial cells passed in any slab. These failures are recorded rather than silently discarded. The original natural interval expression overestimates frequency/amplitude dependency; shrinking a scalar bracket alone is not a proof of a finite-width regime. No complete new outer/inner partition was generated after this failed prerequisite diagnostic. The finite-width result is **CERTIFICATION-LIMITED**, not certified and not mathematically refuted. A predictor-centered parameterized Krawczyk/Taylor model and a new spatial partition would be needed to continue this attempt.

## Location of the minimum historical margin

The attaining physical-frequency cell is

    m∈[0.0337475773334841,0.03681553890925539]
    h∈[0.07363107781851078,0.07669903939428206].

Its center is u≈(−0.83755351019,2.31937895128), and its coordinate enclosure is

    u1∈[−0.90198070328,−0.77312631709],
    u2∈[2.25495175819,2.38380614437].

The historical lower bound is 0.8118137698122954, only about 5.8978946e−6 above the incumbent. The actual center cost is 1.27373900320810, about 0.46193 higher. Its center is at L∞ distance 3.93489 from A and 10.14404 from B, or 1.93489 and 8.14404 outside the corresponding radius-two boxes. A multiprecision unconstrained stationary solve from this center converges to A outside the cell. More decisively, a 31-node directed interval-gradient subdivision excludes stationary points throughout the cell. It is **enclosure slack, not a third basin**.

All attempted slab boxes, failed cells, bounds, midpoint roots, and the gradient-exclusion result are in `finite_epsilon_certificate.json`. The pre-existing narrow crossing bracket remains valid; no finite-width claim has been substituted for it.
''')
    write('projector_bound_derivation.md',r'''
# Explicit projector displacement lemma used by the certificates

Let t be the centered integer sample index, L=max|t|=(N−1)/2, T2=(sum t⁴)^1/2, T3=(sum t⁶)^1/2. Frequencies m,h are in physical radians/sample. Set w1(h)=cos(ht), w2(h)=t sinc(ht), and qj(m,h)=exp(imt)wj(h)/||wj(h)||. The two columns are orthonormal: w1 is even and w2 odd. Write Q=[q1,q2], P=QQ*.

For a rectangle of radii rm,rh around (mc,hc), evaluate

    cj=||wj(hc)||,
    a1=||−t sin(hc t)||,
    a2=||t² sinc′(hc t)||,
    δ1=a1 rh + (T2/2) rh²,
    δ2=a2 rh + (T3/6) rh².

If cj>δj for both j, the exact implemented bound is

    θ = sqrt(2) L rm
         + sqrt[(δ1/(c1−δ1))²+(δ2/(c2−δ2))²].

All quantities in the exclusion implementation are outward interval enclosures, and the radii include an additional 10^-12 radians. Replacing θ by min(1,θ) would be valid but is not the archived formula.

## Derivation

1. sinc(x)=∫_0^1 cos(ax) da gives |sinc″(x)|≤1/3 for all real x. Thus ||w1″||≤T2 and ||w2″||≤T3/3. Taylor's integral remainder yields ||wj(h)−wj(hc)||≤δj.
2. Normalize along the straight segment wj(hc)+τ[wj(h)−wj(hc)]. Its norm is at least cj−δj. Differentiating normalization projects the increment onto the orthogonal complement of the normalized vector, so its derivative norm is at most δj/(cj−δj). Integration from τ=0 to1 bounds ||qj(mc,h)−qj(mc,hc)|| by that ratio. Taking the Frobenius norm of the two-column difference gives the square-root term above. This normalization argument does not introduce a factor of two.
3. At fixed h, ∂m qj=i diag(t)qj, hence ||∂m Q||F≤sqrt(2)L. Integration over the m displacement gives the first term. For equal-rank orthogonal projectors, ||P−Pc||2=||(I−Pc)Q||2≤||Q−Qc||2≤||Q−Qc||F. A triangle decomposition through (mc,h) proves the stated θ.
4. The residual reverse triangle inequality gives ||(I−P)x||≥||(I−Pc)x||−||P−Pc|| ||x||, and hence the implemented cell lower bound [max(0,sqrt(Jc_lower)−||x||upper θupper)]². For uncertain ε, the same formula is uniform when every operand encloses the full ε interval; dependency inflation can make it uninformative without making it invalid.

## Source correspondence

`three_to_two_tone_stage2_outer_global.py:verify_cell` uses N=21,L=10. `stage7_n31_certificate.py:directed_cell` computes L,T2,T3 from N=31. `D1,D2` in those routines are δ1,δ2 above, not Dirichlet derivatives. `center_projected_J`/`center_projected` evaluates Jc from the two orthogonal basis columns. The float proposal in `three_to_two_tone_stage1_global.py` uses the same formula but is never a validity predicate.

## Every removable-singularity case

- D_N and every derivative used by interval Gram/Taylor/root calculations are finite cosine/sine sums, including zero and periodic aliases. The displayed sine quotient is not used across a zero denominator.
- The archived outer routines evaluate sinc(z)=sin(z)/z and sinc′(z)=(z cos(z)−sin(z))/z² **only at strictly positive center h and positive t**. t=0 is handled explicitly by the initial constant terms. Every recorded h center is positive, even for a cell whose lower edge is h=0. Consequently those point intervals do not contain zero.
- Coverage of h=0 is provided by the analytic sinc extension and the uniform derivative bound above; no quotient is evaluated on an interval crossing h=0. The original helper functions would not safely evaluate such an interval and should not be advertised as doing so.
- For a new center exactly at zero or a direct interval sinc implementation, use the integral/entire-series definition: sinc(0)=1, sinc′(0)=0; sum (-1)^k z^(2k)/(2k+1)! with a factorial remainder (and differentiated remainder). Do not use quotient division across zero.
- Inner Gram-quotient boxes lie in separated A/B neighborhoods with positive determinant. Near-coalescent outer cells use the centered basis instead. The new tangent calculation additionally deflates A/v⁴ and B/v² at v=0, as documented in its report.
''')
    write('inner_taylor_bound_derivation.md',r'''
# Exact inner interval Taylor and exclusion predicates

Let X=X1×X2 be a normalized-frequency box, u0 its midpoint, Y=X−u0, and E the amplitude point/bracket. Evaluate the interval jet J0=J(u0;E), g0=∇uJ(u0;E), H=∇u²J(X;E). The exact implemented scalar enclosure is

    L = J0 + sum_i g0_i Y_i + (1/2) sum_{i,k} H_ik Y_i Y_k.

It follows by the second-order Taylor integral formula along u0+t(u−u0); the averaged Hessian lies in H. Repeated factors are ordinary interval products, including Y_i Y_i, not dependency-aware squared intervals. This may lower the bound unnecessarily but remains safe. Exclude a cell by objective only if lower(L)>upper(incumbent).

The centered gradient enclosure is

    G_i = g0_i + sum_k H_ik Y_k.

If either G_i excludes zero, there is no stationary point in X. Otherwise choose a floating inverse midpoint Hessian C, convert each entry to a decimal interval point, and compute

    K = u0−C g0 +(I−C H)Y.

If any coordinate of K is disjoint from X, there is no stationary root. A numerical C is a permitted preconditioner: its accuracy is not a premise; the interval inclusion/disjointness inequality is the premise. Root existence/uniqueness requires the separate inclusion K⊂int X and contraction test, not merely failure of exclusion. Positive H11 and det H certify a minimum on a root box.

This matches `three_to_two_tone_stage2_inner_global.py:process` term for term. An inner cell entirely inside a separately validated large root box is accepted by containment; all other cells must pass one of the above exclusions or be subdivided. Objective exclusion alone may retain a nonstationary low-bound region, so the gradient and Krawczyk tests are essential. Global completeness also requires the outer torus cover and periodic boundary identification, not just these local predicates.

The 31-node check of the minimum-margin outer cell reused the centered gradient formula; no stationary point remains inside that cell. For finite-width ε, all jets must carry the entire parameter interval. The attempted width-0.005 slabs failed their local prerequisite; they are not certified by the historical point calculations.
''')
    write('misspecification_resolution_related_work.md',r'''
# Misspecification and resolution: verified additions

Checked 26 September 2026. These are focused positioning additions, not a new exhaustive novelty claim.

| Reference | Established result and relevance | Distinction of this paper |
| --- | --- | --- |
| H. White, “Maximum Likelihood Estimation of Misspecified Models,” Econometrica 50(1), 1–25, 1982. DOI 10.2307/1912526. [Publisher record](https://www.jstor.org/stable/1912526) | Quasi-likelihood targets and asymptotic inference under a wrong model. | We do not introduce pseudo-true parameters or claim local covariance theory is generally invalid. We examine coexistence of two finite-record best lower-order approximations and the crossing of their costs. |
| S. Fortunati, F. Gini, M. S. Greco, C. D. Richmond, “Performance Bounds for Parameter Estimation under Misspecified Models: Fundamental findings and applications,” IEEE SPM 34(6), 142–157, 2017. DOI 10.1109/MSP.2017.2738017. [Author preprint](https://arxiv.org/abs/1709.08210), [author institutional record](https://scholars.duke.edu/publication/1645361) | Misspecified bounds, pseudo-true targets and robustness to model error are established SP theory. | The two-mode noise observation is supplementary evidence about selection between competing pseudo-true centers, not a replacement for those bounds. |
| S. T. Smith, “Statistical Resolution Limits and the Complexified Cramér–Rao Bound,” IEEE TSP 53(5), 1597–1609, 2005. DOI 10.1109/TSP.2005.845426. [Paper](https://www.researchgate.net/profile/Steven-Smith-56/publication/3319243_Statistical_Resolution_Limits_and_the_Complexified_Cramer-Rao_Bound/links/57431b8008ae9f741b3875d2/Statistical-Resolution-Limits-and-the-Complexified-Cramer-Rao-Bound.pdf) | Resolution based on comparing source separation with uncertainty in estimating that separation, including unknown complex amplitudes. | A DFT-bin comparison is geometric and does not establish this statistical criterion. Our noiseless question holds fitted order at two for a three-tone record. |
| D. Batenkov, G. Goldman, Y. Yomdin, “Super-resolution of near-colliding point sources,” Information and Inference 10(2), 515–572, 2021. [Publisher](https://academic.oup.com/imaiai/article/10/2/515/5820889), [author preprint](https://arxiv.org/abs/1904.09186) | Minimax recovery rates for clustered nodes and amplitudes depend on separation and bandwidth. | Recovery of the true cluster under perturbations differs from selecting the globally best deliberately lower-order representation. Neither result is superseded by the branch-exchange example. |

The introduction now positions misspecification first, and distinguishes statistical/source-recovery resolution from representation switching. No universal sub-Rayleigh threshold, failure of CRLB theory, or new general misspecification principle is claimed.
''')
    # Audit actual N31 dependency closure and its archived copies; keep source immutable.
    entries=[
      ('three_to_two_tone_stage3_full_crossings.py','module ref / baseline imports','stage3_tangent_reference.json; stage2_reference.json','harmless plotting/metadata','baseline unused; ref used for defaults and comparison fields, separated below'),
      ('three_to_two_tone_stage3_full_crossings.py','seed is None branch','N21 v_A,v_B,kappa,lambda','initialization only','N31 supplies explicit N31 seed; this branch is bypassed'),
      ('three_to_two_tone_stage3_full_crossings.py','lam,lam4 after stationary solve','N21 lambda_N,lambda_quartic','harmless plotting/metadata','incorrect N31 metadata; no feedback to root solve or predicate'),
      ('three_to_two_tone_stage2_inner_global.py','module ref','N21 stage2_reference.json','initialization only','initial known_center assigned then overwritten before in_tiny by known_center_override'),
      ('three_to_two_tone_stage2_inner_global.py','module local','N21 stage2_local_interval_checks.json','harmless plotting/metadata','loaded dictionary is unused by process'),
      ('three_to_two_tone_stage2_inner_global.py','box radii','2.00000001; 0.000299999999; 0.009999999','proof predicate input','shared tuning constants, not transferred proof results; N31 large boxes and accepted-cell containment independently verified'),
      ('three_to_two_tone_stage2_interval.py','module default N=21','21','proof predicate input','explicit interval.N=31 before all N31 validity calls; finite sums use current N'),
      ('stage7_n31_certificate.py','seed,roots,epsilon','N31-specific explicit seed and freshly solved crossing','initialization only','all roots and crossing signs independently interval verified'),
      ('stage7_n31_certificate.py','local_check inverse Hessian','new C from N31 midpoint Hessian','proof predicate input','no cached N21 preconditioner; converted to interval points and checked anew'),
      ('stage7_n31_certificate.py','adaptive(...,N) and directed_cell','full-torus spatial partition and L,T2,T3','proof predicate input','fresh N31 partition; no N21 cells read; all exclusion and containment checks repeated'),
      ('stage7_n31_inner_log.py','N31 certificate roots and epsilon','N31_global_certificate.json','proof predicate input','reads N31 u_A/u_B and bracket only; both known_center and radius overrides supplied'),
      ('stage10_replay_runner.py','three copied N21 JSON files','import-time dependencies','harmless plotting/metadata','not certificates accepted by the N31 verifier; dependencies documented, tested by poisoning'),
      ('crossN_scaling.csv','N=31 rows','lambda31,c31 and N31 tangent roots','harmless plotting/metadata','independently computed numerical N31 coefficients; correct comparison, no local N31 theorem'),
    ]
    stale=['epsilon_leading','epsilon_through_quartic','remainder','remainder_over_d4','relative_error','remainder_after_quartic','remainder_after_quartic_over_d6','relative_error_quartic']
    for k in stale:entries.append(('N31_global_certificate.json','high_precision_numerical_root.'+k,'N21 coefficients carried into N31 metadata','harmless plotting/metadata','scientifically mislabeled and superseded; not used in proof predicates'))
    audit=[dict(file=a,location=b,quantity=c,classification=d,data_flow=e) for a,b,c,d,e in entries]
    # List each replay duplicate of an audited source/file, not just its pattern.
    for a,b,c,d,e in entries:
        copy=ROOT/'stage10_replay'/'n31'/a
        if copy.is_file():audit.append(dict(file=str(copy.relative_to(ROOT)),location=b,quantity=c,classification=d,data_flow='archived replay copy: '+e))
    with (ROOT/'N31_leakage_audit.csv').open('w',newline='') as ff:
        w=csv.DictWriter(ff,fieldnames=audit[0]);w.writeheader();w.writerows(audit)
    poison=read('stage11_N31_poison_test.json')
    write('N31_leakage_audit.md',r'''
# N=31 leakage audit

## Finding

**No N=21 certified result was found to enter an N=31 validity predicate.** There is genuine N=21-to-N=31 comparison-metadata leakage, broader than the single quartic field. The shared crossing routine unconditionally computes eight comparison fields from the N=21 tangent file after solving the actual N=31 equations. All eight are listed in `N31_leakage_audit.csv` and must not be read as N31 tangent estimates.

The stale metadata field was generated by carrying the N=21 quartic coefficients into the preliminary N=31 comparison metadata.

## Dependency audit

The audit covers the two N31 runners, their shared crossing/interval/inner/proposal dependencies, imported JSON defaults, generated N31 JSON/CSVs, cross-N comparison table, replay runner and N31 replay copies. The CSV explicitly classifies root defaults, coefficient fields, bracket origins, partition widths, preconditioners and copied files. Shared box widths are algorithmic tuning choices that enter containment predicates; their validity is re-established with N31 roots, Hessians, Gram determinants and containment tests. They are not imported N21 conclusions. No N21 cached partition or inverse Hessian is accepted as a proof premise.

Two stale defaults do exist at import time. The inner routine first assigns an N21 root and then overwrites it with the explicit N31 override before testing any cell. The interval objective module starts with N=21 but the N31 runner assigns interval.N=31 before evaluation. This design is fragile and should be refactored in a future clean implementation; the frozen proof path nevertheless supplies the overrides explicitly.

## Adversarial data-flow check

`stage11_leakage_audit.py` replaced N21 λ,c with 123,456 and replaced tangent root guesses with unrelated values. With the explicit N31 seed, the freshly recomputed crossing, four root coordinates and both objective values were bit-for-bit unchanged as decimal strings; only the eight comparison fields changed. It also poisoned both N21 default inner root centers with values near 1000, reran the N31 inner exclusion using the documented overrides, and compared every serialized terminal row with the archive.

'''+f"- Root-field invariance: {poison['root_fields_exactly_equal']}.\n- Poisoned-default inner rows identical: {poison['inner_rows_identical']} ({poison['inner_row_count']} rows).\n- N31 frozen SHA-256: `{poison['frozen_sha256'].upper()}`.\n\n"+r'''
This supplements the existing clean full replay and selected independent MPFR checks. It is a targeted dependency test, not a new independent backend proof of every outer cell. No frozen certificate or generating script was edited.

## Correct comparison

'''+f"For N21, 4λ21+16c21={pred}; εc={base['epsilon_cross'][:36]}. Absolute error={err}; relative error={err/Decimal(base['epsilon_cross'])} ({100*err/Decimal(base['epsilon_cross']):.9f}%).\n\n"+r'''
The corrected numerical N31 prediction remains 0.24857675463874848, with absolute discrepancy 0.00056633521375329 and relative discrepancy 0.227313%. It is secondary: only N21 has the proved local quartic law.
''')
    write('ARCHIVE_NOTE_N31_CERTIFICATE_REVISED.md',r'''
# Revised provenance note — frozen N31 certificate

This note supersedes the characterization in `ARCHIVE_NOTE_N31_CERTIFICATE.md`. The frozen JSON is retained unchanged, SHA-256 E4E32DC3595BCAF18EE48DF8C95152187C3944D6A61E6980B6EFCC3C87D57CAF.

**The stale metadata field was generated by carrying the N=21 quartic coefficients into the preliminary N=31 comparison metadata.** Its value 0.24763765574844870339879795294973264925… is exactly 4λ21+16c21. Calling it merely a preliminary extrapolation omitted the source of the mismatch.

The other seven fields derived from the same N21 coefficients—epsilon_leading, remainder, remainder_over_d4, relative_error, remainder_after_quartic, remainder_after_quartic_over_d6, relative_error_quartic—are likewise superseded for N31 comparison purposes. The actual N31 crossing, fitted roots, costs, root boxes, Hessians, Gram determinants, bracket signs, slopes and cell predicates are independently evaluated at N31 and are unaffected.

The corrected N31 numerical comparison is 4λ31+16c31=0.24857675463874848; it is not an N31 local theorem. `N31_leakage_audit.md`, its CSV, and the poisoned-input test document the data flow and independent predicate checks. No silent alteration of the historical evidence file has been made.
''')
    write('major_revision_summary.md',r'''
# Major scientific revision

## Verdict: TSP-STRENGTHENED

This is a scientific strengthening of the narrow family, not a prediction of editorial acceptance.

1. **Small-z globality succeeds.** Directed interval isolation exhausts the tangent frequency period (18 minima,18 maxima); A/B beat all other separated minima. A recurrence-based classification also covers every central/coalescent scaling and diverging-amplitude limit. The best such competitor has coefficient 0.07681023 versus A/B 0.05876262. Compactness then gives an existential small-z global exchange theorem. No explicit radius or certified bridge to z=2 is claimed.
2. **Phase is mechanistically explained.** The exact arbitrary-phase tangent loss has in-phase and quadrature projection terms. Representative interval checks verify crossings at φ/π=0.4,0.5,0.75,1. Numerical continuation identifies a cusp near 0.33928149π and a central-chart collision near 0.46495169π. The complete phase interval is not globally interval-certified.
3. **Exchange is not confined to sub-bin spacing.** All seven requested points have numerical A/B global ordering, including z=π and near2π. Local roots at z=5 and near2π are interval checked. No disappearance boundary was found. At z=5 the omitted amplitude exceeds one; this qualification is explicit.
4. **The actual projector and Taylor predicates are now explicit.** Reports and supplement give the implemented θ, constants, norm-normalization argument, singularity handling and inner exclusion formulas.
5. **N31 proof dependencies pass the leakage audit.** Eight comparison metadata fields use N21 coefficients and are explicitly superseded. Poisoned N21 coefficients leave N31 roots/costs unchanged; poisoned N21 inner defaults reproduce the archived N31 terminal cells. The independently evaluated proof predicates remain valid.
6. **Finite-width epsilon certification does not succeed.** Historical ±0.035 cases are certified isolated points. Twelve width-0.005 parameter slabs fail the tested local/outer enclosure predicates. The smallest historical outer margin comes from a cell with no stationary point, not a third basin. The finite-width interval remains certification-limited.
7. **Noise is reframed around pseudo-true branches.** The main error summary is total variance = within-branch variance + between-branch selection variance. MSE to the generating strong pair is removed. No new Monte Carlo study was added.

## Evidence and files

The requested analytical reports, CSV maps, machine-readable certificates, revised manuscript and supplement are in this workspace. New executable scripts have the `stage11_` prefix. `smallz_global_certificate.json` separates analytical proof dependencies from computed interval predicates. The old N21/N31 finite-z certificate JSONs remain frozen.

## Remaining limitations

No numerical small-z radius; no connected global bridge to z=2; no full-domain interval proof of the large-z numerical cases; no certified full phase boundary; no finite-width amplitude certificate. The result remains N21,b10 with equal strong amplitudes, plus the independently certified finite N31 example. These are scope limitations, not evidence of a losing tangent competitor or contaminated proof chain.
''')
    write('response_to_hostile_review.md',r'''
# Response to the substantive review

| Objection | Scientific revision | Remaining limit |
| --- | --- | --- |
| Two local minima do not prove a global small-z mechanism. | Exhaustive interval tangent isolation; exact strong-pair control; recurrence proof covering all coalescent limits; strict margins and compactness global theorem. | Radius remains existential. |
| The direct fit to the strong tones is omitted. | Optimized-amplitude coefficient 0.08550039 and fully reoptimized coalescent coefficient 0.07681023 both exceed A/B0.05876262. | Specific to the stated family and scaled λ neighborhood. |
| Phase dependence is unexplained. | General-phase projection formula, 65-point map, four interval-verified crossings, numerical cusp and central collision. | Cusp location and full phase existence boundary are numerical. |
| The example may only reflect sub-Rayleigh spacing. | Explicit bin conversion; requested sweep through two bins; exchange survives all seven cases. | At z5 the omitted amplitude is slightly larger than1; fixed b changes other pair spacings. No universal resolution law. |
| Global epsilon certificate is infinitesimal. | Audited all historical cases; attempted 12 finite-width slabs and preserved failed predicates. | Wider interval remains uncertified. This objection is not claimed resolved. |
| Tiny exclusion margin could hide a competitor. | Exact attaining cell reported; center loss is much higher; interval-gradient subdivision excludes any stationary point. | Bound looseness remains a computational issue. |
| θ cannot be reconstructed. | Explicit normalized-basis lemma, δ1/δ2 constants and exact implemented θ; separate inner Taylor predicate. | Existing code only evaluates sinc quotients at nonzero centers, now explicitly documented. |
| N31 may be contaminated by N21. | Eight inherited comparison fields identified; full dependency table and poisoned-input checks; corrected archive wording. | Frozen metadata retained for provenance, clearly superseded. |
| Noise MSE uses an inappropriate target. | Focus shifted to branch probabilities, branch-conditioned covariance, and pseudo-true center separation. | Frozen-projector law remains an approximation. |
| Missing misspecification/resolution literature and notation inconsistencies. | White, Fortunati, Smith and clustered-superresolution references; Δ,φ,ω,ν,u and λ21 usage fixed. | No general estimator or universal resolution claim. |

The word “Critical” is removed from the title. The abstract names interval-verified nondegeneracy and tangent classification. The main theorem now states only the global conclusion justified by the complete limiting-chart argument. The paper does not present failed finite-width validation as a success.
''')
if __name__=='__main__':main()
