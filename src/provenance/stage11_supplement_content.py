from pathlib import Path
R=Path(__file__).resolve().parent/'supplement'
texts={
'S1_complete_local_proof.tex':r'''
\section*{S1. Small-spacing global proof and quartic law}\label{app:local}
Fix $N=21,b=10,\phi=\pi$, centered integer $t$, $s=t/N$, $p=z^2$ and $\epsilon=\lambda p$. All sums are finite; no continuum or orthogonality approximation is used. The conclusion has an existential positive radius.

\subsection{Exact projection and tangent elimination}
Write $D(q)=\sum_t e^{iqs_t}$ and $S_k=\sum_t s_t^k$. Then $D^{(k)}(q)=\sum_t(is_t)^ke^{iqs_t}$, with $S_2=(N^2-1)/(12N)$ and $S_4=(N^2-1)(3N^2-7)/(240N^3)$. The sixth moment $S_6=(N^2-1)(3N^4-18N^2+31)/(1344N^5)$ enters the products in the next-order residual coefficient and the projector derivative norm; it is not an independent hypothesis. For noncoincident fitted normalized frequencies,
\begin{equation}\label{eq:app-gram2}
G=\begin{pmatrix}N&D(u_2-u_1)\\D(u_2-u_1)&N\end{pmatrix},\qquad
J=\|x\|^2-h^*G^{-1}h,\quad h_j=D(u_j-z)+D(u_j+z)-\epsilon D(u_j-b).
\end{equation}
The desingularized fit is $y=(2+\alpha p)e^{i\kappa ps}+\beta pe^{ivs}$, with $\alpha,\beta$ complex and $\kappa,v$ real. The divided residual extends analytically to $p=0$ with
\begin{equation}\label{eq:app-r0}
r_0=-s^2-\lambda e^{ibs}-\alpha-2i\kappa s-\beta e^{ivs}.
\end{equation}
For the real inner product $\Re\langle\cdot,\cdot\rangle$, the five columns $[1,e^{ivs},2is,i,ie^{ivs}]$ have block Gram matrices
\begin{equation}
M_1=\begin{pmatrix}N&D(v)&0\\D(v)&N&-2D'(v)\\0&-2D'(v)&4S_2\end{pmatrix},\qquad
M_2=\begin{pmatrix}N&D(v)\\D(v)&N\end{pmatrix}.
\end{equation}
At phase $\pi$ the second block has zero right-hand side. Define $q_0=-S_2-\lambda D(b)$, $q_1=2\lambda D'(b)$, $q_3=D''(v)-\lambda D(b-v)$,
\begin{align}
A&=N-D(v)^2/N-D'(v)^2/S_2,&B&=q_3-D(v)q_0/N+D'(v)q_1/(2S_2),\\
R&=S_4-2\lambda D''(b)+N\lambda^2-q_0^2/N-q_1^2/(4S_2)-B^2/A,\nonumber\\
\beta&=B/A,&\alpha&=(q_0-D(v)\beta)/N,\quad\kappa=(q_1+2D'(v)\beta)/(4S_2).\nonumber
\end{align}

\subsection{Regular local wells and the analytic crossing}
The three equations $R_v(v_A,\lambda)=R_v(v_B,\lambda)=R(v_A,\lambda)-R(v_B,\lambda)=0$ have an interval-verified root near
\begin{equation}
v_A=-4.15344048664877,\quad v_B=12.4680948443009,\quad\lambda_{21}=0.0659804111269914.
\end{equation}
The recorded lower bounds are $R_{vv,A}>0.0027512694$, $R_{vv,B}>0.0197033302$, $\partial_\lambda(R_A-R_B)>2.2286882219$, $A(v_A)>5.24629611$, $A(v_B)>19.30370897$, $\beta_A>0.1225281271$ and $\beta_B<-0.0638766035$. The original Krawczyk radii are $10^{-9},10^{-9},10^{-10}$ with contraction below $2.224\times10^{-6}$. The new full-period classification revalidates the crossing at radii $10^{-24},10^{-24},10^{-25}$.

For the full divided loss $L(p,\lambda,w)$, with six real fit coordinates $w$, the linear-coordinate Hessian is twice the positive Gram matrix and its Schur complement is $R_{vv}>0$. Hence $L_{ww}>0$. The analytic implicit-function theorem produces unique local minimum branches $w_j(p,\lambda)$; for $p>0$ their ordinary coordinates are nonsingular since the satellite amplitude is $p\beta_j+O(p^2)\ne0$ and satellite frequency tends to $v_j\ne0$. The central frequency is $\kappa_jp+O(p^2)$. After shrinking the neighborhood, the physical within-fit separation is at least $|v_j|/(2N)$. A second implicit-function theorem applied to $L_A-L_B$ yields the unique analytic crossing $\lambda_c(p)$, since its $\lambda$ derivative is positive.

\subsection{Complete tangent classification and the strong-pair control}
Directed interval isolation covers $[-66,66]$, containing $[-21\pi,21\pi]$. There are exactly 18 minima and 18 maxima in one period. A/B have common loss $0.05876262329535$; every other stationary minimum is higher by more than $0.06444970$. The terminal cover and all root/cost/curvature intervals are in \path{stage11_tangent_partition.csv}, \path{stage11_tangent_roots.csv}, and \path{smallz_global_certificate.json}; 32,406 cells were visited with none unresolved. Each terminal cell either excludes zero from $R_v$ or has fixed-sign $R_{vv}$ and opposite endpoint signs for $R_v$, proving exactly one root. The redundant end pieces contain no root.

For $|v|\le1$, factor $A=v^4\widetilde A$, $B=v^2\widetilde B$ before division. Degree-40 finite-moment Taylor polynomials enclose $\widetilde A,\widetilde B$ and their first two derivatives, with an added outward remainder $10^{-35}$. The exponential coefficient majorant $10^4/k!$ gives derivative tails below $10^4\sum_{k\ge41}k^2/(k+2)!<10^{-43}$, using $2\max|s|<1$. Other cells use finite trigonometric sums. Thus the central singularity is included rather than discarded. Independent five-column multiprecision least squares agrees with all 36 stationary cost enclosures and six deflated near-zero checks.

Fixing fitted frequencies at $\pm z$ and optimizing their amplitudes gives
\begin{equation}
C_{\rm strong}=\lim_{z\to0}J/z^4=\lambda^2\{N-D(b)^2/N-D'(b)^2/S_2\}=0.08550038682621.
\end{equation}
This follows because the strong pair is exactly in the fitted span and that span tends to complex $\operatorname{span}\{1,s\}$. Allowing both frequencies to change yields the limiting class $\alpha+\beta s+\gamma s^2$, $\alpha,\beta\in\mathbb C$, $\gamma\in\mathbb R$, whose best loss is
\begin{equation}\label{eq:app-coal}
C_{\rm coal}=C_{\rm strong}-\lambda^2\frac{\{D''(b)+S_2D(b)/N\}^2}{S_4-S_2^2/N}=0.07681023119888.
\end{equation}
The optimum has $\gamma=-0.725519235\ldots$ and is attainable with frequencies $\pm\sqrt{-\gamma}\,z$ and amplitude corrections. Its margin above A/B exceeds $0.01804760$. A satellite fixed at $b$ is already covered, with tangent loss $0.10362265533211$.

\subsection{Why these charts exhaust all global competitors}
Take any sequence $p\downarrow0$ and fits with residual $O(p)$; A/B trials ensure global minimizers satisfy this bound. Then $y=2+ph$ with bounded $h$. Extract limits of $h$ and nodes $\xi_j=e^{iu_j/N}$. Every fit satisfies
\begin{equation}\label{eq:app-recurrence}
y_{t+2}-(\xi_1+\xi_2)y_{t+1}+\xi_1\xi_2y_t=0.
\end{equation}
The limiting constant record forces $(1-\xi_1)(1-\xi_2)=0$: at least one node is central. If exactly one is central, the nonconfluent Gram is boundedly invertible and the five real derivative columns above are independent. A Taylor expansion with a bounded local left inverse forces central-frequency and coefficient perturbations to be $O(p)$; the limiting residual belongs to the enumerated satellite tangent class.

If both nodes are central, substitute $y=2+ph$ into \eqref{eq:app-recurrence}. Bounded $h$ forces $(1-\xi_1)(1-\xi_2)=O(p)$, hence $u_1u_2/p$ is bounded for local representatives tending to zero. Extract its real limit $\gamma$. The recurrence gives $\Delta_t^2h=2\gamma/N^2$, so $h=\alpha+\beta s+\gamma s^2$. This reasoning allows arbitrary relative rates, vanishing or diverging amplitudes, and exact confluent limits. Both nodes approaching a noncentral frequency are impossible at $O(p)$ residual.

Compactify the fitted frequency subspaces by adjoining confluent subspaces; their projectors are continuous and the compactified loss attains a minimum. A sequence of global minimizers escaping disjoint A/B neighborhoods would, by the preceding exhaustion, have either coalescent leading loss at least \eqref{eq:app-coal}, or a satellite tangent limit outside both wells. The former contradicts the strict margin; the latter contradicts the complete classification and compactness outside the wells. These exclusions persist for $\lambda$ in a sufficiently small neighborhood by continuity. Inside each well the positive desingularized Hessian gives the unique local branch. The closure minimizer is therefore an ordinary separated fit. This proves the global assertion with existential $p_0,\eta>0$. The displayed stationary margin is not a uniform gap to points arbitrarily close to A/B.

\subsection{Quartic coefficient and orientation}
At fixed fit coordinates the next divided-residual term is $r_1=s^4/12-i\kappa\alpha s+\kappa^2s^2$. Stationarity removes coordinate-derivative terms from the envelope, giving
\begin{equation}
J_j=p^2\{R_j+pT_j+O(p^2)\},\qquad T_j=2\Re\sum_t\overline{r_{0,j,t}}r_{1,j,t},\qquad
c_{21}=-\frac{T_A-T_B}{\partial_\lambda(R_A-R_B)}.
\end{equation}
Here $T_A=-0.00017117693705507491\ldots$, $T_B=-0.00243942282576626480\ldots$, and $c_{21}\in[-0.001017749297469898,-0.001017749297469697]$. Thus $\epsilon_c=\lambda_{21}z^2+c_{21}z^4+O(z^6)$. The fundamental theorem of calculus gives $J_A-J_B=p(\epsilon-\epsilon_c)\{\partial_\lambda(R_A-R_B)+O(p+|\lambda-\lambda_{21}|)\}$, whose last factor is positive locally. Analyticity supplies finite remainder constants on a smaller neighborhood, but no numerical radius or bound through $z=2$ is supplied.
''',
'S2_interval_certificate_protocol.tex':r'''
\section*{S2. Finite-spacing interval predicates}\label{app:global}
The N=21 and N=31 cases at $z=2,b=10,\phi=\pi$ use independent roots and full-torus partitions. Let fitted physical frequencies be $\nu_j$, normalized frequencies $u_j=N\nu_j$, mean $m$, and shortest-arc half-separation $h$. Up to periodic boundary identification, the unordered torus is covered by $m\in[-\pi,\pi]$, $h\in[0,\pi/2]$.

\subsection{Explicit centered projector displacement}
Use $c_{m,h}(t)=e^{imt}\cos(ht)$ and $d_{m,h}(t)=e^{imt}t\operatorname{sinc}(ht)$. Even/odd parity makes these columns orthogonal. Their norms are positive at $h=0$, where the span is $\{e^{imt},te^{imt}\}$. Put $w_1=\cos(ht)$, $w_2=t\operatorname{sinc}(ht)$, $q_j=e^{imt}w_j/\|w_j\|$, $Q=[q_1,q_2]$, $P=QQ^*$. For a rectangle centered at $(m_c,h_c)$ with radii $r_m,r_h$, define
\begin{align}
c_j&=\|w_j(h_c)\|,&a_1&=\|t\sin(h_ct)\|,&a_2&=\|t^2\operatorname{sinc}'(h_ct)\|,\\
L&=(N-1)/2,&T_2&=(\sum t^4)^{1/2},&T_3&=(\sum t^6)^{1/2},\nonumber\\
\delta_1&=a_1r_h+T_2r_h^2/2,&\delta_2&=a_2r_h+T_3r_h^2/6.&&\nonumber
\end{align}
If $c_j>\delta_j$, the exact implemented bound is
\begin{equation}\label{eq:app-theta}
\|P(m,h)-P(m_c,h_c)\|_2\le\theta
=\sqrt2 Lr_m+\left[\sum_{j=1}^2\left(\frac{\delta_j}{c_j-\delta_j}\right)^2\right]^{1/2}.
\end{equation}
Indeed $\operatorname{sinc}(x)=\int_0^1\cos(ax)\,da$ implies $|\operatorname{sinc}''(x)|\le1/3$. Hence $\|w_1''\|\le T_2$, $\|w_2''\|\le T_3/3$ and Taylor's integral remainder bounds each column displacement by $\delta_j$. Normalize along the straight segment between the unnormalized vectors. Its norm is at least $c_j-\delta_j$, while the derivative of normalization is the orthogonal component of the increment divided by that norm. Integration bounds the normalized displacement by $\delta_j/(c_j-\delta_j)$, giving the square-root term. For the mean-frequency displacement, $\partial_m q_j=i\operatorname{diag}(t)q_j$ and $\|\partial_mQ\|_F\le\sqrt2L$. Finally equal-rank projector geometry gives $\|P-P_c\|_2=\|(I-P_c)Q\|_2\le\|Q-Q_c\|_F$. Combining the two paths proves \eqref{eq:app-theta} without an extra factor of two.

The reverse triangle inequality gives the actual exclusion predicate, with directed endpoints,
\begin{equation}\label{eq:app-cell}
\left[\max\{0,\sqrt{\underline J_c}-\overline{\|x\|_2}\,\overline\theta\}\right]^2>U,
\end{equation}
where $U$ is a feasible incumbent upper bound. The archived N21 routine \path{three_to_two_tone_stage2_outer_global.py} uses $L=10$; \path{stage7_n31_certificate.py} computes all norms at N31. Both inflate radii by $10^{-12}$ physical radians. Floating point proposes the partition only. Every nonexcluded terminal cell must be wholly contained in an A/B inner neighborhood after inflation.

\subsection{Removable singularities and arithmetic}
All interval Dirichlet derivatives use finite sine/cosine sums; the sine quotient is never divided across a zero. The archived outer sinc and sinc-prime quotient helpers are evaluated only at positive cell-center $h_c$ and positive integer $t$. The $t=0$ terms are supplied separately. Even cells touching $h=0$ have strictly positive centers. Coverage of their edge follows from the entire sinc extension and the uniform derivative proof, not interval quotient evaluation across zero. These helpers would require finite Taylor sums and remainder bounds if called on a zero-containing argument; at zero, sinc is one and its derivative zero. Inner Gram quotients are used only in separated neighborhoods; outer confluent cells use the centered basis. The new tangent chart has the additional deflation described in S1.

\subsection{Root boxes and inner Taylor predicates}
For each crossing endpoint and the whole narrow bracket, interval Krawczyk inclusion/contraction and positive Hessian principal minors isolate A/B strict minima, while a positive Gram determinant verifies finite amplitudes. Opposite endpoint signs and a positive interval envelope slope imply a unique crossing in that bracket. Table~\ref{tab:app-cert} records the bounds.

For an inner box $X$ in normalized $u$ coordinates, midpoint $u_0$, $Y=X-u_0$, and amplitude interval $E$, evaluate $J_0=J(u_0;E)$, $g_0=\nabla J(u_0;E)$, and $H=\nabla^2J(X;E)$. The exact implemented Taylor enclosures are
\begin{equation}\label{eq:app-inner}
\mathcal L=J_0+\sum_i g_{0i}Y_i+\tfrac12\sum_{i,k}H_{ik}Y_iY_k,\qquad
\mathcal G_i=g_{0i}+\sum_kH_{ik}Y_k.
\end{equation}
The integral remainder places the path-averaged Hessian in $H$. Ordinary interval products $Y_iY_i$ are used, which may be loose but are safe. A box is excluded if $\underline{\mathcal L}>U$ or either $\mathcal G_i$ excludes zero. Otherwise, with a floating inverse midpoint Hessian converted to interval-point preconditioner $C$, form
\begin{equation}
K=u_0-Cg_0+(I-CH)Y.
\end{equation}
A disjoint coordinate of $K$ excludes roots. Strict inclusion and contraction, checked separately, establish root existence and uniqueness; positive $H_{11},\det H$ establish a minimum. An inner cell is accepted only when wholly inside a previously validated root box. Complete logs are \path{stage2_inner_cells.csv} and \path{N31_inner_interval_cells.csv}.

\subsection{Wider amplitude audit and minimum-margin cell}
The historical N21 values $\epsilon_c\pm0.035$ are separately interval-certified global points. Their outer exclusion counts are 12,288 and 10,961, with zero unresolved inner cells; the intervening amplitude interval is not thereby certified. The new twelve width-$0.005$ slabs spanning $\epsilon_c\pm0.03$ failed the tested parameter-interval root/outer predicates. Failed cells and tested radii are retained in \path{finite_epsilon_certificate.json}; no finite-width claim is made.

The smallest historical crossing margin, $5.8978946\times10^{-6}$, occurs at $m\in[0.03374757733,0.03681553891]$, $h\in[0.07363107782,0.07669903939]$, centered at $u=(-0.83755351,2.31937895)$. Its center cost is $1.27373900$, versus incumbent $0.81180787$; the small lower-bound margin is slack. A 31-node interval-gradient subdivision excludes every stationary point in this cell. It is not a third competitive basin.
''',
'S4_extended_generality.tex':r'''
\section*{S4. Phase, resolution, and extended numerical records}
For arbitrary phase the real Gram blocks in S1 remain unchanged. Put $c=\cos\phi$, $d=\sin\phi$, $q_0=-S_2+\lambda cD(b)$, $q_1=-2\lambda cD'(b)$, $q_3=D''(v)+\lambda cD(b-v)$. Use the corresponding $A,B$ from S1 and set $C=N-D(v)^2/N$, $H=\lambda d[D(b-v)-D(v)D(b)/N]$. Both Gram right-hand sides now contribute:
\begin{equation}
R=S_4+2\lambda cD''(b)+N\lambda^2-q_0^2/N-q_1^2/(4S_2)-\lambda^2d^2D(b)^2/N-B^2/A-H^2/C.
\end{equation}
The 65-point map \path{phase_crossing_map.csv} tracks equal-cost wells and scans the full tangent period at each successful crossing. Interval checks at $\phi/\pi=0.4,0.5,0.75,1$ validate positive curvatures and slope, proving local crossings on open phase neighborhoods. Their full boxes are in \path{stage11_phase_validation.json}. They do not certify globality for the whole phase interval. At $\pi$, the strict global tangent margins also persist on some open phase neighborhood by continuity.

Numerical solving of $R_v=R_{vv}=R_{vvv}=0$ gives $(v_*,\lambda_*,\phi_*/\pi)=(5.11554281869,0.0581287218835,0.339281485776)$ with $R_{vvvv}\approx0.00101860386>0$: the two wells and intervening maximum merge at a cusp. Separately $v_A=0$ occurs at $\phi/\pi=0.464951686220$, $\lambda=0.066813862548$, $v_B=8.33021724227$. The noncoalescent segment connected to $\pi$ therefore numerically ends at that chart collision; a second segment extends toward the cusp. Failed continuation below the cusp is not a global nonexistence proof. The exact in-phase $B$ cross term and quadrature $H$ term explain why phase reshapes the landscape.

At $z=2,2.5,\pi,3.5,4,5,2\pi$, \path{resolution_sweep.csv} records the numerical equal-cost points, direct strong-pair control, positive Hessians and independent full-torus ordering. The search used a $512\times512$ grid, up to 200 grid-minimum refinements and 32 random starts per point; the best candidate was checked by a multiprecision stationary solve. Interval checks additionally isolate the local A/B minima at $z=5$ and the decimal value $6.283185307179586$. Only $z=2$ has full-domain interval globality. Exchange survives through two DFT bins; no family-specific disappearance boundary was found. At $z=5$ the omitted amplitude exceeds one and the fixed location $b=10$ does not keep every pairwise spacing fixed.

The previous cross-length and one-factor checks remain in \path{crossN_scaling.csv}, \path{crossN_globality_checks.csv}, \path{asymmetry_sweep.csv}, and \path{existence_region.csv}. They include $N=11,15,21,31,41$, strong-amplitude ratios $0.8$--$1.2$, $b=8$--$12$, and strong-frequency shifts $-0.4$--$0.4$. These records supply numerical persistence checks, not a certified parameter rectangle or an all-$N$ theorem.
''',
'S7_archive_notes.tex':r'''
\section*{S7. N=31 metadata and proof-dependency audit}
The frozen N31 JSON retains SHA-256
\begin{center}\footnotesize\texttt{E4E32DC3595BCAF18EE48DF8C95152187C3944D6A61E6980B6EFCC3C87D57CAF}\end{center}
The stale metadata field was generated by carrying the N=21 quartic coefficients into the preliminary N=31 comparison metadata. Its value $0.2476376557484487\ldots$ is exactly $4\lambda_{21}+16c_{21}$. Seven associated leading/remainder/relative-error fields share the same provenance and are superseded for N31 comparisons. The actual N31 crossing, roots, objectives, brackets, Hessians, Gram bounds and cell predicates do not use those fields.

The dependency audit \path{N31_leakage_audit.csv} traces all shared defaults, radii, partition inputs, preconditioners and replay copies. N31 supplies its own seed and explicit inner-center overrides and sets the interval objective's N to31 before evaluation. Shared radius choices are independently verified by N31 inclusion, positivity and containment predicates; no N21 certificate result is a premise. Deliberately replacing the N21 coefficients by 123 and456 changes only the eight comparison fields, leaving recomputed N31 roots/costs identical. Poisoning N21 default inner centers reproduces all825 archived N31 terminal rows. These are targeted data-flow tests, supplementing the earlier full clean replay.

For the proved N21 law, the prediction at $z=2$ differs from its certified crossing by $0.00055297442383$ ($0.2228023\%$). The corrected numerical N31 prediction is $0.24857675463874848$ ($0.227313\%$ discrepancy); no N31 local theorem is asserted. See \path{ARCHIVE_NOTE_N31_CERTIFICATE_REVISED.md}. No frozen certificate was edited.
'''
}
for name,text in texts.items():(R/name).write_text(text.strip()+'\n')
