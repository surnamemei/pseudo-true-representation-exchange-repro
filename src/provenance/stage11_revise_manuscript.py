"""Scientific manuscript revision. Run once against the pre-revision source."""
from pathlib import Path
import shutil
R=Path(__file__).resolve().parent
def replace_between(s,a,b,new):
    i=s.index(a);j=s.index(b,i);return s[:i]+new+s[j:]
def main():
    paths=['paper/main_tsp_reviewfriendly.tex','paper/sections/body_reviewfriendly.tex','supplement/supplement.tex','supplement/S1_complete_local_proof.tex','supplement/S2_interval_certificate_protocol.tex','supplement/S4_extended_generality.tex','supplement/S7_archive_notes.tex']
    backup=R/'stage11_before_revision';backup.mkdir(exist_ok=True)
    for p in paths:
        out=backup/p;out.parent.mkdir(parents=True,exist_ok=True)
        if not out.exists():shutil.copy2(R/p,out)
    p=R/paths[0];p.write_text(p.read_text().replace('Critical Branch Exchange','Branch Exchange'))
    p=R/paths[1];s=p.read_text()
    s=replace_between(s,r'\begin{abstract}',r'\end{abstract}',r'''\begin{abstract}
A finite multitone record fitted with fewer undamped tones can change its best lower-order representation as an omitted component grows. For a specified symmetric three-tone family with $N=21$, finite-sum analysis with interval-verified nondegeneracy and complete tangent classification proves a global exchange between two separated fits for sufficiently small strong-tone half-spacing $\Delta$. With $z=N\Delta$, the crossing obeys $\epsilon_c(z)=\lambda_{21}z^2+c_{21}z^4+O(z^6)$. A recurrence argument covers competing coalescent limits, including diverging fitted amplitudes. An arbitrary-phase tangent formula explains how coherent projection terms reshape the competing minima. Separate full-domain interval calculations certify global exchange at $z=2$ for $N=21,31$. Numerical exchange persists to two DFT bins; the entire connecting region is not certified. Noise selects between distinct branch pseudo-true centers. The small-spacing global theorem has an existential radius, and finite-width amplitude certification remains unresolved.
''')
    s=s.replace('a local $N=21$ analytical crossing law through $z^4$','an existential small-spacing global $N=21$ crossing law through $z^4$')
    s=s.replace('(iii) numerical cross-length, asymmetry, and noise consequences.','(iii) an exact phase-dependent tangent mechanism, targeted resolution tests, and numerical noise consequences.')
    s=s.replace('The small-$z$ theorem and $z=2$ certificates are not joined by a certified global continuation.','The global small-$z$ neighborhood has no explicit radius and is not joined to $z=2$ by a certified continuation.')
    s=s.replace('Variable projection of linear amplitudes is standard',r'''Pseudo-true targets and inference under misspecification are established by White and subsequent signal-processing theory \cite{white1982,fortunati2017}. They supply the appropriate interpretation of a deliberately lower-order fit; the generating strong pair need not be its estimation target. Variable projection of linear amplitudes is standard''')
    old='Classical frequency-estimation thresholds and ambiguity are established'
    s=s.replace(old,r'''Statistical resolution criteria compare source separation with estimation uncertainty \cite{smith2005}; clustered super-resolution theory quantifies recovery of the generating nodes \cite{batenkov2021}. Here the noiseless task is to choose a globally best representation at a deliberately smaller order. A DFT-bin separation does not establish a statistical resolution limit. Classical frequency-estimation thresholds and ambiguity are established''')
    # Remove a repetitive positioning paragraph to make space for actual results.
    s=replace_between(s,'The approximation setting also fixes a distinction',r'\section{Finite-record formulation}', '')
    s=s.replace('Our analytical family fixes $b=10$ and the weak-tone phase at $\pi$:',r'''Let $\Delta>0$ be the true physical strong-tone half-spacing, with $\omega_1=-\Delta$, $\omega_2=\Delta$, and $\omega_3=b/N$. Our phase-parametrized record is $x_t=2\cos(zs_t)+\epsilon e^{i\phi}e^{ibs_t}$, where $z=N\Delta$ and $\phi\in[0,\pi]$. The global small-spacing theorem fixes $b=10$ and $\phi=\pi$:''')
    s=s.replace('The strong tones have equal unit amplitudes at physical frequencies $\pm z/N$.','The strong tones have equal unit amplitudes at physical frequencies $\omega_{1,2}=\pm z/N$; fitted physical frequencies remain $\nu_j$ and normalized fitted coordinates $u_j=N\nu_j$.')
    s=s.replace(r'\section{Local small-spacing branch exchange}',r'\section{Global small-spacing branch exchange}')
    s=s.replace('[Two separated local branches]','[Two globally competing small-spacing branches]')
    s=s.replace('This is a local statement; other fits may have lower loss.','Every globally optimal fit in this scaled neighborhood is one of A or B, modulo permutation. At equal cost both are global; away from it the lower-cost branch is the unique global fit. The assertion uses the interval tangent classification and the limiting-chart exhaustion below.')
    s=s.replace('[Transverse local critical law]','[Transverse small-spacing global law]')
    s=s.replace('This continuation establishes two regular competing branches, not global optimality for every small $p$.','Globality additionally requires the following exhaustion of competing limits.')
    at='To see the quartic term, let $R_j(\lambda)$'
    ins=r'''\subsection{Excluding every leading competitor}
Directed interval root isolation over one complete frequency period finds 18 minima and 18 maxima of $R_N(v,\lambda_{21})$. A/B alone have cost $0.05876262330\ldots$; every other stationary minimum exceeds this by more than $0.06444970$. The removable central singularity is evaluated after factoring $A(v)=v^4\widetilde A(v)$ and $B(v,\lambda)=v^2\widetilde B(v,\lambda)$. The finite interval cover leaves no unresolved cell.

Fitting the generating strong frequencies $\pm z$ with optimized amplitudes gives the leading coefficient
\begin{equation}\label{eq:strong-control}
C_{\rm strong}=\lambda^2\{N-D(b)^2/N-D'(b)^2/S_2\}
\end{equation}
equal to $0.08550038683\ldots$. Allowing both fitted frequencies to approach zero at arbitrary rates gives the larger limiting class $\alpha+\beta s+\gamma s^2$, with $\alpha,\beta\in\mathbb C$ and $\gamma\in\mathbb R$. Its best coefficient is
\begin{equation}\label{eq:coal-control}
\begin{aligned}
C_{\rm coal}=C_{\rm strong}
-\lambda^2\frac{\{D''(b)+S_2D(b)/N\}^2}{S_4-S_2^2/N}
\end{aligned}
\end{equation}
or $0.07681023120\ldots$, strictly above A/B. This attainable control includes a reoptimized strong-pair scale.

To see why no divergent-amplitude chart is omitted, every two-tone fit satisfies the order-two recurrence with nodes $\xi_j=e^{iu_j/N}$. A fit with $O(p)$ residual obeys $y=2+ph$ with bounded $h$. The limiting recurrence forces $(1-\xi_1)(1-\xi_2)=0$. With only one central limiting node, the rank-five real tangent derivative gives the central-plus-satellite class already classified. With both nodes central, the recurrence forces $u_1u_2=O(p)$ and $\Delta_t^2h=2\gamma/N^2$ in the limit, yielding precisely the quadratic class above. No bound on fitted amplitudes is assumed. Compactness of the confluent frequency closure and the strict margins force all small-$p$ global minimizers into the A/B wells; their positive desingularized Hessians give unique minima there. Supplementary Sec.~S1 supplies the full argument. This is an existential global result, not an explicit small-$z$ radius.

\subsection{Crossing coefficient}
'''
    s=s.replace(at,ins+at)
    s=s.replace('The solid local $N=21$ asymptotic law','The solid small-spacing $N=21$ asymptotic law')
    # Explicit code predicate in the main paper, full proof in supplement.
    marker='Every terminal cell that fails this strict test'
    s=s.replace(marker,r'''For cell radii $r_m,r_h$, put $c_j=\|w_j(h_c)\|$, $a_j=\|w'_j(h_c)\|$, $T_k=\|t^k\|$, and $L=(N-1)/2$, where $w_1=\cos(ht)$ and $w_2=t\operatorname{sinc}(ht)$. The implemented bound is
\begin{equation}\label{eq:explicit-theta}
\begin{aligned}
\delta_1&=a_1r_h+T_2r_h^2/2,\quad
\delta_2=a_2r_h+T_3r_h^2/6,\\
\theta&=\sqrt2 Lr_m+
\left[\sum_{j=1}^2\{\delta_j/(c_j-\delta_j)\}^2\right]^{1/2},
\end{aligned}
\end{equation}
provided $c_j>\delta_j$. This follows from normalized-column displacement and $|\operatorname{sinc}''|\le1/3$; Supplementary Sec.~S2 derives each term and the inner Taylor predicate. '''+marker)
    s=replace_between(s,'The bracket centers are $0.2481906301722774$',r'\begin{figure*}[t]',r'''The historical N=21 cases $\epsilon_c\pm0.035$ also have complete pointwise interval certificates, as stated in Fig.~\ref{fig:geometry}. They do not certify the amplitude interval between them. A new attempt on $[\epsilon_c-0.03,\epsilon_c+0.03]$, split into twelve width-$0.005$ slabs, failed the tested parameter-interval enclosure predicates. The finite-width extension remains open. The smallest historical outer margin comes from a cell centered at $u\approx(-0.83755,2.31938)$: its actual center cost is $1.27374$, far above the incumbent $0.81181$, and an interval-gradient subdivision excludes every stationary point in the cell. The small margin measures enclosure slack rather than a third basin.

''')
    s=replace_between(s,r'\section{Numerical generality}',r'\section{Noise consequence}',r'''\section{Phase mechanism and resolution regime}\label{numerical-generality-checks}
\subsection{Arbitrary weak-tone phase}
For $q=-s^2+\lambda e^{i\phi}e^{ibs}$, put $c=\cos\phi$, $d=\sin\phi$, $q_0=-S_2+\lambda cD(b)$, $q_1=-2\lambda cD'(b)$, and $q_3=D''(v)+\lambda cD(b-v)$. With $A,B$ as in \eqref{eq:review-AB} using these right-hand sides, define
\begin{equation}\label{eq:phase-H}
C=N-D(v)^2/N,\quad
H=\lambda d\{D(b-v)-D(v)D(b)/N\}.
\end{equation}
Eliminating both real Gram blocks gives the exact general-phase formula
\begin{equation}\label{eq:phase-R}
\begin{aligned}
R={}&S_4+2\lambda cD''(b)+N\lambda^2
-q_0^2/N-q_1^2/(4S_2)\\
&-\lambda^2d^2D(b)^2/N-B^2/A-H^2/C.
\end{aligned}
\end{equation}
The in-phase cross term mixes the pair's curvature with $D(b-v)$; the quadrature term $H^2/C$ reshapes that gain. At $\phi=\pi$, the latter vanishes and the weak component opposes the central curvature. Thus phase changes the finite-record competition itself, rather than merely relabelling a fit.

A 65-point tangent map and interval Krawczyk checks at $\phi/\pi=0.4,0.5,0.75,1$ find separated positive-curvature crossings, with $\lambda_c\approx0.06216541,0.06928102,0.07207353,0.06598041$. Regularity implies open phase neighborhoods, excluding a measure-zero explanation. Numerically the two wells merge at $\phi_*/\pi\approx0.33928149$, $v_*\approx5.11554$, $\lambda_*\approx0.05812872$, where $R_v=R_{vv}=R_{vvv}=0$ and $R_{vvvv}>0$. The tracked A satellite passes through zero at $\phi/\pi\approx0.46495169$; the separated chart is invalid at that isolated phase. Hence the numerical noncoalescent phase family connected to $\pi$ ends at this central collision, while a second separated segment reaches the cusp. Neither complete interval is globally certified. Supplementary Sec.~S4 distinguishes the numerical endpoint solves from representative interval root checks.

\begin{figure*}[t]
\centering\includegraphics[width=\textwidth]{figures/fig7_phase_resolution.pdf}
\caption{(a) General-phase tangent equal-cost curve (numerical); squares mark interval-verified local crossings, and the dashed line marks the numerical central-chart collision. (b) Finite-spacing crossing at seven requested spacings. Only the filled $z=2$ marker has full-domain interval globality; other markers show numerical full-torus ordering. The omitted amplitude slightly exceeds one at $z=5$.}\label{fig:persistence}
\end{figure*}

\subsection{Beyond the sub-bin baseline}
The strong-tone separation is $z/\pi$ DFT bins, or $0.63662$ bins at the certified baseline $z=2$. At $z=2,2.5,\pi,3.5,4,5,2\pi$, numerical crossings occur at $\epsilon_c\approx0.24819,0.37470,0.56005,0.67053,0.82742,1.05236,0.98343$. A $512\times512$ full-torus grid, up to 200 grid-minimum refinements, 32 random starts, a strong-pair initialization, and multiprecision stationary solves found A/B as the lowest two distinct minima at every point. The observed third-minimum margins range from $0.88296$ to $10.5029$. Additional local interval checks verify positive A/B roots at $z=5$ and the decimal approximation to $2\pi$; they are not full-domain proofs.

No disappearance boundary is observed through two bins. A direct continuation jump from four to five initially failed; smaller steps recovered the branch across a region of small B curvature. The solver failure therefore did not locate a physical boundary. Refinement initialized at the generating strong pair reaches A throughout this sweep. At $z=5$ the designated omitted component is no longer weaker than the unit strong tones, and at fixed $b=10$ the positive strong tone approaches it as $z$ grows. The results concern this specified family and are not a universal resolution law.

\subsection{Other record lengths and the quartic comparison}
Table~\ref{tab:crossN} retains the five sampled record lengths at $z=2$; only $N=21,31$ have full-domain finite-spacing certificates. The N=21 quartic prediction at $z=2$ is $0.2476376557484487$, compared with certified $0.2481906301722774$: absolute error $0.00055297442$, relative error $0.2228023\%$. This pointwise accuracy does not bound the remainder on an interval. The secondary numerical N=31 prediction is $0.24857675463874848$ ($0.227313\%$ error). Its frozen archive had carried N=21 coefficients into eight comparison fields; the corrected provenance and clean proof-dependency audit are in Supplementary Sec.~S7. Prior one-factor amplitude/location tests remain in the data archive.
\input{sections/crossN_table.tex}

''')
    s=s.replace('Crosses are within-branch means; the generating strong pair is shown for reference.','Crosses are within-branch means; black plus signs mark the two noiseless branch optima.')
    s=replace_between(s,'For the ordered fitted-frequency pair',r'\section{Discussion}',r'''For the ordered fitted-frequency pair $\widehat u$, branch label $C\in\{A,B\}$, probabilities $p_j$, and conditional means $\mu_j$, the law of total covariance gives
\begin{equation}\label{eq:main-15}
\begin{aligned}
\operatorname{tr}\operatorname{Cov}(\widehat u)
={}&\sum_j p_j\operatorname{tr}\operatorname{Cov}(\widehat u\mid C=j)\\
&+p_Ap_B\|\mu_A-\mu_B\|^2.
\end{aligned}
\end{equation}
At $\eta=0$, the within-branch variance trace falls from $0.2094$ at 20 dB to $0.00248$ at 40 dB, while the between-branch term remains about $40.1$--$40.6$ and the fitted centers stay roughly $12.7$ normalized units apart. The appropriate deterministic centers are the two noiseless branch optima. A covariance conditioned on one selected branch omits the selection term; this observation does not contradict misspecified estimation bounds \cite{white1982,fortunati2017}. No MSE relative to the generating strong pair is used to interpret this experiment.

''')
    s=s.replace('The analytical theorem uses equal strong amplitudes, $N=21$, weak location $b=10$, phase $\pi$, and an unquantified small-$z$ radius.','The global small-spacing theorem uses equal strong amplitudes, $N=21$, weak location $b=10$, phase $\pi$, and an unquantified small-$z$ radius.')
    s=s.replace('A connected global certificate in $z$ is the clearest extension; arbitrary $(K+1)\to K$ underfit, damped tones, and other phase configurations require separate analysis.','The phase cusp and resolution sweep have no complete global interval proof, and the finite-width amplitude attempt failed its enclosure checks. Arbitrary $(K+1)\to K$ underfit and damped tones require separate analysis.')
    s=s.replace('we proved a local noncoalescent branch exchange and its quadratic critical-amplitude law with quartic correction.','we proved an existential small-spacing global noncoalescent branch exchange and its quadratic amplitude law with quartic correction, after classifying all leading limiting fits.')
    s=s.replace('Numerical generality and noise studies show the associated two-mode estimator behavior within their tested settings.','An exact arbitrary-phase tangent formula and numerical continuation explain the coherent phase mechanism and show exchange beyond one-bin spacing. The wider finite-amplitude certificate and the connected bridge to the small-spacing theorem remain open.')
    p.write_text(s)
    p=R/'paper/sections/evidence_table.tex';s=p.read_text().replace("three evidence categories","evidence categories")
    s=s.replace('Two separated strict local branches; transverse crossing','Global A/B exchange for sufficiently small spacing')
    s=s.replace('Desingularized tangent problem, interval-verified nondegeneracy, analytic implicit-function theorem.','Full-period interval tangent classification, recurrence-based limiting-chart exhaustion, analytic implicit-function theorem.')
    s=s.replace('Cross-length and one-factor asymmetry persistence','Phase/resolution and cross-length persistence')
    s=s.replace('Multiprecision solves and full-torus grid refinements.','General-phase formula; representative interval roots; numerical full-torus grid refinements.')
    p.write_text(s)
    p=R/'supplement/supplement.tex';s=p.read_text().replace('Critical Branch Exchange','Branch Exchange')
    s=s.replace('The analytical proof concerns the symmetric $N=21$ family and an existential small-$z$ neighborhood.','The analytical global proof concerns the symmetric $N=21$ family and an existential small-$z$ neighborhood, using complete interval tangent classification.')
    p.write_text(s)
    p=R/'supplement/S3_certificate_tables.tex';s=p.read_text().replace('in A/B boxes;','in the A box first, then the B box; derivatives are with respect to $u=N\\nu$ (normalized angular-frequency coordinates), so $H_{11}$ and $\\det H$ have units of loss per $u^2$ and per $u^4$;')
    p.write_text(s)
    p=R/'supplement/S5_noise_search_details.tex';s=p.read_text();i=s.index('Relative to $u_0=');j=s.index('These are numerical outcomes',i);s=s[:i]+'The noiseless branch optima define the relevant pseudo-true centers; the main paper uses total variance, not error to the generating strong pair. '+s[j:];p.write_text(s)
if __name__=='__main__':main()
