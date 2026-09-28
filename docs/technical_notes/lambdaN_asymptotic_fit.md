# The large-$N$ crossing coefficient

The continuum stationary equal-cost solve gives

$$
\lambda_\infty=0.06650697039239553545591632113\ldots,
\quad
v_{A,\infty}=-4.06003021618010689\ldots,
\quad
v_{B,\infty}=12.42390445039620145\ldots .
$$

These are high-precision candidates; stage14_continuum_local_interval.json independently encloses a unique local crossing with positive curvatures and nonzero transverse slope. Globality is recorded separately in continuum_global_certificate.json.

## Derivation of the rate

The centered grid is a midpoint rule. Its exact kernel is

$$
d_N(q)=d(q)\frac{q/(2N)}{\sin(q/(2N))}
=d(q)\left(1+\frac{q^2}{24N^2}
+\frac{7q^4}{5760N^4}+O(N^{-6})\right).
$$

The finite moments obey $\mu_{2,N}=(1-N^{-2})/12$ and
$\mu_{4,N}=1/80-1/(24N^2)+7/(240N^4)$. Thus the reduced tangent loss has an even-power expansion in $t=N^{-2}$ on the regular A/B neighborhoods:
$R_N/N=R_\infty+tR_2+t^2R_4+O(t^3)$.
There is no $1/N$ term.

Let $F(\lambda,t)$ be the optimized A-minus-B tangent cost difference, with the branch stationarity equations solved locally. Its simple continuum zero and nonzero $\partial_\lambda F$ imply an analytic crossing $\lambda(t)$. Differentiating $F(\lambda(t),t)=0$ gives

$$
a=\lambda'(0)=-F_t/F_\lambda
=-0.23156010864980660008316350618\ldots .
$$

At second order the branch-frequency changes matter. If $G_j=R_{v\lambda,j}a+R_{vt,j}$ and $H_j=R_{vv,j}>0$, then

$$
b=-\frac{\tfrac12F_{\lambda\lambda}a^2+F_{\lambda t}a
+\tfrac12F_{tt}-G_A^2/(2H_A)+G_B^2/(2H_B)}{F_\lambda}
=-0.28314327490967780278161639227\ldots .
$$

Consequently, conditional on the globally isolated continuum crossing,

$$
\lambda_N=\lambda_\infty-\frac{0.2315601086498066\ldots}{N^2}
-\frac{0.2831432749096778\ldots}{N^4}+O(N^{-6}).
$$

The derivation is analytic; the displayed coefficients are high-precision numerical evaluations, not yet interval-enclosed.

## Validation and limits

lambdaN_asymptotics.csv combines the existing interval-certified $N=11,15,21,31,41$ crossings with four **numerical validation** solves at $N=101,201,501,1001$. At $N=1001$, $N^2(\lambda_N-\lambda_\infty)=-0.2315603912296\ldots$, approaching the derived $a$. The $N^{-4}$ prediction errors decay rapidly, reaching about $-2.0\times10^{-18}$ at $N=1001$. These data validate the expansion and do not certify a numerical all-$N$ threshold.

The figure continuum_vs_finiteN_figure.pdf distinguishes existing certified finite-$N$ anchors from the four numerical large-$N$ validations. An effective $N_0$ still requires uniform quantitative exclusion bounds.
