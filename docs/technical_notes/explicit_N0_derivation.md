# Stage 15 effective-threshold certificate

The continuum theorem alone gives an existential $N_0$. This Stage 15 certificate instead treats $t=N^{-2}$ as an **interval parameter** and checks every old full-line competitor predicate for all $t$ in one interval $[0,N_0^{-2}]$. A passing test therefore covers every odd $N\ge N_0$ at once. The certificate is for $b=10$, $\phi=\pi$ and the small-$z$ tangent problem only.

## Exact parameterization and rigorous truncation

The finite centered kernel is $d_t(q)=d(q)H(q^2t/4)$, with $H(w)=\sqrt w/\sin\sqrt w$. The checker obtains the first nine $H$ coefficients by an exact rational recurrence and encloses all omitted terms. For $|q|\le110$ and $k\le46$, it uses a directed $\pm10^{-30}$ remainder for $d_t^{(k)}$. Here is a conservative analytic justification for the allowed $N$ range.

On the complex circle $|u|=1/2$,

$$
\left|\frac{\sin u}{u}-1\right|
\le \frac{\sinh(1/2)}{1/2}-1=:\eta<0.043,
\qquad |H(u^2)|\le(1-\eta)^{-1}.
$$

Cauchy gives $|h_m|\le(1-\eta)^{-1}4^m$. For a complex-$q$ circle of radius 55 about any real $q\in[-110,110]$, $|q|\le165$ and $|d(q)|\le e^{27.5}$. Thus the kernel truncation tail is bounded on the circle by

$$
B(N)=\frac{e^{27.5}}{1-\eta}\,
\frac{(165/N)^{18}}{1-(165/N)^2}.
$$

Cauchy's derivative factor $k!/55^k\le1$ for $0\le k\le46$. Consequently $B(N)<10^{-30}$ suffices for the $k\le46$ derivative enclosure. The first allowed odd $N$ from this **chosen** analytic remainder inequality is $35{,}377$ (the real threshold is approximately $35{,}375.016$). The checker may certify a larger threshold because its global predicates impose additional conditions.

On $|v|\le1$, the checker uses Taylor models for the analytic deflations $A_t(v)/v^4$ and $B_t(v)/v^2$. Their coefficients use the same exact rational $h_m$ recurrence. A complex-$v$ circle of radius 10 bounds each deflated analytic function by less than $10^6$ for the certified $t$ and $\lambda$ intervals: $|d_t^{(k)}|$ is bounded using $|s|\le1/2$ in the centered integral/analytic-kernel formula, $\mu_{2,t}$ stays away from zero, and $|\lambda|<0.067$. The Taylor tail after degree 40, including two derivatives on $|v|\le1$, is below $10^{-25}$. The $m\ge9$ omitted $H$ contributions are many orders smaller. The code widens all three deflated outputs by $\pm10^{-25}$.

## Uniform root and global predicates

For each $t$ in the certified interval, a parameterized directed Krawczyk map encloses a unique joint A/B stationary and equal-cost root. Its signs preserve positive branch curvatures and Gram factors, nonzero satellite amplitudes, and transversality. The Krawczyk image lies inside the two root-containing cells in the frozen Stage 14 partition.

The exact-decimal Stage 14 partition covers $[-100,100]$ with 27,708 terminal cells. Each cost cell must retain cost greater than a feasible A trial; each gradient cell must retain a fixed derivative sign; each monotonicity cell must retain curvature and endpoint-derivative signs. The two remaining cells contain the isolated A/B roots. This is stronger than comparing just the narrowest cost margin: it preserves the actual exhaustive global proof predicate by predicate. The coalescent loss is checked separately.

For the expanding finite principal period $|v|\le\pi N$, discrete summation by parts gives conservative bounds

$$
|d_N(v)|\le\frac{2\pi}{|v|},\qquad
|d_N'(v)|\le\frac{\pi}{|v|},\qquad
|d_N''(v)|\le\frac{\pi}{2|v|}.
$$

For $|v|\ge100$, the $d_N(b-v)$ bound uses the principal circular separation at least $|v|-b$. Inserting these bounds into the Schur-complement loss yields the certified positive finite-period tail margin. The interval checker records this margin and the coalescent margin separately.

## Threshold and scope

The machine-readable result is `explicit_N0_certificate.json`; its supporting all-$t$ run is `stage15_interval_trial_N35377.json`, produced by `stage15_uniform_interval.py`. The selected threshold is

$$N_{0,\mathrm{interval}}=35{,}377.$$

At this floor the all-$t$ Krawczyk contraction is below $0.057402$, all 27,708 partition predicates pass, and the coalescent and expanding-tail lower margins are respectively $0.0008511670$ and $0.0031356601$. The smallest cost-cell margin is $5.08278\times10^{-9}$; the smallest derivative-sign margin is $1.02446\times10^{-9}$. The preceding odd $N=35{,}375$ violates the **chosen** analytic $10^{-30}$ truncation inequality; no claim is made that the theorem itself fails there. A stand-alone closed-form $N_{0,\mathrm{analytic}}$ from global rational-loss sup-norm constants was not obtained, so the JSON records it as null. A numerical fit of $\lambda_N$ is not used in this certificate.

The following entries are conservative lower bounds from the continuum certificate and from the uniform finite-$N$ interval replay. Their difference $E_i^*$ bounds erosion of the **certified inequalities**, which is sufficient for transfer. Define $E_i(N)=E_i^*$ for every odd $N\ge35{,}377$; this is an explicit uniform error budget (not a sharp pointwise rate). Each finite lower bound exceeds half the corresponding continuum margin.

| Inequality | Continuum lower $m_i$ | Finite uniform lower | $E_i^*$ |
|---|---:|---:|---:|
| branch curvature | $1.33682\times10^{-4}$ | $1.30049\times10^{-4}$ | $3.63255\times10^{-6}$ |
| $|\beta|$ | $0.0644229455$ | $0.0644223994$ | $5.46149\times10^{-7}$ |
| crossing slope | $0.1078175136$ | $0.1078144651$ | $3.04849\times10^{-6}$ |
| regular competitor cost | $6.28143\times10^{-9}$ | $5.08278\times10^{-9}$ | $1.19864\times10^{-9}$ |
| coalescent cost | $0.000851168167$ | $0.000851167020$ | $1.14761\times10^{-9}$ |

Existing individual certificates cover $N=11,15,21,31,41$; no claim about every intervening odd $N$ is made.

The threshold is a sufficient proof threshold, not a physical regime change or an estimate of the smallest true $N$ with A/B exchange. A direct interval sweep through all smaller odd $N$ would be disproportionate; the existing five anchors remain the only individual finite-$N$ tangent certificates below it.
