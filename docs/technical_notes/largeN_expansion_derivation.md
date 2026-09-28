# Exact centered-grid expansion for the tangent loss

Fix the paper's $b=10$, $\phi=\pi$ family. Write $t=N^{-2}$, $d(q)=2\sin(q/2)/q$, and $d_t(q)=D_N(q)/N$. All limits at $q=0$ are removable. The exact identity

$$
d_t(q)=d(q)H(q^2t/4),\qquad H(w)=\frac{\sqrt w}{\sin\sqrt w}
$$

gives, with rational coefficients,

$$
d_t=d+t d_2+t^2d_4+O(t^3),\qquad
d_2(q)=\frac{q^2d(q)}{24},\quad
d_4(q)=\frac{7q^4d(q)}{5760}.
$$

The next coefficient is $d_6(q)=31q^6d(q)/967680$. The coefficients of $H$ follow exactly from $h_0=1$ and $h_m=-\sum_{j=1}^m(-1)^j h_{m-j}/(2j+1)!$. For $k=0,1,2,3,4$ (and any fixed higher derivative), differentiation gives

$$
d_t^{(k)}(q)=d^{(k)}(q)+\frac{t}{24}
\sum_{j=0}^{\min(k,2)}{k\choose j}(2)_j q^{2-j}d^{(k-j)}(q)
+\frac{7t^2}{5760}\sum_{j=0}^{\min(k,4)}{k\choose j}(4)_j q^{4-j}d^{(k-j)}(q)
+O(t^3).
$$

Here $(m)_j=m!/(m-j)!$. This formula explicitly covers $D_N'/N$, $D_N''/N$, $D_N'''/N$ and the fourth derivative needed for $R_{vv}$. In particular, the $t$ coefficient is $(q^2d)'/24=(2qd+q^2d')/24$ for the first derivative; $(2d+4qd'+q^2d'')/24$ for the second; and $(6d'+6qd''+q^2d''')/24$ for the third.

## Uniform compact remainder

For any fixed real $|q|\le Q$, take the complex-$t$ Cauchy radius $r=(Q+1)^{-2}$ and a complex-$q$ circle of radius one. Then $|q\sqrt t/2|\le1/2$. Since $|\sin u/u-1|\le \sinh(1/2)/(1/2)-1=:\eta<0.043$, $|H|\le(1-\eta)^{-1}$ there; the integral representation gives $|d(q)|\le e^{1/2}$ on the $q$ circle. Cauchy's estimates therefore give, for $N\ge\sqrt2(Q+1)$,

$$
\sup_{|q|\le Q}\left|d_t^{(k)}-d^{(k)}-t d_2^{(k)}-t^2d_4^{(k)}\right|
\le \frac{2k!e^{1/2}}{1-\eta}\left(\frac{Q+1}{N}\right)^6.
$$

This is a conservative but explicit derivative bound for every fixed $k$. The kernel is analytic in $t$ and centered symmetry removes odd powers of $N^{-1}$.

## Exact moments

Let $\mu_{j,N}=S_j/N=N^{-1}\sum s_n^j$. Direct centered power sums give

$$
\begin{aligned}
\mu_{2,N}&=\frac1{12}-\frac{t}{12},\\
\mu_{4,N}&=\frac1{80}-\frac{t}{24}+\frac{7t^2}{240},\\
\mu_{6,N}&=\frac1{448}-\frac{t}{64}+\frac{7t^2}{192}-\frac{31t^3}{1344}.
\end{aligned}
$$

Odd moments vanish identically. These exact polynomials independently confirm the even-power rule.

## Explicit loss coefficients

Use index $j=0,1,2$ for powers $t^j$: $d_j=(d,d_2,d_4)$, $m_j=(1/12,-1/12,0)$, $n_j=(1/80,-1/24,7/240)$, and $u_j=12$ (the first three coefficients of $1/\mu_{2,N}=12/(1-t)$). Set

$$
q_{0j}=-m_j-\lambda d_j(b),\qquad q_{1j}=2\lambda d_j'(b).
$$

Every sum below is truncated to nonnegative indices and total index $j$. The exact coefficient operations are

$$
\begin{aligned}
A_j={}&\mathbf1_{j=0}-\sum_{a+c=j}d_a(v)d_c(v)
-\sum_{a+c+e=j}u_a d_c'(v)d_e'(v),\\
B_j={}&d_j''(v)-\lambda d_j(b-v)-\sum_{a+c=j}d_a(v)q_{0c}
+\frac12\sum_{a+c+e=j}u_a d_c'(v)q_{1e},\\
C_j={}&n_j-2\lambda d_j''(b)+\mathbf1_{j=0}\lambda^2
-\sum_{a+c=j}q_{0a}q_{0c}
-\frac14\sum_{a+c+e=j}u_aq_{1c}q_{1e}.
\end{aligned}
$$

Writing $W_j=\sum_{a+c=j}B_aB_c$, the reduced loss expansion $R_N/N=R_\infty+tG_2+t^2G_4+O(t^3)$ has

$$
\begin{aligned}
R_\infty&=C_0-W_0/A_0,\\
G_2&=C_1-W_1/A_0+W_0A_1/A_0^2,\\
G_4&=C_2-W_2/A_0+W_1A_1/A_0^2
+W_0(A_2/A_0^2-A_1^2/A_0^3).
\end{aligned}
$$

These formulas hold directly on compact regular sets with $A_0>0$. At $v=0$, apply them to the analytic deflations $A/v^4$ and $B/v^2$; this gives the continuous endpoint loss and the same uniform expansion. Differentiating the analytic coefficient functions gives $R_v/N$, $R_{vv}/N$, and $R_\lambda/N$ expansions with coefficients $\partial_vG_j$, $\partial_{vv}G_j$, and $\partial_\lambda G_j$ and compact $O(t^3)$ remainders.

The expansion alone does not provide an effective $N_0$: one must also bound the rational loss and its derivatives on every competitor cell and the expanding finite-frequency tail.
