# Analytic continuation of the continuum crossing

Fix $b=10$, $\phi=\pi$, let $t=N^{-2}$, and write $r(v,\lambda,t)=R_N(v,\lambda)/N$ through the analytic midpoint-kernel continuation. At $t=0$, the certified continuum system has a unique local solution $(v_A,v_B,\lambda_\infty)$ and an invertible Jacobian. Its diagonal branch curvatures are strictly positive and its equal-cost slope is strictly positive; see `stage14_continuum_local_interval.json`.

Define

$$
F(v_A,v_B,\lambda,t)=
\begin{pmatrix}
r_v(v_A,\lambda,t)\\r_v(v_B,\lambda,t)\\
r(v_A,\lambda,t)-r(v_B,\lambda,t)
\end{pmatrix}.
$$

The analytic implicit-function theorem gives a unique nearby root $x(t)=x_0+x_1t+O(t^2)$, where $x_1=-F_x(x_0,0)^{-1}F_t(x_0,0)$. The scalar coefficient formulas are especially transparent. Put $S=\partial_\lambda(R_{\infty,A}-R_{\infty,B})\ne0$ and use $G_2$ from `largeN_expansion_derivation.md`. At stationary branch roots,

$$
\begin{aligned}
a_\lambda&=-\frac{G_2(v_A,\lambda_\infty)-G_2(v_B,\lambda_\infty)}{S},\\
a_A&=-\frac{R_{v\lambda,A}a_\lambda+\partial_vG_2(v_A,\lambda_\infty)}{R_{vv,A}},\\
a_B&=-\frac{R_{v\lambda,B}a_\lambda+\partial_vG_2(v_B,\lambda_\infty)}{R_{vv,B}}.
\end{aligned}
$$

Analytic differentiation, followed by high-precision evaluation at the certified continuum root, gives

$$
\begin{aligned}
\lambda_N&=\lambda_\infty-0.23156010864980659834\ldots N^{-2}
-0.28314327490967740462\ldots N^{-4}+O(N^{-6}),\\
v_{A,N}&=v_{A,\infty}-40.49563646218550239\ldots N^{-2}+O(N^{-4}),\\
v_{B,N}&=v_{B,\infty}+19.15838762577246328\ldots N^{-2}+O(N^{-4}).
\end{aligned}
$$

For the second crossing coefficient, let $F_{\rm fix}(\lambda,t)=r(v_{A,\infty},\lambda,t)-r(v_{B,\infty},\lambda,t)$ hold the branch frequencies at their continuum values, and put $H_j=R_{vv,j}$ and $G_j=R_{v\lambda,j}a_\lambda+R_{vt,j}$. Then

$$
b_\lambda=-\frac{\tfrac12F_{{\rm fix},\lambda\lambda}a_\lambda^2+F_{{\rm fix},\lambda t}a_\lambda
+\tfrac12F_{{\rm fix},tt}-G_A^2/(2H_A)+G_B^2/(2H_B)}{F_{{\rm fix},\lambda}}.
$$

This is the source of the displayed $N^{-4}$ coefficient; it includes branch-frequency drift. These coefficients are **derived by differentiation**, not fitted from the nine validation rows. The decimals are high-precision evaluations, while the certified theorem is the analytic local expansion. The finite-$N$ data in `lambdaN_largeN_data.csv` and `lambdaN_residual_scaling.csv` distinguish the five existing interval-certified anchors from four numerical large-$N$ validations.

The corrected residual $N^4[\lambda_N-\lambda_\infty-a_\lambda/N^2]$ changes from $-0.3011078249$ at $N=11$ to $-0.2831452734$ at $N=1001$, approaching the derived $b_\lambda=-0.2831432749\ldots$. The $N^{-4}$ remainder at $N=1001$ is about $-1.99\times10^{-18}$; this numerical agreement is validation, not a uniform proof or an effective threshold.
