# Stage 13: the proposed positive lower transition does not exist

For the tangent problem $N=21$, $b=10$, $\phi=\pi$, the proposed inference fixes a separated satellite frequency while taking $\lambda\downarrow0$. The optimizing frequency moves toward the coalescent endpoint. Positive cost at a *fixed* separated frequency at $\lambda=0$ therefore does not imply a positive crossing with the coalescent cost.

## Exact endpoint calculation

Put $S_j=\sum_t s_t^j$, $E_2=s^2-S_2/N$, $U_3=s^3-(S_4/S_2)s$, $V_2=\|E_2\|^2$, and $V_3=\|U_3\|^2$. Let

$$
M=\langle U_3,\sin(bs)\rangle=D'''(b)+(S_4/S_2)D'(b),\quad
Q=D''(b)+S_2D(b)/N.
$$

After removing the central complex amplitude and real central-frequency tangent, the satellite vector divided by its leading $-v^2/2$ factor is $E_2+i(v/3)U_3+O(v^2)$. Scalar projection gives

$$
\lim_{v\to0}R_{21}(v,\lambda)=C_{\rm coal}(\lambda),\qquad
\partial_vR_{21}(0,\lambda)=\frac23\gamma(\lambda)\lambda M,\quad
\gamma(\lambda)=-1+\lambda Q/V_2.
$$

Directed intervals in the machine certificate give $M=-0.24033790138997\ldots<0$, $Q=0.47984461105929\ldots>0$, $V_2,V_3>0$, and $\lambda_\gamma=V_2/Q=0.24038264103490\ldots$. Thus $v=0$ is not a minimizer for any $0<\lambda<\lambda_\gamma$. The existing interval-certified A/B crossing lies at $\lambda_{21}=0.06598041112699\ldots<\lambda_\gamma$.

There is a global-in-$\lambda$ witness. The coalescent cost is $C_{\rm coal}(\lambda)=17.64367452454865\ldots\lambda^2$. The regular fit with satellite fixed at $v=b$ has constant cost $R(b,\lambda)=0.10362265533211\ldots$ because it absorbs the weak-tone term. Hence $R(b,\lambda)<C_{\rm coal}(\lambda)$ for $\lambda>0.07663600191508\ldots$. This threshold is below $\lambda_\gamma$. The two arguments cover **all $\lambda>0$**; there is no positive coalescent global regime. Outward interval bounds for every quoted comparison are in lower_transition_certificate.json.

At $\lambda=0$, the endpoint is the unique global tangent minimizer. For any regular $v\ne0\pmod{2\pi N}$, a quadratic sequence cannot lie in $\operatorname{span}\{1,s,e^{ivs}\}$: second finite differences would make a constant equal a nonconstant exponential. Near the endpoint, $R_{vv}(0,0)=2V_3/9>0$. Compactness of the frequency circle and the implicit-function theorem yield a unique global regular branch for all sufficiently small positive $\lambda$:

$$
v_0(\lambda)=\frac{3M}{V_3}\lambda+O(\lambda^2)
=-99.26134693376595\ldots\lambda+O(\lambda^2),
$$
$$
C_{\rm coal}(\lambda)-R_0(\lambda)
=\frac{M^2}{V_3}\lambda^2+O(\lambda^3)
=7.95208793706778\ldots\lambda^2+O(\lambda^3)>0.
$$

At $\lambda_{21}$, the existing tangent interval certificate gives $R_A=R_B\approx0.05876262329535$, while $C_{\rm coal}\approx0.07681023119888$, with endpoint gap $0.01804760790353\ldots$. This endpoint is handled by C7, **not** the C6 margin over other regular stationary minima.

The certificate excludes a positive coalescent crossing. It does not certify which *regular* branch is globally best throughout $(0,\lambda_{21})$. In particular, connecting the near-zero branch $v_0$ to A at the certified crossing remains numerical.
