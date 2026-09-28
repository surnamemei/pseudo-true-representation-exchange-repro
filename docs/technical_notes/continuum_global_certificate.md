# Continuum global A/B crossing certificate

Model: $b=10$, $\phi=\pi$, normalized $L^2[-1/2,1/2]$ tangent loss $R_\infty(v,\lambda)$. The frequency variable $v$ ranges over the **entire real line**. The sinc sidelobes have indefinitely many small extrema, so isolating every stationary point individually is unnecessary for globality; a rigorous tail bound excludes all sufficiently distant extrema at once.

## Local crossing

The high-precision candidate is

| Quantity | Value |
|---|---:|
| $\lambda_\infty$ | $0.0665069703923955354559\ldots$ |
| $v_{A,\infty}$ | $-4.0600302161801068928\ldots$ |
| $v_{B,\infty}$ | $12.4239044503962014517\ldots$ |
| common cost $R_*$ | $0.0028878765288308924689\ldots$ |
| $\beta_A,\beta_B$ | $0.1274999426\ldots,\ -0.0644229456\ldots$ |
| $R_{vv,A},R_{vv,B}$ | $0.0001336818966\ldots,\ 0.0009603691055\ldots$ |
| $\partial_\lambda(R_A-R_B)$ | $0.1078175136496\ldots$ |

The directed interval Krawczyk calculation in stage14_continuum_local_interval.json proves a unique solution of the two stationarity equations and equal-cost equation in its product box. Its contraction upper bound is approximately $7.65\times10^{-12}<1$; both Gram factors, curvatures, and the crossing slope have strictly positive lower bounds, and the two satellite coefficients have opposite nonzero signs.

## Full real-line exclusion

The interval partition covers $[-100,100]$. The replay file stage14_continuum_global_replay.json verifies **27,708 terminal cells**, with no failed predicate or uncovered gap:

| Predicate | Cells |
|---|---:|
| cost strictly above a feasible A trial | 25,274 |
| derivative excludes stationarity | 1,420 |
| fixed-sign curvature and same-sign endpoint derivatives exclude a root | 1,012 |
| unique A/B stationary-root boxes | 2 |

The positive directed cost-cell margin is at least $6.2814\times10^{-9}$. This deliberately conservative number is a proof margin, not the numerical third-best gap. Only the two *competitive* stationary roots are isolated; higher-loss stationary roots may lie in cost-excluded cells. The Krawczyk root boxes lie strictly inside the corresponding global root boxes, linking the local equality certificate to the full-domain exclusion.

At the coalescent endpoint, $C_{\rm coal,\infty}=0.003739044696188905\ldots$, with a certified gap above the A/B cost of at least $0.0008511681673580$. For $|v|\ge100$, the analytic sinc derivative bounds give
$R_\infty(v,\lambda_\infty)\ge0.00663716576118\ldots$, leaving a certified tail margin of at least $0.00374928923235$. This covers all $v\in\mathbb R$, including every distant sidelobe minimum.

The companion endpoint certificate stage14_continuum_endpoint_interval.json shows $v=0$ is improvable for every positive $\lambda$. Thus the continuum crossing satisfies the analogue of C1–C7. The finite-$N$ limiting-chart exhaustion C8 is a separate argument and is not inferred from continuum optimization.

## Reproduction and scope

Run stage14_continuum.py, stage14_continuum_interval.py, stage14_continuum_global_hybrid.py, stage14_replay_global.py, and stage14_continuum_endpoint_interval.py in that order. The global JSON and partition CSV are machine-readable proof objects. The certificate is for this $b,\phi$ family at the unique crossing in the stated local lambda box; it is not a global phase diagram for every lambda. It does not give an effective finite-$N$ threshold.
