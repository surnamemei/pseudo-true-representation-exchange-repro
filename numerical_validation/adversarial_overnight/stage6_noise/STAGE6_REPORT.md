# Stage 6 — noisy estimator: local vs global behaviour

N = 21, z = 2, b = 10, φ = π; ε = ε_c(1 + η), ε_c = 0.2481906301722774; SNR 20/30/40 dB with the frozen noise convention; **1000 records per (η, SNR)** (39,000 records; the pre-declared reduction rule was not triggered). Evidence level: GLOBAL_NUMERICAL (numerical global search per record; not certified).

## 6A. Solver invariance (before trusting the Monte Carlo)

600 frozen records (η ∈ {−0.1, −0.025, 0, 0.025, 0.1} × 3 SNR × 40). Production (blind global search, no oracle seeds) vs reference (4× denser grid, 4× more random starts, oracle A/B seeds): **0/600 disagreements** in objective (> 1e-8) or label. Archived two-oracle-seed pipeline vs reference: **0/600**. Reference labels: A 294, B 294, OUT 12. No production upgrade was needed.

## 6B. Selection probabilities and errors

| SNR | η | P(A) | P(B) | P(out) | global MSE to noiseless winner | A-conditioned MSE (own root) | B-conditioned MSE (own root) | local sandwich trace A | local sandwich trace B |
|---|---|---|---|---|---|---|---|---|---|
| 20 | -0.200 | 0.790 | 0.111 | 0.099 | 18.37 | 0.2542 | 0.09104 | 0.3993 | 0.09276 |
| 20 | -0.150 | 0.758 | 0.175 | 0.067 | 28.63 | 0.2436 | 0.09406 | 0.3889 | 0.08901 |
| 20 | -0.100 | 0.670 | 0.260 | 0.070 | 42.3 | 0.2521 | 0.09764 | 0.3787 | 0.08537 |
| 20 | -0.075 | 0.600 | 0.350 | 0.050 | 56.84 | 0.2457 | 0.07932 | 0.3737 | 0.0836 |
| 20 | -0.050 | 0.584 | 0.359 | 0.057 | 58.55 | 0.2537 | 0.07201 | 0.3688 | 0.08187 |
| 20 | -0.025 | 0.500 | 0.442 | 0.058 | 72 | 0.2419 | 0.08161 | 0.364 | 0.08017 |
| 20 | +0.000 | 0.464 | 0.489 | 0.047 | 83.25 | 0.2327 | 0.07254 | 0.3592 | 0.07851 |
| 20 | +0.025 | 0.394 | 0.570 | 0.036 | 69.76 | 0.2201 | 0.08278 | 0.3545 | 0.07688 |
| 20 | +0.050 | 0.348 | 0.620 | 0.032 | 61.78 | 0.2293 | 0.07787 | 0.3499 | 0.07529 |
| 20 | +0.075 | 0.299 | 0.678 | 0.023 | 52.25 | 0.2304 | 0.0761 | 0.3454 | 0.07374 |
| 20 | +0.100 | 0.254 | 0.728 | 0.018 | 44.02 | 0.2453 | 0.06906 | 0.3409 | 0.07223 |
| 20 | +0.150 | 0.182 | 0.807 | 0.011 | 31.23 | 0.2393 | 0.06753 | 0.3323 | 0.06931 |
| 20 | +0.200 | 0.090 | 0.898 | 0.012 | 16.42 | 0.2208 | 0.06618 | 0.3239 | 0.06654 |
| 30 | -0.200 | 1.000 | 0.000 | 0.000 | 0.0391 | 0.0391 | — | 0.03993 | 0.009276 |
| 30 | -0.150 | 1.000 | 0.000 | 0.000 | 0.0396 | 0.0396 | — | 0.03889 | 0.008901 |
| 30 | -0.100 | 0.976 | 0.024 | 0.000 | 3.927 | 0.04109 | 0.0137 | 0.03787 | 0.008537 |
| 30 | -0.075 | 0.943 | 0.057 | 0.000 | 9.316 | 0.03842 | 0.0102 | 0.03737 | 0.00836 |
| 30 | -0.050 | 0.830 | 0.170 | 0.000 | 27.62 | 0.03329 | 0.009842 | 0.03688 | 0.008187 |
| 30 | -0.025 | 0.687 | 0.313 | 0.000 | 50.85 | 0.03633 | 0.008894 | 0.0364 | 0.008017 |
| 30 | +0.000 | 0.510 | 0.490 | 0.000 | 82.95 | 0.03584 | 0.007824 | 0.03592 | 0.007851 |
| 30 | +0.025 | 0.299 | 0.701 | 0.000 | 48.58 | 0.03345 | 0.006692 | 0.03545 | 0.007688 |
| 30 | +0.050 | 0.157 | 0.843 | 0.000 | 25.55 | 0.04047 | 0.007431 | 0.03499 | 0.007529 |
| 30 | +0.075 | 0.068 | 0.932 | 0.000 | 11.01 | 0.03338 | 0.006887 | 0.03454 | 0.007374 |
| 30 | +0.100 | 0.032 | 0.968 | 0.000 | 5.177 | 0.04676 | 0.007451 | 0.03409 | 0.007223 |
| 30 | +0.150 | 0.003 | 0.997 | 0.000 | 0.4872 | 0.0573 | 0.006916 | 0.03323 | 0.006931 |
| 30 | +0.200 | 0.000 | 1.000 | 0.000 | 0.006552 | — | 0.006552 | 0.03239 | 0.006654 |
| 40 | -0.200 | 1.000 | 0.000 | 0.000 | 0.003949 | 0.003949 | — | 0.003993 | 0.0009276 |
| 40 | -0.150 | 1.000 | 0.000 | 0.000 | 0.003805 | 0.003805 | — | 0.003889 | 0.0008901 |
| 40 | -0.100 | 1.000 | 0.000 | 0.000 | 0.003966 | 0.003966 | — | 0.003787 | 0.0008537 |
| 40 | -0.075 | 1.000 | 0.000 | 0.000 | 0.003742 | 0.003742 | — | 0.003737 | 0.000836 |
| 40 | -0.050 | 0.999 | 0.001 | 0.000 | 0.1662 | 0.004202 | — | 0.003688 | 0.0008187 |
| 40 | -0.025 | 0.939 | 0.061 | 0.000 | 9.93 | 0.003559 | 0.001197 | 0.00364 | 0.0008017 |
| 40 | +0.000 | 0.508 | 0.492 | 0.000 | 82.74 | 0.003742 | 0.0007576 | 0.003592 | 0.0007851 |
| 40 | +0.025 | 0.066 | 0.934 | 0.000 | 10.75 | 0.004481 | 0.0007245 | 0.003545 | 0.0007688 |
| 40 | +0.050 | 0.001 | 0.999 | 0.000 | 0.1633 | — | 0.0007266 | 0.003499 | 0.0007529 |
| 40 | +0.075 | 0.000 | 1.000 | 0.000 | 0.0006914 | — | 0.0006914 | 0.003454 | 0.0007374 |
| 40 | +0.100 | 0.000 | 1.000 | 0.000 | 0.0007807 | — | 0.0007807 | 0.003409 | 0.0007223 |
| 40 | +0.150 | 0.000 | 1.000 | 0.000 | 0.0006277 | — | 0.0006277 | 0.003323 | 0.0006931 |
| 40 | +0.200 | 0.000 | 1.000 | 0.000 | 0.0006725 | — | 0.0006725 | 0.003239 | 0.0006654 |

## 6C. Local theory vs Monte Carlo (pre-registered hypotheses)

Local theory = frequency block of the misspecified sandwich covariance (σ²/2)H⁻¹Re(JᴴJ)H⁻¹ at each noiseless branch root (H includes residual curvature; FD-validated to 4.2e-11). It is called a *local Hessian (sandwich) covariance approximation*, not a CRB/MCRB: the MCRB's regularity assumption (unique interior pseudo-true parameter) fails at η = 0 by construction.

* **H1 (branch-conditioned agreement away from the exchange; |η| ≥ 0.10, 30/40 dB, ≥ 100 selections): True** — trace ratios A@30dB,η=-0.20: 0.980, A@30dB,η=-0.15: 1.016, A@30dB,η=-0.10: 1.085, B@30dB,η=+0.10: 1.030, B@30dB,η=+0.15: 0.998, B@30dB,η=+0.20: 0.986, A@40dB,η=-0.20: 0.990, A@40dB,η=-0.15: 0.975, A@40dB,η=-0.10: 1.048, B@40dB,η=+0.10: 1.081, B@40dB,η=+0.15: 0.906, B@40dB,η=+0.20: 1.012.
* **H2 (global failure at the exchange; MSE to A ≥ 10 × local trace of A at η = 0): True** — ratios 20 dB: 221×, 30 dB: 2217×, 40 dB: 22287×.
* **H3 (bimodality at η = 0; min(P_A, P_B) ≥ 0.2): True** — 20 dB: P_A=0.464, P_B=0.489, 30 dB: P_A=0.510, P_B=0.490, 40 dB: P_A=0.508, P_B=0.492.
* Pre-registered 'meaningful TSP consequence' (H1 ∧ H2): **True**.

## Transition widths and comparison with the archived experiment

| SNR | 10–90% probit width, pre-registered labels (OUT counts as not-A) | POST_HOC: nearest-root labels | archived (manuscript) |
|---|---|---|---|
| 20 dB | 0.4790 | 0.4181 | 0.385 [0.366, 0.404] |
| 30 dB | 0.1295 | 0.1295 | 0.129 [0.123, 0.136] |
| 40 dB | 0.0419 | 0.0419 | 0.039 [0.0368, 0.0419] |

POST_HOC_SENSITIVITY: all 580 OUT records occur at 20 dB; 99.8% are nearest to A at distance 1.00–2.65 (median 1.21) — displaced A-family fits, not a third representation. The archived label rule (lower of two seeded refinements) assigns these to A.

Findings: (i) the 30 and 40 dB widths reproduce the archived values; (ii) the 20 dB width is **not robust**: 0.479 (pre-registered labels), 0.418 (nearest-root labels) vs archived 0.385 with a model-conditional CI of ±0.02 that excludes both — the archived interval understates label-rule and η-design sensitivity at 20 dB; (iii) at 20 dB the local sandwich trace of A over-predicts the A-conditioned MSE (0.23–0.25 vs 0.33–0.40) because conditioning on the 1.0-radius label truncates the A-family tail — local theory is only reliable here at ≥ 30 dB.

Figures: `stage6_noise/selection_probability.png|pdf`, `stage6_noise/bimodality.png`, `stage6_noise/global_vs_conditional_error.png|pdf`, `stage6_noise/local_theory_comparison.png|pdf`. Trial-level data: `stage6_noise/trials.csv` (39,000 rows).