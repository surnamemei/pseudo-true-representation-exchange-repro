# Derivation and deterministic predictions

Sealed with `SEAL_2` **before any Monte Carlo outcome was read**. The plan is `00_THEORY_PLAN.md` (`SEAL_1`, digest `6051209b333985ac5c79d7e6699703b3db4c669f56caaddc8350de9d88e7783b`).

The code is `code/transition_theory.py` and `code/second_order_diagnostic.py`. The outputs are in `results/transition_theory/deterministic/`: `predictions.csv`, `deterministic.json`, `checks.json` and `second_order_diagnostic.json`.

## 1. Derivation

**Model.**
- The parameters are θ = (Re c₁, Im c₁, Re c₂, Im c₂, u₁, u₂), and the model is m(θ) = V(u)c with V(u) = [e^{i u₁ s}, e^{i u₂ s}].
- The branch-j cost of a record y is J_j(y) = min ‖y − m(θ)‖², taken over the basin of branch j.
- At the noiseless record x(ε):
  - θ_j* is the minimizer and r_j = x − m(θ_j*) the residual;
  - M_j = ∂m/∂θ at θ_j*, and stationarity gives Re(M_jᴴ r_j) = 0.

**(a) Slope in ε.** The record is x(ε) = x₀ + ε a, with a = e^{iπ}e^{i b s}. By the envelope theorem, the minimizer's motion does not enter the first derivative:

dJ_j/dε = 2 Re⟨r_j, a⟩, where ⟨p, q⟩ = pᴴq.

Hence

g = dΔ₀/dε at ε_c = 2 Re⟨r_A − r_B, a⟩.

This is the finite-spacing slope at N = 21, z = 2. It is not the small-spacing tangent Γ.

**(b) Perturbation by noise.** Write y = x + w and θ = θ_j* + δ. Then

‖r_j + w − M_jδ − ½ δᵀ(∂²m)δ‖² = ‖r_j + w‖² − 2δᵀ Re(M_jᴴ w) + δᵀ H_j δ + O(‖w‖‖δ‖², ‖δ‖³),

where H_j = Re(M_jᴴ M_j) − Re Σ_t r̄_{jt} ∂²m_t. This is half the Hessian of the loss, as in `s6_noise.local_cov`.

Minimizing over δ gives δ* = H_j⁻¹ b_j, with b_j = Re(M_jᴴ w). Therefore

J_j(x + w) = J_j(x) + 2 Re⟨r_j, w⟩ + ‖w‖² − q_j(w) + O(‖w‖³), with q_j(w) = b_jᵀ H_j⁻¹ b_j ≥ 0.

The ‖w‖² terms cancel between the branches. With d = r_A − r_B,

Δ_w = Δ₀ + 2 Re⟨d, w⟩ − [q_A(w) − q_B(w)] + O(‖w‖³).

The first-order theory keeps Δ₀ + 2 Re⟨d, w⟩.

**(c) Distribution of the linear term.** The noise is circular (E wwᴴ = σ²I, E wwᵀ = 0). So dᴴw is circular complex Gaussian with variance σ²‖d‖², and its real part has half of that:

Var[2 Re⟨d, w⟩] = 4 · σ²‖d‖²/2 = 2σ²‖d‖², and s_Δ = √2 σ ‖d‖.

**(d) Selection.** A is selected when Δ_w < 0. This assumes the noisy global minimum is one of the two continued branches. Then

P_A = Φ(−Δ₀/s_Δ), and P_B = 1 − P_A.

This is version (B): exact Δ₀(η), d(η) and σ(η) at each ε = ε_c(1 + η), taken from the frozen continuation.

**(e) Local law.** Since Δ₀(η) = g ε_c η + O(η²),

P_A ≈ Φ(−Kη), with K = g ε_c / (√2 σ_c ‖d_c‖).

This gives a 10–90% width W = 2Φ⁻¹(0.9)/|K| ∝ σ_c ∝ 10^{−SNR/20}. The slope is positive, g > 0, so A is selected below the crossing and B above: the certified orientation.

**(f) Global MSE.** Suppose that, conditional on the selected branch j, the fitted pair is approximately N(u_j, C_j), with C_j the local sandwich covariance, and that this is independent of the selection event. Then E‖û − u_W‖² = Σ_j P_j [tr C_j + ‖u_j − u_W‖²]. This splits into:
- the within-mode term Σ P_j tr C_j;
- the between-mode term Σ P_j ‖u_j − u_W‖².

The identity Σ_j P_j ‖u_j − u_W‖² = P_A P_B ‖u_A − u_B‖² + ‖Σ_j P_j u_j − u_W‖² links this to the total-variance form, which measures error about the mixture mean.

**(g) Neglected terms.**
- **Re-optimization bias.** E[q_j] = (σ²/2) tr(H_j⁻¹ G_j), with G_j = Re(M_jᴴ M_j). The mean gap is therefore Δ₀ − E[q_A − q_B], which shifts the effective crossing by δη = E[q_A − q_B]/(g ε_c). This is computed below as a diagnostic only.
- **Higher-order spread.** The variance of q_A − q_B is O(σ⁴).
- **Selection correlation.** The noise that decides the selection also moves the within-branch error. It is correlated with it, which is not modelled.
- **Escape solutions.** OUT solutions are not modelled.

## 2. Deterministic results

These come from float64 computations with the frozen code; the crossing quantities were also re-derived in 50-digit mpmath.

**Crossing point.**

| Quantity | Value |
|---|---|
| ε_c (Stage-6 float) | 0.2481906301722774 |
| 50-digit root of Δ₀ = 0 | 0.24819063017227739820…; η = −4.7 × 10⁻¹⁷ |
| J_A = J_B at ε_c | 0.8118078719175 (Δ₀ = 1.2 × 10⁻¹⁴) |
| u_A, u_B at ε_c | (−4.772447, 0.605170), (−0.057654, 12.463419) |
| ‖u_A − u_B‖² | 162.8473507 |
| ‖x(ε_c)‖² | 64.894024695 |

**Slope and residual gap.**

| Quantity | Value |
|---|---|
| **g = dΔ₀/dε at ε_c** | **8.259051045764851**; the branch slopes are dJ_A/dε = 7.177553 and dJ_B/dε = −1.081498 |
| g, 50-digit | 8.2590510457648658139644744… |
| **‖r_A − r_B‖ at ε_c** | **1.3073721889711298**; 50-digit 1.3073721889711320440… |
| ‖r_A − r_B‖ over the 13 η | 1.2119 at η = −0.20, rising monotonically to 1.4313 at η = +0.20 |

**Per SNR, at the crossing.**

| SNR | σ_c | s_Δ at crossing | **K** | **W₁₀₋₉₀, closed form** | W₁₀₋₉₀, exact curve |
|---|---:|---:|---:|---:|---:|
| 20 dB | 0.175789 | 0.325018 | **6.30679** | **0.406404** | 0.409847 |
| 30 dB | 0.055589 | 0.102780 | **19.94383** | **0.128516** | 0.128623 |
| 40 dB | 0.017579 | 0.032502 | **63.06793** | **0.040640** | 0.040644 |

**Width ratios.**
- Predicted W₃₀/W₄₀ = √10 = 3.162278 in the closed form, and 3.164644 on the exact curve.
- Predicted W₂₀/W₃₀ = 3.162278.

**Neglected second-order mean term** (diagnostic only; not used):
- tr(H⁻¹G) is 5.6041 for A and 5.7531 for B. A correctly specified fit would give 6.
- The implied crossing shifts are δη = −1.1 × 10⁻³, −1.1 × 10⁻⁴ and −1.1 × 10⁻⁵ at 20, 30 and 40 dB, i.e. −0.28%, −0.09% and −0.03% of W.
- So the omitted mean term is negligible at every SNR. The leading unmodelled risks are OUT solutions (20 dB) and the selection–error correlation (MSE).

The per-η predictions, P_A, P_A,local, tr C_A, tr C_B, and the within-mode, between-mode and total MSE, are in `predictions.csv`.

**Predicted P_A at 30 and 40 dB** (version B):

| η | −0.20 | −0.15 | −0.10 | −0.075 | −0.05 | −0.025 | 0 | 0.025 | 0.05 | 0.075 | 0.10 | 0.15 | 0.20 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 30 dB | 1.0000 | 0.9986 | 0.9771 | 0.9329 | 0.8410 | 0.6911 | 0.5000 | 0.3092 | 0.1599 | 0.0681 | 0.0236 | 0.0015 | 0.0000 |
| 40 dB | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 0.9992 | 0.9427 | 0.5000 | 0.0576 | 0.0008 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |

**Predicted global MSE** (u units, winner-referenced; B at η = 0):

| η | −0.20 | −0.10 | −0.05 | −0.025 | 0 | 0.025 | 0.05 | 0.10 | 0.20 |
|---|---|---|---|---|---|---|---|---|---|
| 30 dB | 0.0456 | 3.78 | 26.0 | 50.4 | 81.4 | 50.3 | 26.0 | 3.85 | 0.0134 |
| 40 dB | 0.00399 | 0.00379 | 0.134 | 9.34 | 81.4 | 9.38 | 0.135 | 0.00072 | 0.00067 |

Wherever 0.01 < P_A < 0.99, the between-mode term is at least 97.9% of the predicted MSE, at every SNR. The within-mode term ranges from 6.7 × 10⁻⁴ to 0.37 and scales as σ². It dominates only in the tails.

At η = 0 the A-referenced and B-referenced predictions are equal, because P_A = 1/2 there: 81.643, 81.446 and 81.426 at 20, 30 and 40 dB.

## 3. Implementation checks (all pass; `checks.json`)

| Check | Result | Tolerance |
|---|---|---|
| g: envelope vs Richardson finite difference, float64 | relative difference 1.2 × 10⁻¹⁰ | ≤ 10⁻⁶ |
| g: envelope vs finite difference, 50 digits (h = 10⁻¹² ε_c) | relative difference 8.3 × 10⁻²⁶ | ≤ 10⁻²⁰ |
| g: float64 vs mpmath | 1.7 × 10⁻¹⁵ | — |
| ‖r_A − r_B‖: float64 vs mpmath | 1.7 × 10⁻¹⁵ | ≤ 10⁻⁹ |
| Spread formula: synthetic Gaussian, fixed d_c, M = 2 × 10⁶ | variance ratio − 1 = 4.8 × 10⁻⁴; mean/s_Δ = 2.9 × 10⁻⁴ | ≤ 5 × 10⁻³; ≤ 3.5 × 10⁻³ |
| Orientation | g > 0; P_A > 0.5 below and < 0.5 above at every SNR, in both versions; monotone | — |
| Frozen continuation reproduced | J relative 4.9 × 10⁻¹⁴; u 1.4 × 10⁻¹³ | ≤ 10⁻⁹; ≤ 10⁻⁷ |
| Fitted parameters | none | — |

The mpmath roots were Newton-polished to gradient norms of 1.5 × 10⁻⁵⁰ and 3.1 × 10⁻⁵¹.

**Validation script.** `code/validate.py` was run only on synthetic stand-in files with the frozen schema. On a mock with forced failures, it returned FAILS_PREDICTIVELY with F3 at 30 dB and F2 at 40 dB, which confirms the classification logic. It is sealed unchanged before the one-shot run.
