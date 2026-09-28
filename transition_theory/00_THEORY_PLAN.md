# Predictive theory of branch selection and global MSE near the certified exchange: plan

Written and sealed on 2026-09-28, **before any Monte Carlo outcome was read**. The seal record is in `SEAL_1.json` and `SEAL_1.sha256`, next to this file. The machine-readable form is `transition_spec.json`; where the two differ, this plan governs.

This is one theory with one frozen validation. There are no new Monte Carlo records, no new SNR, η, N, z, signal family or model order, no retuned labels, and no fitted transition parameter. No manuscript file changes until the author approves the result.

## 0. Blinding record

**Read before this plan was sealed** (all method-only or deterministic):
- **Code:**
  - `validation/adversarial_overnight/s6_noise.py`, the Stage-6 driver: noise, labels, sandwich covariance, the analysis and the probit width;
  - `advcore.py`, the model and VarPro objective;
  - `advrun.py`, its header only;
  - `tools/posthoc_stage6.py`.
- **Plans:** `PREREGISTRATION.md` §Stage 6 and `FREEZE.md` §7–10.
- **Deterministic data:** `validation/adversarial_overnight/stage6_noise/noiseless_roots.json`, the noiseless A/B continuation. At η = 0, J_A − J_B = 1.25 × 10⁻¹⁴.
- **Conventions from a blinded subagent.** A subagent read `supplement/S7_noise_table.tex`, `supplement/S7_noise_methods.tex` and the noise section of `paper/sections/body.tex`. It reported conventions only, with every outcome number replaced by "#".
- **File names only.** A name-only search located the files containing the two widths quoted in the task.
- **Fingerprints.** The SHA-256 of each outcome file was computed without displaying its content.

**Not read before seal 2:**
- `stage6_noise/summary.csv`, `trials.csv`, `blocks/`, `DONE.json`, `POSTHOC_labels_and_widths.json` and `STAGE6_REPORT.md`, plus the Stage-6 figures;
- `S7_noise_table.tex`, `S7_noise_methods.tex` and the noise section of `body.tex`, directly;
- `noise_branch_trials.csv` and the other archived noise outputs;
- `MORNING_REPORT.md`, `FINAL_SUMMARY.txt`, `CLAIM_AUDIT.md` and `HOSTILE_REVIEW.md`.

**Outcome information known before sealing.** The task statement quotes the preregistered probit widths, about 0.1295 at 30 dB and 0.0419 at 40 dB, a ratio of about 3.09. **The width-scaling check (F1 below) is therefore not blind.** The absolute predicted widths, every P(A) curve and every MSE prediction are computed without any other outcome.

**Found during reconnaissance: overlap with the manuscript.** The noise section of `body.tex` already contains two subsections, whose numbers were not read:
- a "frozen-branch cost-gap approximation", Pr(A) ≈ Φ(−ΔJ_det/σ_Δ), with predicted widths, evaluated against the older archived 400-record experiment;
- an exact two-mode error decomposition.

The present analysis is therefore:
- a blind, preregistered test of that kind of first-order theory on the preregistered 1000-record Stage-6 data;
- plus the closed-form local width law with the exact finite-spacing slope;
- plus the winner-referenced two-branch MSE prediction.

This overlap matters for Section 9, on what may enter the paper.

## 1. Setting (frozen conventions, from code)

- **Record.** x_t(ε) = e^{−i z s_t} + e^{+i z s_t} + ε e^{iπ} e^{i b s_t}, with s_t = (t − (N−1)/2)/N for t = 0, …, N−1, and N = 21, z = 2, b = 10 (`advcore.make_signal`).
  - Frequencies u are in normalized units: a tone is e^{i u s_t}.
  - The strong pair is at u = ±2 and the weak tone at u = 10.
- **Fit.** J(u; y) = min over c ∈ ℂ² of ‖y − V(u)c‖², the unweighted residual sum of squares (`advcore.Objective`). Here V(u) = [e^{i u₁ s}, e^{i u₂ s}].
- **Grid.** ε = ε_c(1 + η) with ε_c = 0.2481906301722774, the Stage-6 float of the certified N = 21, z = 2 crossing. The 13 η values are
  η ∈ {−0.20, −0.15, −0.10, −0.075, −0.05, −0.025, 0, 0.025, 0.05, 0.075, 0.10, 0.15, 0.20}.
- **Branches.** A and B are the noiseless local minimizers continued from the frozen crossing roots (`s6_noise.noiseless_roots`). The certified orientation is: A is globally optimal below ε_c, and B above.
- **Noise.** w_t are i.i.d. circular complex Gaussian, with Re w_t and Im w_t each N(0, σ²/2). The variance is σ² = (‖x(ε)‖²/N)·10^{−SNR/10}, for SNR ∈ {20, 30, 40} dB.
- **Labels (Monte Carlo).** A record is labelled A if its fitted pair lies within max-norm distance < 1 of the noiseless A root (unordered pairs, on the torus); B likewise; otherwise OUT.
  - Empirical **P(A) = n_A/n**, with OUT counted as not A and n = 1000.
- **Global MSE (frozen S7 table and Stage-6 `MSE_to_winner`).** The mean over all n records, OUT included, of ‖û − u_W‖². Pairs are sorted, and the distance is squared Euclidean in u.
  - u_W = u_A for η < 0, and u_W = u_B for η ≥ 0: the tie at η = 0 is referenced to B.
  - The table caption also reports A-referenced values at η = 0 (`MSE_to_A`).
- **Local covariance C_j.** The frequency block of the frozen sandwich covariance `s6_noise.local_cov`, (σ²/2)·H⁻¹Re(JᴴJ)H⁻¹, with 6 real parameters. H includes the residual curvature. It is evaluated at each branch root with the σ² above.

## 2. First-order branch-selection theory

**Branch costs.** J_j(y) is the branch-j locally minimized cost, for j ∈ {A, B}. For the noiseless record, θ_j* is the minimizer (two frequencies and two complex amplitudes), and r_j = x − V(u_j*)c_j* is the residual.

**Noiseless gap.** Δ₀(η) = J_A(x) − J_B(x).

**Perturbation.** Evaluating the noisy cost at the noiseless minimizer, and using stationarity (the envelope theorem):

J_j(x + w) = J_j(x) + 2 Re⟨r_j, w⟩ + ‖w‖² − q_j(w) + o(‖w‖²),

where q_j(w) ≥ 0 is the quadratic reduction from re-optimizing the branch.

**Noisy gap.** The ‖w‖² terms cancel in the difference. To first order,

Δ_w = Δ₀ + 2 Re⟨d, w⟩, with d = r_A − r_B.

**Spread.** Under the frozen noise model, Re⟨d, w⟩ ~ N(0, σ²‖d‖²/2). So

Var[2 Re⟨d, w⟩] = 2σ²‖d‖², and s_Δ = √2 σ ‖d‖.

**Selection probability.** The global fit selects A when Δ_w < 0, assuming the noisy global minimum is one of the two continued branches (P_B = 1 − P_A). Then:

- **(B, finite-range and primary):** P_A(η) = Φ(−Δ₀(η) / s_Δ(η)), with Δ₀(η), d(η) and σ(η) evaluated exactly at each ε = ε_c(1 + η).
- **(A, local closed form):** P_A,local(η) = Φ(−K η), where
  - K = g ε_c / (√2 σ_c ‖d_c‖);
  - g = dΔ₀/dε at ε_c, the finite-spacing envelope slope;
  - σ_c and d_c are evaluated at ε_c.

**The slope g.** By the envelope theorem, dJ_j/dε = 2 Re⟨r_j, ∂x/∂ε⟩ with ∂x/∂ε = e^{iπ}e^{i b s}. Hence

g = 2 Re⟨r_A − r_B, e^{iπ}e^{i b s}⟩ at ε_c.

This is not the small-spacing tangent Γ.

**Orientation.** A below and B above requires g > 0. P_A must then decrease through 0.5 as η increases through 0, in both versions.

**Cost scale.** The cost is the unnormalized residual sum of squares. K, W and P_A do not depend on any rescaling of J, because g and s_Δ scale together.

## 3. Predicted transition width

The primary prediction is the closed form:

W₁₀₋₉₀ = 2Φ⁻¹(0.9)/|K| = 2Φ⁻¹(0.9) · √2 σ_c ‖d_c‖ / (|g| ε_c), with Φ⁻¹(0.9) = 1.2815515655446004.

**Scaling with SNR.** σ_c² = (‖x(ε_c)‖²/N)·10^{−SNR/10}, so W ∝ σ_c ∝ 10^{−SNR/20}.
- A +10 dB change multiplies W by exactly 10^{−1/2} in the closed form, because g, ‖d_c‖ and ε_c do not depend on SNR.
- The predicted 30→40 dB ratio is **W₃₀/W₄₀ = √10 = 3.1623**. The exact-curve version adds only the weak ε-dependence of x and d.

**Secondary: the exact-curve width.** It is the difference between the η values at which Φ(−Δ₀(η)/s_Δ(η)) equals 0.1 and 0.9. They are found by root-finding on the deterministic continuation. No Monte Carlo is involved, and the Stage-6 grid itself is unchanged.

**Predicted midpoint.** η₅₀,pred is the root of Δ₀(η) = 0, which is 0 to float precision, since ε_c is the crossing. The local form has η₅₀ = 0 exactly.

## 4. Two-branch global-MSE theory

**Prediction.** With P_B = 1 − P_A and P_A from the primary version (B):

MSE_pred(η, SNR) = P_A [tr C_A + ‖u_A − u_W‖²] + P_B [tr C_B + ‖u_B − u_W‖²].

**Decomposition.**
- within-mode = P_A tr C_A + P_B tr C_B;
- between-mode = P_A ‖u_A − u_W‖² + P_B ‖u_B − u_W‖².

**The reference u_W.** It follows the frozen table: A for η < 0 and B for η ≥ 0.
- At η = 0, u_W = u_B is the primary quantity, compared with `MSE_to_winner`.
- The A-referenced prediction at η = 0 is also reported, as secondary and descriptive, against `MSE_to_A`, the values in the table caption.

**Not used.** The P_A P_B ‖u_A − u_B‖² form measures error about the mixture mean. It is not the manuscript's winner-referenced MSE. The identity E‖û − u₀‖² = Σ p_j tr Cov_j + p_A p_B ‖μ_A − μ_B‖² + ‖Eû − u₀‖² shows the two are consistent once the squared-bias term is included.

**Descriptive only.** MSE_pred is also reported with P_A,local.

## 5. Second-order caution (status fixed now)

The first-order theory neglects:
- **Re-optimization terms.** The O(σ²) terms q_A − q_B shift the effective crossing by δη ≈ E[q_A − q_B]/(g ε_c).
- **Selection correlation.** The noise that decides the selection is correlated with the within-branch frequency error, which shifts conditional means.
- **Escape solutions.** OUT fits (displaced A-family, coalescent or other solutions) are ignored; the theory assumes P_B = 1 − P_A.
- **Near the exchange.** The local-covariance accuracy is not guaranteed close to the exchange.

**Regimes.** **30 and 40 dB are the primary validation regime. 20 dB is a stress test only.** It is reported but never used in the classification. This status does not change after the results are seen.

## 6. Deterministic implementation checks (before any outcome is read)

1. **Slope g, two independent methods.**
   - (i) The envelope formula.
   - (ii) A centered finite difference of the re-optimized branch costs:
     - in float64, steps h = 10⁻⁴ ε_c and h/2 with Richardson extrapolation;
     - at the crossing, in 50-digit mpmath with Newton-polished roots, a single centered difference with h = 10⁻¹² ε_c. Its truncation error is O(h²) and its rounding error O(10⁻⁵⁰/h), both far below the tolerance.

   Pass: relative agreement ≤ 10⁻⁶ in float64 and ≤ 10⁻²⁰ in mpmath.
2. **‖r_A − r_B‖** at the crossing and at all 13 η, in float64. The crossing value is also computed in mpmath; agreement ≤ 10⁻⁹ relative.
3. **Spread formula.** A synthetic Gaussian check on the fixed crossing vector d_c, with no optimizer and no branch outcome.
   - M = 2,000,000 draws with the frozen noise convention; the seed is `SeedSequence([20260928, 99, 1])`, distinct from every Stage-6 seed.
   - Pass: |sample variance / 2σ²‖d‖² − 1| ≤ 5·√(2/M), and |sample mean|/s_Δ ≤ 5/√M.
4. **Orientation.** At every SNR, in both versions, predicted P_A > 0.5 at η < 0 and < 0.5 at η > 0.
5. **Reproduced roots.** The float64 continuation is recomputed. J must agree with `noiseless_roots.json` to ≤ 10⁻⁹ relative, and u to ≤ 10⁻⁷.
6. **No fitted parameter.** Every input is deterministic and listed in `transition_spec.json`.

**Seal 2.** Before any outcome file is opened, the following are hashed into `SEAL_2.json` and made read-only:
- this plan's seal;
- `DERIVATION.md`;
- the implementation, including the validation script;
- the deterministic inputs;
- the predicted curves;
- the check results;
- the fingerprints of the outcome files.

The project is not a git repository, so the seal is a SHA-256 record in the project's freeze style, not a git commit.

## 7. One-shot validation on the frozen Stage-6 Monte Carlo

**Inputs,** read only after seal 2:
- `stage6_noise/summary.csv`: n, P_A, P_B, P_OUT, `MSE_to_winner`, `MSE_to_A`, `MSE_to_B`, `local_trace_A` and `local_trace_B`;
- `stage6_noise/DONE.json`, the `probit_widths` produced by the preregistered analysis: `width_10_90`, `slope`, `intercept`, `eta50`.

Both are hash-checked against the fingerprints in seal 2, which equal the reviewer-archive manifest. `trials.csv` is not needed. No new noisy record is generated, and nothing is re-fitted.

**Reported for each SNR (20, 30, 40):**
- **Selection:**
  - predicted P(A) and P(A)_local against empirical P(A) at all 13 η;
  - maximum and mean absolute error for each version;
  - the Brier score of the theoretical probabilities, over records, from the counts;
  - for reference, the Brier score of the empirical frequencies, which is the in-sample floor;
  - W_pred, closed form and exact curve, against W_emp, the preregistered probit width;
  - the ratio W_pred/W_emp;
  - η₅₀,emp against η₅₀,pred;
  - the 30→40 dB width ratios, predicted and empirical;
  - P_OUT, which the theory treats as zero.
- **MSE:** MSE_pred, MSE_emp (`MSE_to_winner`), the within- and between-mode terms, emp/pred and log₁₀(emp/pred), and at η = 0 also the A-referenced pair.
- **Consistency check.** The recomputed tr C_A and tr C_B against the frozen `local_trace_A` and `local_trace_B`.

## 8. Pre-declared interpretation (fixed now; not revised after results)

**Quantities used.**
- **Widths:** the closed-form W_pred (Section 3).
- **P(A) errors and MSE:** version (B) with the frozen table reference.
- **Empirical width and midpoint:** the preregistered probit in `DONE.json`.
- **Central region:** the 9 points with |η| ≤ 0.1.

**Criteria, evaluated separately at 30 dB and at 40 dB.**
- **C1, orientation and crossing.**
  - Empirical P(A) at η = −0.20 is > 0.5, and at η = +0.20 is < 0.5.
  - The preregistered probit slope is < 0.
  - |η₅₀,emp − η₅₀,pred| ≤ W_pred/2.
- **C2, width:** |W_pred/W_emp − 1| ≤ 0.25.
- **C3, selection error:** the mean over the 13 η of |P_pred(A) − P_emp(A)| is ≤ 0.08.
- **C4, MSE order of magnitude everywhere:** 0.1 ≤ MSE_emp/MSE_pred ≤ 10 at all 13 η.
- **C5, central MSE accuracy:** 0.5 ≤ MSE_emp/MSE_pred ≤ 2.0 at a strict majority, at least 5, of the 9 central points.

**Failure conditions.**
- **F1, width scaling wrong:** |(W₃₀/W₄₀)_emp / (W₃₀/W₄₀)_pred − 1| > 0.25. This check is not blind; see Section 0.
- **F2, orientation or crossing wrong, at 30 or at 40 dB:** the empirical P(A) orientation test fails, or the probit slope is ≥ 0, or |η₅₀,emp − η₅₀,pred| > W_pred.
- **F3, MSE missed, at 30 or at 40 dB:** MSE_emp/MSE_pred lies outside [0.1, 10] at a strict majority, at least 5, of the 9 central points.

**Outcome.**
- **STRONG_PREDICTION:** C1–C5 all hold at both 30 and 40 dB.
- **FAILS_PREDICTIVELY:** F1, or F2, or F3.
- **USEFUL_FIRST_ORDER:** otherwise. That is, width scaling and the qualitative selection curve are right, the MSE prediction is within an order of magnitude over most of the centre, and at least one quantitative criterion of STRONG fails.

If version (B) succeeds but the local g·η law degrades away from η = 0, that is expected. It is reported as "local asymptotic law + finite-range deterministic refinement".

## 9. What may enter the paper (the author decides; no manuscript edit before approval)

The recommendation follows the task's rules:
- **STRONG:** a concise proposition or derivation, with the predicted and empirical widths and the two-branch MSE, and a reframed Section VIII.
- **USEFUL:** only a compact interpretation: W ∝ σ, and that the 30→40 dB narrowing agrees with first-order theory. It is not called a fully predictive model.
- **FAILS:** archive the analysis and leave the paper unchanged.

In every case, the report must say how the result relates to the manuscript's existing frozen-branch cost-gap subsection (Section 0).

## 10. Not done, regardless of outcome

No:
- extra SNR values, Monte Carlo records, η or ε points, prevalence scans or application examples;
- new signal family, N, z or model order;
- change to labels or the one-unit radius;
- second predictive model.

## 11. Outputs

- **`paper/transition_theory/`:** `00_THEORY_PLAN.md`, `transition_spec.json`, the `SEAL_1` and `SEAL_2` records, `DERIVATION.md`, `VALIDATION_REPORT.md`, and `code/`, which holds `transition_theory.py` and `validate.py`.
- **`results/transition_theory/`:** `deterministic/`, with the inputs, the predicted curves and the checks, sealed at seal 2; and `validation/`, with the comparison tables and figures.
