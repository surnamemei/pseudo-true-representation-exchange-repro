# Validation report: predictive theory of branch selection and global MSE near the certified exchange

**Outcome: STRONG_PREDICTION.** The criteria were fixed in `00_THEORY_PLAN.md` §8 and applied once by the sealed `code/validate.py` on 2026-09-28 at 18:19:55 +10:00.

No manuscript file has been changed. Any change awaits the author's approval.

## 0. Provenance

**Order of events.**

| Step | Time (+10:00) | Content |
|---|---|---|
| SEAL_1 | 18:12:43 | Plan and spec. Digest `6051209b333985ac5c79d7e6699703b3db4c669f56caaddc8350de9d88e7783b`. |
| SEAL_2 | about 18:19 | 16 files: the derivation, the three scripts, the deterministic inputs and outputs, and the checks. Chained to SEAL_1. Digest `11ff99d6615eeb89e6957aa08ffc038ba573054308bd53964e928740d3b67f47`. |
| Validation | 18:19:55 | One run of `code/validate.py`. The script refuses a second run. |

The project is not a git repository. The seals are SHA-256 records with read-only permissions, and both digests were stated in the working session before the validation ran.

**Monte Carlo inputs.** Two files, both hash-identical to their SEAL_2 fingerprints and to the reviewer-archive manifest:
- `validation/adversarial_overnight/stage6_noise/summary.csv`;
- `validation/adversarial_overnight/stage6_noise/DONE.json`.

No other outcome file was needed or read. No noisy record was generated, and nothing was re-fitted.

**Blinding caveat.** The task statement quoted the empirical 30 and 40 dB widths, so the width-scaling check (F1) is not blind. Nothing in the theory was adjustable: the widths follow from g, ‖r_A − r_B‖, ε_c and σ.

**Correction to the sealed plan, §0.** The plan says the manuscript's noise section already contains a "frozen-branch cost-gap approximation". That is true of `paper/sections/body.tex` and `body_10page.tex`, but **the manuscript does not input those files.** The active files are:
- `main_tsp_reviewfriendly.tex` and `tsp_submission_main.tex`, which input `body_reviewfriendly.tex`;
- `supplement/supplement.tex`, which inputs `S7_noise_methods.tex`.

Neither the active main text nor the active supplement contains a selection-probability theory. The frozen-projector approximation, with its predicted widths 0.4101, 0.1286 and 0.04064, survives only in inactive files: an earlier draft, and the unused `supplement/S5_noise_search_details.tex` and `S7_noise_uncertainty.tex`. There it was compared, not blind, with the archived 400-record experiment.

The present result is therefore new relative to the manuscript as it stands, and consistent with that earlier draft. The predictions, the criteria and the outcome are unaffected by this correction.

## 1. Answers

| # | Item | Result |
|---|---|---|
| 1 | Exact finite-spacing slope g = dΔ₀/dε at ε_c, with Δ₀ = J_A − J_B the unnormalized residual sums of squares | **8.259051045764851**. Envelope formula; the float64 Richardson finite difference agrees to 1.2 × 10⁻¹⁰, and 50-digit mpmath gives 8.25905104576486581… (finite difference agrees to 8 × 10⁻²⁶). Branch slopes: dJ_A/dε = 7.1776, dJ_B/dε = −1.0815. |
| 2 | ‖r_A − r_B‖ at the crossing | **1.3073721889711298** (mpmath agrees to 1.7 × 10⁻¹⁵). Over the 13 η it rises from 1.2119 at η = −0.2 to 1.4313 at η = +0.2. |
| 3 | K = g ε_c / (√2 σ_c ‖r_A − r_B‖) | 20 dB: **6.3068**; 30 dB: **19.9438**; 40 dB: **63.0679** |
| 4 | W₁₀₋₉₀ = 2Φ⁻¹(0.9)/K = 2.3119 σ_c = 4.0640 × 10^{−SNR/20} | 20 dB: **0.4064**; 30 dB: **0.12852**; 40 dB: **0.040640**. The exact curve gives 0.4098, 0.12862 and 0.040644. |
| 5 | Empirical preregistered probit widths | 20 dB: 0.47897; 30 dB: **0.12950**; 40 dB: **0.04190** |
| 6 | P(A), predicted vs empirical, mean and max absolute error over 13 η | 30 dB: **0.0047 / 0.0110**; 40 dB: **0.0016 / 0.0084**; 20 dB (stress): 0.045 / 0.106 |
| 7 | Global-MSE prediction, empirical/predicted | 30 dB: 0.85–1.35 for \|η\| ≤ 0.1 and 0.15–1.92 over all 13. 40 dB: 0.74–1.24 over all 13. 20 dB (stress): 0.94–1.10 over all 13. |
| 8 | Within-mode vs between-mode | Where 0.01 < P(A) < 0.99 in the central region, the between-mode term is ≥ 99.0% (30 dB) and ≥ 99.96% (40 dB) of the predicted MSE. The within-mode term, which scales as σ² (0.022 at 30 dB and 0.0022 at 40 dB, at η = 0), dominates only in the tails. |
| 9 | Outcome | **STRONG_PREDICTION** |
| 10 | Strong enough to change Section VIII? | Yes, by the pre-declared rule. A focused rewrite of the section is justified; see §7. |
| 11 | Suggested manuscript changes | §8, not applied |
| 12 | No new Monte Carlo, no scientific tuning | Confirmed; see §10 |

## 2. Branch selection

| SNR | W_pred, closed form | W_pred, exact curve | W_emp, preregistered probit | Pred/emp | η₅₀,emp (predicted 0) | \|Δη₅₀\|/W_pred | Brier: theory / empirical-frequency floor |
|---|---:|---:|---:|---:|---:|---:|---|
| 20 dB (stress) | 0.4064 | 0.4098 | 0.4790 | 0.848 | −0.0244 | 0.060 | 0.20655 / 0.20366 |
| **30 dB** | **0.12852** | 0.12862 | **0.12950** | **0.992** | 0.00030 | 0.002 | 0.08639 / 0.08635 |
| **40 dB** | **0.04064** | 0.04064 | **0.04190** | **0.970** | 0.00032 | 0.008 | 0.02854 / 0.02853 |

**Width scaling.** The empirical W₃₀/W₄₀ is 3.091, against the predicted √10 = 3.162 (−2.3%). For the stress level, the empirical W₂₀/W₃₀ is 3.70 against 3.16.

**Selection probabilities** (the prediction is version B, the exact deterministic gap):

| η | −0.2 | −0.15 | −0.1 | −0.075 | −0.05 | −0.025 | 0 | 0.025 | 0.05 | 0.075 | 0.1 | 0.15 | 0.2 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 20 dB P(A) emp | 0.790 | 0.758 | 0.670 | 0.600 | 0.584 | 0.500 | 0.464 | 0.394 | 0.348 | 0.299 | 0.254 | 0.182 | 0.090 |
| 20 dB P(A) pred | 0.896 | 0.828 | 0.736 | 0.682 | 0.624 | 0.563 | 0.500 | 0.437 | 0.377 | 0.319 | 0.265 | 0.174 | 0.107 |
| 30 dB P(A) emp | 1.000 | 1.000 | 0.976 | 0.943 | 0.830 | 0.687 | 0.510 | 0.299 | 0.157 | 0.068 | 0.032 | 0.003 | 0.000 |
| 30 dB P(A) pred | 1.000 | 0.999 | 0.977 | 0.933 | 0.841 | 0.691 | 0.500 | 0.309 | 0.160 | 0.068 | 0.024 | 0.002 | 0.000 |
| 40 dB P(A) emp | 1.000 | 1.000 | 1.000 | 1.000 | 0.999 | 0.939 | 0.508 | 0.066 | 0.001 | 0.000 | 0.000 | 0.000 | 0.000 |
| 40 dB P(A) pred | 1.000 | 1.000 | 1.000 | 1.000 | 0.999 | 0.943 | 0.500 | 0.058 | 0.001 | 0.000 | 0.000 | 0.000 | 0.000 |

**The local law performs as well as the exact gap** in this range: mean absolute errors of 0.0047 and 0.0016 in both versions. Δ₀(η) is close to linear for |η| ≤ 0.2, and ‖r_A − r_B‖ varies by only ±9%. The anticipated "local law plus finite-range refinement" split is therefore not needed within the tested offsets.

**Descriptive statistics** (not preregistered; no threshold applied):
- **Binomial z-scores of the observed proportions against the prediction.**
  - 30 dB: all |z| ≤ 1.74, Σz² = 9.5 over 13 points;
  - 40 dB: all |z| ≤ 1.14, Σz² = 1.9 over the 7 points with 0 < p < 1.

  At 30 and 40 dB the prediction is statistically indistinguishable from the Monte Carlo.
- **20 dB.** P(A) is over-predicted, by up to 0.106 at η = −0.2. There, 1–10% of records are OUT: displaced A-family fits beyond the one-unit radius, which the two-branch theory does not represent and the preregistered label counts as "not A". The B-selection probability P(B) = 1 − P(A) − P(OUT) is still predicted within 0.032 at every η; this comparison is post hoc. The predicted 20 dB width, 0.4064, is 0.97 of the post-hoc nearest-root width, 0.418 (manuscript S7).

## 3. Global MSE

Frequency-pair squared error in u units, relative to the noiseless winner: A for η < 0, B for η ≥ 0, following the frozen convention.

| η | −0.2 | −0.15 | −0.1 | −0.075 | −0.05 | −0.025 | 0 | 0.025 | 0.05 | 0.075 | 0.1 | 0.15 | 0.2 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 20 dB emp | 18.4 | 28.6 | 42.3 | 56.8 | 58.5 | 72.0 | 83.3 | 69.8 | 61.8 | 52.3 | 44.0 | 31.2 | 16.4 |
| 20 dB pred | 17.5 | 28.5 | 43.4 | 52.2 | 61.6 | 71.5 | 81.6 | 71.4 | 61.4 | 51.9 | 43.2 | 28.4 | 17.4 |
| 20 dB emp/pred | 1.05 | 1.00 | 0.97 | 1.09 | 0.95 | 1.01 | 1.02 | 0.98 | 1.01 | 1.01 | 1.02 | 1.10 | 0.94 |
| 30 dB emp | 0.0391 | 0.0396 | 3.93 | 9.32 | 27.6 | 50.8 | 83.0 | 48.6 | 25.6 | 11.0 | 5.18 | 0.487 | 0.00655 |
| 30 dB pred | 0.0456 | 0.266 | 3.78 | 11.0 | 26.0 | 50.4 | 81.4 | 50.3 | 26.0 | 11.1 | 3.85 | 0.253 | 0.0134 |
| 30 dB emp/pred | 0.86 | 0.15 | 1.04 | 0.85 | 1.06 | 1.01 | 1.02 | 0.97 | 0.98 | 0.99 | 1.35 | 1.92 | 0.49 |
| 40 dB emp | 0.00395 | 0.00381 | 0.00397 | 0.00374 | 0.166 | 9.93 | 82.7 | 10.7 | 0.163 | 0.000691 | 0.000781 | 0.000628 | 0.000672 |
| 40 dB pred | 0.00399 | 0.00389 | 0.00379 | 0.00391 | 0.134 | 9.34 | 81.4 | 9.38 | 0.135 | 0.000936 | 0.000722 | 0.000693 | 0.000665 |
| 40 dB emp/pred | 0.99 | 0.98 | 1.05 | 0.96 | 1.24 | 1.06 | 1.02 | 1.15 | 1.21 | 0.74 | 1.08 | 0.91 | 1.01 |

**Decomposition.**
- **η = 0.** The prediction is 81.6, 81.4 and 81.4 at 20, 30 and 40 dB. The between-mode term is ½‖u_A − u_B‖² = 81.42; the within-mode term is 0.219, 0.0219 and 0.00219.
- **A-referenced values at η = 0** (the table caption's quantity): empirical/predicted is 0.970, 0.978 and 0.983.
- **Away from the exchange,** the prediction reduces to the local trace, and the observed error follows it. For example, 40 dB at η ≤ −0.075 gives ratios 0.96–1.05.

**Tails at 30 dB.** Here the MSE is set by rare mis-selections, each costing about 163 u²:
- at η = −0.15, 1.4 B selections were expected and 0 observed (Poisson probability 0.25), giving the ratio 0.15;
- at η = +0.15, 1.5 A selections were expected and 3 observed, giving 1.92;
- at η = +0.2, 0.04 were expected and 0 observed, giving 0.49.

These are sampling limits of a 1000-record Monte Carlo, not model misfits. All three points are inside the pre-declared order-of-magnitude band.

**Consistency check.** The recomputed tr C_A and tr C_B equal the frozen `local_trace_A` and `local_trace_B` to 8 × 10⁻¹⁴ relative.

## 4. Pre-declared criteria (00_THEORY_PLAN.md §8)

| Criterion | 30 dB | 40 dB | 20 dB (stress, not used) |
|---|---|---|---|
| **C1** orientation, probit slope < 0, \|Δη₅₀\| ≤ W/2 | pass (slope −19.8; 0.002 W) | pass (slope −61.2; 0.008 W) | pass (0.06 W) |
| **C2** \|W_pred/W_emp − 1\| ≤ 0.25 | pass (−0.8%) | pass (−3.0%) | pass (−15.2%) |
| **C3** mean \|ΔP(A)\| ≤ 0.08 | pass (0.0047) | pass (0.0016) | pass (0.045) |
| **C4** 0.1 ≤ emp/pred ≤ 10 at all 13 η | pass (0.149–1.925) | pass (0.739–1.244) | pass (0.944–1.100) |
| **C5** at least 5 of 9 central points in [0.5, 2] | pass (9/9) | pass (9/9) | pass (9/9) |
| **F2** orientation or crossing wrong | no | no | no |
| **F3** at least 5 of 9 central points outside [0.1, 10] | no (0) | no (0) | no (0) |

**F1**, the width-scaling test (not blind), is not triggered: (W₃₀/W₄₀)_emp/(W₃₀/W₄₀)_pred − 1 = −2.3%.

C1–C5 all pass at both 30 and 40 dB, and no F condition holds. **Outcome: STRONG_PREDICTION.**

## 5. Interpretation and limits

**Result.** A first-order law built only from noiseless quantities predicts, at 30 and 40 dB and within binomial sampling error:
- the whole selection curve, without any fitted parameter;
- its 10–90% width (to 1% and 3%);
- the global error near the exchange (factor 0.74–1.35 for |η| ≤ 0.1).

The noiseless quantities are the finite-spacing cost-gap slope g, the residual difference ‖r_A − r_B‖, the branch minimizers and the local covariances. The law implies W ∝ σ: each 10 dB narrows the transition by √10.

**Why it works.** The neglected re-optimization bias is tiny: tr(H⁻¹G) is 5.60 for A and 5.75 for B, so the crossing shift is at most 0.3% of W (`second_order_diagnostic.json`, sealed). The selection–error correlation, which is not modelled, is not visible at this precision.

**Limits.**
- It is a first-order approximation for one specified configuration: N = 21, z = 2, b = 10, φ = π. It is not a theorem, not a bound, and not a claim about other configurations.
- The two-branch model has no term for displaced A-family fits. At 20 dB, 1–10% of records are such fits, and P(A) is over-predicted there (stress test).
- Tail MSE values are rare-event quantities.

## 6. Relation to the manuscript

**Current manuscript (Section VIII, "Noisy estimation: local versus global behaviour").** It reports:
- the preregistered Monte Carlo;
- the conditional accuracy of the sandwich covariance;
- the dominance of mode selection;
- the empirical widths, compared only with the archived experiment.

It gives no predictive law for P(A), the width or the global MSE.

**What this analysis adds:**
1. a closed-form, parameter-free selection law, with a width law W = 2Φ⁻¹(0.9)√2σ‖r_A − r_B‖/(g ε_c) that scales as σ;
2. a two-branch prediction of the winner-referenced global MSE, showing that the between-mode term supplies at least 99% of it near the exchange;
3. a sealed, blind validation against the preregistered 1000-record run.

An earlier internal draft had a frozen-projector version of (1) with the same predicted widths, compared non-blind with the archived data. That draft is not part of the current manuscript.

## 7. Recommendation

By the pre-declared rule, STRONG_PREDICTION permits a concise derivation in the main paper, together with predicted-vs-empirical widths, the two-branch MSE prediction, and a reframing of Section VIII around prediction.

I recommend doing this as a focused rewrite of Section VIII:
- two displayed equations;
- one extended table;
- theory curves added to the existing Fig. 7;
- matching one-line edits to the evidence table, the contribution list and, optionally, the abstract.

It should not become a new theorem-level claim. The evidence level is "first-order approximation validated numerically, prediction frozen before comparison", and 20 dB stays a qualitative stress case.

## 8. Suggested manuscript changes (exact text; NOT applied)

All edits are to non-frozen copies. The frozen `body_reviewfriendly.tex`, `noise_table.tex`, `evidence_table.tex`, `S7_noise_methods.tex` and figure scripts would be revised only after approval. The frozen PDFs would then be regenerated (see §9).

### 8.1 Section title (`body_reviewfriendly.tex`, Section VIII)

```latex
\section{Noisy estimation: predicting mode selection and global error}\label{sec:noise}
```

### 8.2 New text after the local-approximation paragraph (after `\eqref{eq:noise-local}`)

```latex
\emph{Predicting mode selection.} Let $r_j$ be the residual of the noiseless branch-$j$ fit and
$\Delta J_0(\epsilon)=J_{\mathcal A}-J_{\mathcal B}$ the noiseless cost gap. Each branch cost is stationary in its own
parameters, so noise changes it to first order by $2\Re\langle r_j,w\rangle$, and the noisy gap is approximately Gaussian
with mean $\Delta J_0$ and standard deviation $\sqrt2\,\sigma\|r_{\mathcal A}-r_{\mathcal B}\|_2$. Hence
\begin{equation}\label{eq:noise-select}
P(\mathcal A)\approx\Phi\Big(-\frac{\Delta J_0(\epsilon)}{\sqrt2\,\sigma\|r_{\mathcal A}-r_{\mathcal B}\|_2}\Big)
\approx\Phi(-K\eta),\qquad K=\frac{g\,\epsilon_c}{\sqrt2\,\sigma\|r_{\mathcal A}-r_{\mathcal B}\|_2},
\end{equation}
where, by the envelope theorem, $g=d\Delta J_0/d\epsilon=2\Re\langle r_{\mathcal A}-r_{\mathcal B},\partial x/\partial\epsilon\rangle$
at $\epsilon_c$. The 10--90\% selection width is $W=2\Phi^{-1}(0.9)/K\propto\sigma$, so each 10~dB narrows the
transition by $\sqrt{10}$. At $N=21$, $z=2$, $g=8.259$ and $\|r_{\mathcal A}-r_{\mathcal B}\|_2=1.307$, giving
$W=4.06\times10^{-\mathrm{SNR}/20}$.

\emph{Predicting the global error.} If the fitted pair of the selected branch is distributed approximately as its local
approximation~\eqref{eq:noise-local}, the error relative to the noiseless winner $u_W$ is
\begin{equation}\label{eq:noise-mix}
E\|\widehat u-u_W\|^2\approx\sum_{j\in\{\mathcal A,\mathcal B\}}P(j)\,\big[\operatorname{tr}C_j+\|u_j-u_W\|^2\big].
\end{equation}
Its between-mode part is fixed by the noiseless branch geometry ($\|u_{\mathcal A}-u_{\mathcal B}\|^2=162.8$ at
$\epsilon_c$). Both predictions use only noiseless quantities; they were computed and frozen before the Monte Carlo
outcomes were examined (Supplementary Sec.~S7).
```

### 8.3 Replace the width sentence in the "near the exchange" paragraph

Current text: "The affected range of $\eta$ narrows as the noise decreases: the 10--90\% selection widths are about $0.13$ at 30~dB and $0.04$ at 40~dB, consistent with the archived experiment ($0.129$ and $0.039$)."

Proposed:

```latex
The affected range of $\eta$ narrows as predicted: the observed 10--90\% widths are $0.1295$ at 30~dB and $0.0419$ at
40~dB, against $0.1285$ and $0.0406$ from \eqref{eq:noise-select}, and their ratio $3.09$ is close to $\sqrt{10}$.
Across the 13 offsets the predicted $P(\mathcal A)$ differs from the observed proportions by $0.005$ (30~dB) and
$0.002$ (40~dB) on average, within binomial sampling error, and \eqref{eq:noise-mix} reproduces the global mean squared
error within a factor $0.74$--$1.35$ for $|\eta|\le0.1$, where the between-mode term supplies more than $99\%$ of it.
```

Keep the existing 20 dB sentence and add:

```latex
The two-branch model does not represent these displaced fits and overpredicts $P(\mathcal A)$ by up to $0.11$ at 20~dB,
although its global-error prediction remains within about $10\%$ there.
```

### 8.4 Closing paragraph of the section

Add after "...controls the error.":

```latex
The same noiseless quantities that locate the exchange, the cost-gap slope and the residual difference of the two
branches, predict the selection probabilities, the width of the affected region, and the resulting global error.
```

### 8.5 `noise_table.tex`

Add two columns: "Width (obs./pred.)" and "Global/mixture, |η| ≤ 0.1".

| SNR | P(A), η = 0 | Cond./local | Global/local | Width, obs./pred. | Global/mixture, \|η\| ≤ 0.1 |
|---|---|---|---|---|---|
| 20 dB | 0.46 | 0.63–1.07ᵃ | 221 | 0.42–0.48ᵇ / 0.41 | 0.95–1.09 |
| 30 dB | 0.51 | 0.98–1.09 | 2217 | 0.1295 / 0.1285 | 0.85–1.35 |
| 40 dB | 0.51 | 0.91–1.08 | 22287 | 0.0419 / 0.0406 | 0.74–1.24 |

Caption addition: "Pred.: \eqref{eq:noise-select}; Global/mixture: Monte Carlo global error divided by \eqref{eq:noise-mix}."

### 8.6 Fig. 7

- Panel (a): add the predicted curves Φ(−ΔJ₀/(√2σ‖r_A − r_B‖)) for each SNR.
- Panels (b) and (c): add the prediction of \eqref{eq:noise-mix}.
- Caption addition: "Lines: first-order predictions \eqref{eq:noise-select} and \eqref{eq:noise-mix}, frozen before comparison."

The data for these curves are `results/transition_theory/deterministic/predictions.csv`.

### 8.7 `evidence_table.tex`, last row

```latex
Noisy mode selection and global error & $N=21$, $z=2$; 20--40 dB & Numerical Monte Carlo; first-order prediction
frozen before comparison & First-order approximation, not a bound; 20 dB displaced $\mathcal A$-family fits not modelled.\\
```

### 8.8 Introduction

**Contribution item 5:**

```latex
\item an estimator-level consequence: each branch's local covariance stays accurate conditionally, while near the
exchange a first-order cost-gap law predicts the mode-selection probabilities, their transition width, and the resulting
global error.
```

**Roadmap sentence** (end of the paragraph beginning "Because the certified instances…"): append

```latex
; a first-order law computed from the noiseless branches predicts both the mixture weights and the resulting error
```

before "(Section~\ref{sec:noise})".

### 8.9 Abstract (optional)

The abstract is 200 whitespace tokens; the replacement sentence keeps it at 199. Replace "In noisy Monte Carlo experiments, a branch-conditioned local covariance approximation predicts the error well away from the exchange, whereas near equal global competition the overall error is dominated by mode selection and is not summarized by any single local pseudo-true mode." with:

"In noisy Monte Carlo experiments, a branch-conditioned local covariance approximation predicts the error away from the exchange; near it, a first-order law computed from the noiseless branches predicts the mode-selection probabilities, their noise-proportional transition width, and the resulting global error."

### 8.10 Supplement S7

Add a paragraph "Frozen first-order prediction" after the transition-width paragraph. It should state:
- the derivation (DERIVATION.md §1);
- the pre-declared criteria;
- the seal digests;
- the outcome: STRONG_PREDICTION at 30 and 40 dB, with 20 dB a stress test;
- that the 20 dB overprediction of P(A) comes from displaced A-family fits.

It should point to `paper/transition_theory/` and `results/transition_theory/` in the reproducibility archive.

## 9. Consequences for the submission package, if the changes are approved

Everything frozen on 2026-09-28 would need updating:
- the named main PDF: a new page count, currently 11 of 13, likely still ≤ 12;
- the double-column supplement;
- the SHA-256 lists and the package record;
- possibly the cover letter's contribution 3 and the portal abstract, if §8.9 is taken;
- the reviewer archive, which must be rebuilt: it has to include the new analysis, but its selection rules currently exclude `paper/`;
- the public GitHub repository, which is committed locally but not yet pushed. Pushing it now and adding the analysis in a second commit, or rebuilding it before the first push, are both possible.

## 10. Confirmation: no new Monte Carlo and no tuning

- **No new noisy records.** No SNR, η, ε, N, z, signal family, model order, label or labelling radius was added or changed, and no second predictive model was built.
- **No fitted parameters.** No quantity was fitted to Monte Carlo outcomes; the only fits are the frozen preregistered probits, which were read, not re-run.
- **The only random computation** was the planned synthetic check of Var[2 Re⟨d, w⟩] on a fixed vector: 2 × 10⁶ draws, seed `SeedSequence([20260928, 99, 1])`, with no optimizer and no outcome.
- **The sealed files are unchanged.** SEAL_2 still verifies after the validation.
- **No other project file changed.** Only `paper/transition_theory/` and `results/transition_theory/` were written. A read-only re-hash confirms the frozen state: 766 frozen files (732 unchanged, 34 changed, the same set as the final pre-submission record, 0 missing), all 51 Stage 24 deliverables, and the submission hashes.
- **No manuscript file was modified.**
