# Evidence levels

This table states, for each result, how it is supported and how far it has been replayed. It follows Table I of the manuscript and Supplementary Sections S5 and S8.

The replays run the chains R1–R7 and T1 of the replay tool (`replay_tools/replay_and_compare.py` in the reviewer archive, `replay/replay_and_compare.py` in the repository); chain definitions are in the README.

## Level definitions

- **Certified.** An analytical proof whose computational hypotheses are verified by complete directed-rounding interval computations (mpmath.iv, the *primary* implementation). The computation covers the whole stated domain; nothing is inferred from a sampled grid.
- **Independent replay.** A second implementation re-evaluated every proof obligation of the archived certificate:
  - it uses Arb ball arithmetic through python-flint 0.9.0 at 320-bit precision;
  - a static audit shows that it imports none of the primary evaluators;
  - it reuses the archived candidate partitions (cell coordinates and parameter brackets) and re-checks their coverage exactly.

  It is not a proof-assistant formalization and not third-party review.
- **Original-code replay.** The frozen primary scripts were rerun unmodified in a clean copy of the project, with every expected output deleted first. The rerun outputs were compared with the archived records: on Windows 11 (Stage 0 of the adversarial validation, 2026-09-28) and on Linux/WSL2 (2026-09-28).
- **Global numerical robustness.** Full-torus numerical global searches run under the preregistered adversarial protocol (`validation/adversarial_overnight/PREREGISTRATION.md`). They are not interval certificates.
- **Validated first-order approximation.** An analytical approximation computed only from deterministic (noiseless) quantities and the prescribed noise level.
  - Its plan, pass/fail criteria, code and predictions were sealed with SHA-256 records before the preregistered Monte Carlo outcomes it predicts were read.
  - It was then compared with them once.
  - It is not a theorem, certificate or bound.
- **Numerical, outside the protocol.** Earlier or supplementary numerical computations: not preregistered and not certified.
- **Exploratory.** Results the manuscript itself labels exploratory; context only.
- **Post-hoc sensitivity.** Analyses chosen after seeing preregistered results. They are reported, but they do not change any preregistered classification.

## Results

| Result | Scope | Evidence level | Independent (Arb) replay | Original-code replay | Main records |
|---|---|---|---|---|---|
| Global exchange of two separated, regular two-tone fits at small spacing (checkable criterion C1–C8) | Centred phase-aligned strong pair (δ = 0), weak tone at b = 10 with phase φ = π; N = 11, 15, 21, 31, 41 | **Certified.** C1–C7 are interval-certified; C8 is analytical. The small-spacing radius is existential. | No | Yes (R2) | `stage12_generalN/` |
| Full-line continuum tangent crossing | Continuum limit at fixed normalized geometry (z = NΔ and b = Nω₃ fixed); b = 10, φ = π | **Certified** (27,708-cell directed partition of [−100, 100] plus analytic tail and coalescent gap) | Yes: the partition is re-checked uniformly in ρ = N⁻² ∈ [0, 35,377⁻²], which includes the continuum ρ = 0 (R4) | Yes (R3) | `continuum_global_certificate.json`, `stage14_continuum_*`, `independent_N0_replay.*` |
| Large-N transfer: every odd N ≥ 10,001 | **Fixed normalized geometry only** (z = NΔ, b = Nω₃ fixed); not fixed physical frequencies | **Certified** (analytic transfer plus directed parameter-interval certificate over 0 ≤ N⁻² ≤ 10,001⁻²) | **No** for the 10,001 floor. The Arb replay covers only the earlier, more conservative floor N ≥ 35,377. | Yes (R5) | `tightened_N0_certificate.json`, `largeN_condition_thresholds.csv` |
| Finite-spacing global exchange at z = 2 | N = 21 and N = 31 crossing brackets | **Certified** (full-torus directed-interval certificates) | Yes: all 65,337 outer cells, 2,801 inner cells and 86 root/crossing obligations (68,224 checks), on the archived partitions (R1) | The Arb replay itself was rerun (R1). The original certificates are frozen records. | `validated_globality_certificate.json`, `N31_global_certificate.json`, `stage12_replay/` |
| N = 21 finite-width amplitude interval | ε ∈ [0.213190630172278, 0.283190630172277] at z = 2 | **Certified** (sufficient interval, not proved maximal) | **No** | Yes (R6, starting with `stage17_epsilon_refine.py 0.035 14287`) | `finite_width_epsilon_certificate.json`, `epsilon_interval_obligations.csv` |
| N = 21 coalescent endpoint never globally optimal for λ > 0 | N = 21, b = 10, φ = π | **Certified** constants (directed interval) with an analytic argument | **No** | Yes (R7) | `lower_transition_certificate.json` |
| Critical law ε_c = λ_N z² + c_N z⁴ + O(z⁶) | Only the centred, phase-aligned family (δ = 0, or δ = χz with small \|χ\|). **Not universal.** | Theorem. λ_N and c_N are interval-enclosed in the certificates of the five finite lengths (`local_crossing.lambda_quartic`). The large-N expansion coefficients a₂ and a₄ are high-precision evaluations of analytic formulas, not interval enclosures. | No | via R2 | `stage12_generalN/`, `stage14_asymptotic_coefficients.json` |
| Strong-pair phase-offset stability, δ = χz | Each certified finite N; continuum | Theorem with an **existential** radius χ₀(N); no explicit or N-uniform radius | — | — | Supplement S6 |
| Fixed nonzero strong-pair phase | Fixed δ ≠ 0, z → 0 | Analytical O(z) scale; the tangent coefficient is numerical; the preregistered Stage-2 classification is **inconclusive** | — | — | `validation/adversarial_overnight/stage2_phase/` |
| Weighting windows: Hann, Hamming, DPSS, Hann² | N = 21, z = 2 | **Global numerical robustness.** The crossing amplitude depends on the window. | — | — | `validation/.../stage1_window/` |
| Weak-tone location | b ∈ [9.5, 10.5] (finite N and continuum; preregistered); integer 6 ≤ b ≤ 13 (continuum; earlier numerical map) | **Global numerical robustness only. [9.5, 10.5] is not a certified b interval.** The certified statement is an existential neighbourhood of b = 10. An explicit b-interval certificate was attempted: every local obligation passed, but the attempt stopped at an efficiency limit, with no counterexample. | — | — | `validation/.../stage3_bwidth/`, `lambda_inf_vs_b.csv` |
| Strong-amplitude imbalance | 1 ∓ ι, \|ι\| ≤ 0.2 | **Global numerical robustness**; an existential stability corollary only | — | — | `validation/.../stage4_imbalance/` |
| Phase / record-origin offsets | δ = χz with \|χ\| ≤ 0.2; δ ≤ 0.6 rad at z = 2 | **Global numerical robustness**, plus the existential theorem above | — | — | `validation/.../stage2_phase/` |
| Matching conventions across N | Fixed normalized versus fixed physical frequencies | Global numerical (preregistered) | — | — | `validation/.../stage5_matching/` |
| Noisy estimation: local covariance versus global error | N = 21, z = 2; 20, 30, 40 dB; 39,000 records | Numerical Monte Carlo. A local approximation, not a bound; the 20 dB width is labelling-sensitive. | — | — | `validation/.../stage6_noise/` |
| Noisy branch selection and global error near the exchange: first-order law P(A) ≈ Φ(−ΔJ₀/(√2σ‖r_A − r_B‖₂)) ≈ Φ(−Kη), and two-branch MSE | N = 21, z = 2; validated at 30 and 40 dB; 20 dB is a stress case | **Validated first-order approximation**, outcome STRONG_PREDICTION under the pre-declared criteria. At 30 and 40 dB:<br>• predicted widths 0.1285 and 0.0406 against observed 0.1295 and 0.0419;<br>• mean \|ΔP(A)\| 0.005 and 0.002;<br>• global MSE within a factor 0.74–1.35 for \|η\| ≤ 0.1.<br>Sealed before the per-setting outcomes were read, but the widths were known beforehand, so the width comparison is not blind. At 20 dB, displaced A-family fits (1–10%) fall outside the approximation. | — | The deterministic predictions are recomputed by T1. The one-shot validation is not rerun. | `paper/transition_theory/`, `results/transition_theory/` |
| Solver, precision and tolerance invariance | All preregistered stages | Numerical (preregistered) | — | — | `validation/.../stage7_numinv/` |
| Model-order (BIC) check | 1,080 records, 20–40 dB | Numerical, outside the protocol; order three selected in every record | — | — | `order_selection_sanity.*` |
| Real-valued sinusoids | Two carriers | **Exploratory** | — | — | `validation/.../stage8_real/` |
| Archived spacing sweep (2 ≤ z ≤ 2π) and weak-tone phase map | — | Numerical, outside the protocol (archived maps; a few representative roots are interval-checked) | — | — | `phase_crossing_map.csv`, `resolution_sweep.csv` |
| Numerical z-continuation from small z to z = 2 | 44 points | Numerical, outside the protocol; illustration only; **no certified connected branch** | — | — | `p_bridge_numerical.csv`, `z_continuation_vs_asymptotics.csv` |
| Fixed-δ finite-record exponent fits (Fig. S3); 20 dB nearest-root labelling | — | **Post-hoc sensitivity** | — | — | `validation/.../stage2_phase/POSTHOC_*`, `stage6_noise/` |

## Not claimed

- A certified b interval.
- A certified connected branch from small spacing to z = 2.
- Behaviour at fixed physical frequencies as N grows.
- A universal z² law.
- A formal CRB or MCRB failure.
- A theorem or bound for noisy branch selection: the selection law and the two-branch error are a validated first-order approximation for one configuration.
- A failure of model-order selection.
- Independent Arb replay of every certified result: the N ≥ 10,001 floor, the N = 21 finite-width amplitude interval, the five tangent instances and the lower-endpoint constants have original-code replay only.
