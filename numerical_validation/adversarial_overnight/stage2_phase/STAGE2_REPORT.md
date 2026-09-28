# Stage 2 — record-origin / strong-pair phase invariance of the z² law

Pre-registered classification (z² law only): **INCONCLUSIVE** — reason below. Existence of the exchange in general is assessed separately.

## 2A. Spacing-scaled offset δ = κz (weak tone fixed at φ = π)

Tangent level: independent realified least-squares implementation (not the frozen closed forms), N = 21. Finite level: full-torus crossing analysis at z ∈ {0.25, 0.5, 1, 2}, ε = λz², λ ∈ [0, 1] (+ one boundary extension when no qualifying crossing).

| κ | tangent λ₂₁(κ) | v_A | v_B | coalescent margin | third margin | finite ε_c/z² at z = 0.25 / 0.5 / 1 / 2 | rel. dev. from λ₂₁(κ) at z = 0.25 / 0.5 / 1 | qualifying at z = 0.25 / 0.5 / 1 / 2 |
|---|---|---|---|---|---|---|---|---|
| +0 | 0.06598041 | -4.1534 | 12.4681 | 1.805e-02 | 0.0644 | 0.065917 / 0.065727 / 0.064971 / 0.062048 | -0.096% / -0.385% / -1.529% | yes / yes / yes / yes |
| +0.05 | 0.07091729 | -3.4910 | 12.3484 | 1.532e-02 | 0.0747 | 0.070835 / 0.070591 / 0.069621 / 0.065902 | -0.115% / -0.461% / -1.827% | yes / yes / yes / yes |
| +0.1 | 0.08551103 | -2.3937 | 12.0707 | 1.114e-02 | 0.1096 | 0.085386 / 0.085010 / 0.083527 / 0.077863 | -0.147% / -0.586% / -2.321% | yes / yes / yes / yes |
| +0.2 | 0.13222566 | -1.0436 | 11.5605 | 5.063e-03 | 0.2720 | 0.131969 / 0.131201 / 0.128171 / 0.116672 | -0.194% / -0.775% / -3.066% | yes / yes / yes / yes |
| +0.3 | 0.18822386 | -0.3655 | 11.2497 | 1.202e-03 | 0.5710 | 0.187760 / 0.186375 / 0.180925 / 0.160373 | -0.246% / -0.982% / -3.878% | no / no / yes / yes |
| +0.3882 | no regular A/B root (A at the central chart) | — | — | — | — | 0.239646 / 0.237445 / 0.228774 / 0.195821 | — / — / — | no / no / no / no |
| +0.4 | 0.24744681 | 0.0389 | 11.0571 | 2.267e-05 | 1.0134 | 0.246650 / 0.244305 / 0.235078 / 0.200027 | -0.322% / -1.270% / -4.998% | no / no / no / no |
| +0.5 | 0.30806419 | 0.3076 | 10.9282 | 2.130e-03 | 1.5995 | 0.306818 / 0.302833 / 0.286956 / 0.212345 | -0.404% / -1.698% / -6.852% | no / no / no / no |

Negative κ gives identical values (exact conjugation/time-reversal symmetry): tangent relative differences ≤ 2e-16, finite-record relative differences ≤ 2.6e-11 (M2A.3 pass).
Quadratic shift (M2A.4): λ₂₁(0.05) − λ₂₁(0) = 4.936882e-03 vs frozen prediction 1.92838866512018·κ² = 4.820972e-03 (rel. error 2.35%; 1.26% at κ = 0.1; -16.4% at κ = 0.2 where O(κ⁴) matters).

**2A (pre-registered M2A.1–M2A.4): pass for every |κ| ≤ 0.20.** A qualifying separated global exchange exists at every tested z, and ε_c/z² converges to the independently computed tangent λ₂₁(κ) (deviation 0.10–0.19% at z = 0.25, shrinking ∝ z²). For |κ| = 0.3 the small-z exchanges (z = 0.25, 0.5) fail the pre-registered regularity rule (A's within-fit separation 0.45–0.63 < 1.0 normalized units: A is near-coalescent) although they remain global switches; for |κ| ≥ 0.388 the A participant sits exactly on the confluent line (h = 0, unbounded two-tone amplitudes), matching the frozen analytical collision κ ≈ 0.3882 at N = 21.

## 2B. Fixed nonzero offset δ (weak tone at φ = π)

| δ/π | tangent ε/z | z = 0.05 | z = 0.1 | z = 0.2 | z = 0.4 | z = 0.8 (outside fit window) |
|---|---|---|---|---|---|---|
| 0.05 | 0.09815 | ε/z = 0.09752; A sep 1.35 (regular) | ε/z = 0.09701; A sep 1.09 (regular) | ε/z = 0.09633; A sep 0.58 (coalescent-class) | ε/z = 0.09645; A sep 0.00 (coalescent-class) | ε/z = 0.10220; A sep 1.60 (regular) |
| 0.1 | 0.19389 | ε/z = 0.19325; A sep 1.44 (regular) | ε/z = 0.19265; A sep 1.25 (regular) | ε/z = 0.19155; A sep 0.85 (coalescent-class) | ε/z = 0.18977; A sep 0.00 (coalescent-class) | ε/z = 0.18833; A sep 0.00 (coalescent-class) |
| 0.2 | 0.36880 | ε/z = 0.36826; A sep 1.41 (regular) | ε/z = 0.36770; A sep 1.18 (regular) | ε/z = 0.36650; A sep 0.54 (coalescent-class) | ε/z = 0.36309; A sep 0.00 (coalescent-class) | ε/z = 0.35275; A sep 0.00 (coalescent-class) |
| 0.3 | 0.50762 | ε/z = 0.50728; A sep 1.29 (regular) | ε/z = 0.50688; A sep 0.84 (coalescent-class) | ε/z = 0.50554; A sep 0.00 (coalescent-class) | ε/z = 0.49903; A sep 0.00 (coalescent-class) | ε/z = 0.48047; A sep 0.00 (coalescent-class) |
| 0.4 | 0.59674 | ε/z = 0.59674; A sep 0.78 (coalescent-class) | ε/z = 0.59609; A sep 0.00 (coalescent-class) | ε/z = 0.59100; A sep 0.00 (coalescent-class) | ε/z = 0.57956; A sep 0.00 (coalescent-class) | ε/z = 0.70223; A sep 3.74 (regular) |

Pre-registered 2B metric: exponent fit on the four in-window z using the *qualifying* crossing nearest the tangent pair. Qualifying crossings exist only at z ≤ 0.1 for δ/π ≤ 0.2, at z = 0.05 for δ/π = 0.3, and at no in-window z for δ/π = 0.4, because the lower-ε participant (A) has within-fit separation < 1.0 or lies on the confluent line. **The pre-registered exponent fit is therefore not evaluable → Stage-2 classification INCONCLUSIVE.**

### POST_HOC_SENSITIVITY (does not replace the pre-registered result)

Using every global switch into the B-like fit (satellite near 10.3), regardless of A's regularity, over z ∈ {0.05, 0.1, 0.2, 0.4}:

| δ/π | α (log–log) | rel. RMS residual ε_c = C₁z | rel. RMS residual ε_c = C₂z² | α without smallest z | α without largest z |
|---|---|---|---|---|---|
| 0.05 | 0.9942 | 0.0049 | 0.5827 | 0.9959 | 0.9912 |
| 0.1 | 0.9913 | 0.0069 | 0.5850 | 0.9891 | 0.9936 |
| 0.2 | 0.9934 | 0.0055 | 0.5845 | 0.9909 | 0.9965 |
| 0.3 | 0.9925 | 0.0066 | 0.5853 | 0.9887 | 0.9975 |
| 0.4 | 0.9861 | 0.0118 | 0.5879 | 0.9797 | 0.9930 |

The global switch location scales as ε_c ∝ z (α = 0.986–0.994) and ε_c/z matches the tangent prediction 0.6274474|sin δ| to better than 1% at z = 0.05. The z² law does **not** hold at fixed δ; the fixed-δ exchange is a first-order (ε ∝ z) phenomenon whose lower participant becomes near-coalescent already at z ≈ 0.1–0.2.

## 2C. Reinterpretation as window placement

A shift of the record centre by τ samples gives κ = τ/N (fraction of record = κ). A *true* shift also rotates the weak tone by bκ rad.

| κ | fraction of record | τ (N=21) | τ (N=31) | τ (N=101) | τ (N=1001) | weak-tone rotation bκ (rad) |
|---|---|---|---|---|---|---|
| 0.00 | 0.00 | 0.00 | 0.00 | 0.0 | 0 | 0.00 |
| 0.05 | 0.05 | 1.05 | 1.55 | 5.1 | 50 | 0.50 |
| 0.10 | 0.10 | 2.10 | 3.10 | 10.1 | 100 | 1.00 |
| 0.20 | 0.20 | 4.20 | 6.20 | 20.2 | 200 | 2.00 |
| 0.30 | 0.30 | 6.30 | 9.30 | 30.3 | 300 | 3.00 |
| 0.40 | 0.40 | 8.40 | 12.40 | 40.4 | 400 | 4.00 |
| 0.50 | 0.50 | 10.50 | 15.50 | 50.5 | 500 | 5.00 |

Secondary (pre-registered) tangent check of a *true* record shift (δ = κz and φ = π + bκ): a valid global A/B tangent crossing exists for κ = 0.05 (λ₂₁ = 0.07294), 0.10 (λ₂₁ = 0.09087), 0.20 (λ₂₁ = 0.13895), 0.30 (λ₂₁ = 0.20804); for |κ| ≥ 0.4 the seeded solver returned no valid global A/B pair (degenerate or non-global roots) — recorded as FAILED for this secondary diagnostic.

## Interpretation for the z² claim

* The quadratic law is **locally stable in the manuscript's κ-scaled sense**: it survives record-centre displacements up to 0.2N samples (4.2 samples at N = 21) with λ(κ) = λ(0) + 1.93κ² + O(κ⁴), numerically at finite z to 0.1% agreement.
* It is **alignment-specific**: it needs |δ| ≲ 0.2–0.3·z. For a fixed relative phase — the generic case when the beat period far exceeds the record (small z) — the crossing scales as ε ∝ z. The fraction of relative phases compatible with the z² regime shrinks ∝ z as z → 0.
* The pre-registered label stays **INCONCLUSIVE** only because the fixed-δ fit required *regular* (separation ≥ 1.0) participants, which fixed-δ exchanges do not provide beyond z ≈ 0.1. The post-hoc evidence points to `QUADRATIC_LAW_LOCALLY_STABLE_BUT_ALIGNMENT_SPECIFIC`; it does not override the pre-registered label.

Evidence levels: tangent rows GLOBAL_NUMERICAL (1-D full-period scan); finite rows GLOBAL_NUMERICAL; post-hoc fits POST_HOC_SENSITIVITY. Figures: `stage2_phase/phase_scaling.png|pdf` (pre-registered; its fixed-δ legend shows α = n/a because the pre-registered fit was not evaluable) and `stage2_phase/POSTHOC_fixed_delta_scaling.png` (POST_HOC_SENSITIVITY; filled markers = regular A participant, open = near-coalescent).