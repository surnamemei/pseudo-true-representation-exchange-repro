# Sealed first-order transition theory (Section VIII of the manuscript)

This folder holds the one approved post-validation analysis. It is a first-order law that predicts, from noiseless branch quantities and the prescribed noise level alone:
- the probability of selecting branch A near the certified N = 21, z = 2 exchange,
  P(A) ≈ Φ(−ΔJ₀/(√2σ‖r_A − r_B‖₂)) ≈ Φ(−Kη);
- its 10–90% transition width, W = 2Φ⁻¹(0.9)/K, which scales as σ;
- the global mean squared error, through a two-branch mixture.

**Evidence level: validated first-order approximation.** It is not a theorem, certificate or bound.

## Provenance

The plan, pass/fail criteria, code and every predicted value were sealed with SHA-256 records before the per-setting outcomes of the preregistered Stage-6 Monte Carlo (`numerical_validation/adversarial_overnight/stage6_noise/`) were read. The comparison was then run once.

- **Not blind.** The empirical widths were already known, so the width comparison is not a blind test.
- **No new data.** No new Monte Carlo record was generated. The two-million-draw Gaussian check verified the variance formula of the implementation only.
- **Outcome under the pre-declared criteria:** STRONG_PREDICTION at 30 and 40 dB.
  - Predicted widths are 0.1285 and 0.0406, against observed 0.1295 and 0.0419.
  - The mean |ΔP(A)| is 0.005 and 0.002.
  - The global MSE is within a factor 0.74–1.35 for |η| ≤ 0.1.
- **20 dB** was declared a stress case in advance. There, displaced A-family fits fall outside the two-branch approximation.

## Contents

| Path | What it is |
|---|---|
| `00_THEORY_PLAN.md`, `transition_spec.json` | The plan and machine-readable spec, sealed first (`SEAL_1.*`) |
| `DERIVATION.md` | Derivation, constants and implementation checks, sealed with the code and predictions (`SEAL_2.*`) |
| `VALIDATION_REPORT.md`, `RESULT_RECORD.*` | The one-shot comparison and the hash record of its outputs |
| `code/transition_theory.py`, `code/second_order_diagnostic.py` | Deterministic predictions and checks |
| `code/validate.py` | The one-shot comparison with the Monte Carlo. It refuses to run a second time. |
| `results/deterministic/` | Sealed predictions (`predictions.csv`), constants (`deterministic.json`) and checks |
| `results/validation/` | Per-setting comparison (`comparison.csv`), criteria evaluation (`validation.json`) and figure |

The seal files record project-relative paths, for example `paper/transition_theory/...` and `results/transition_theory/...`. `replay/workspace_map.csv` maps them to this folder.

**Correction.** Section 0 of the sealed plan says the manuscript's noise section already contained a similar law. That statement referred to an inactive draft file, and it is corrected in `VALIDATION_REPORT.md`.

## Replay

```
python replay/replay_and_compare.py all T1
```

This recomputes every deterministic prediction in a scratch workspace and compares it with the sealed files. Roundoff-level diagnostic fields, such as finite-difference residuals and synthetic-sample moments, are reported but not held to the 1e-10 tolerance.
