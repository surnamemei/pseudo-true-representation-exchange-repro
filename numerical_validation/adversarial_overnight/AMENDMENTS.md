# Protocol amendments (each recorded before the affected run)

## A1 — scale-consistent C3/C4 thresholds (recorded ≈2026-09-28 02:32 +10:00, before the first Stage-2 launch at 02:33:31 in overnight.log and before any Stage 2–8 computation)

**Reason (technical, not outcome-driven).** P0 criteria C3 (|dΔJ/dε| ≥ 1e-3) and C4 (third margin ≥ max(1e-4, 0.02·J*)) contain absolute floors calibrated on the z = 2 baseline (J* ≈ 0.81, E ≈ 84). In Stage 2 the objective scales as J ~ z⁴ (κ-family, ε = λz²) and J ~ z² (fixed-δ family, ε = μz); at z = 0.25 the whole A/B/third structure lives at J ~ 2e-4 and at z = 0.05 (fixed δ) at J ~ 1e-5. The absolute floors would then fail *mechanically*, independent of the phenomenon, and would turn scale into a spurious "DOES_NOT_PERSIST"/"INCONCLUSIVE". No Stage 2–8 result had been computed when this was noticed (only Stages 0–1, all far above both floors).

**Amended criteria (used from Stage 2 onward; Stages 0–1 re-evaluated, see below):**
* C3′: dimensionless transversality ε_c·|dΔJ/dε| / J* ≥ 1e-3 **and** |dΔJ/dε| above the numerical-noise floor 1e-10·E/ε_c.
* C4′: third-competitor margin ≥ max(0.02·J*, 1e-10·E), where E = Σw|y|² (numerical-noise floor instead of 1e-4).
* C1, C2, C5 and all other rules unchanged.
Both the original (C3, C4) and amended (C3′, C4′) booleans are stored for every switch, and every report states whether any verdict differs between the two rule sets.

**Re-evaluation of Stages 0–1 under A1:** all six rectangular/window/Hann² crossings and all eleven homotopy crossings have ε_c·|slope|/J* ≥ 2.5 and third margins ≥ 0.6·J* ≫ 1e-10·E, so every Stage 0–1 verdict is identical under both rule sets (checked programmatically in `stage1_window/amendment_A1_check.json`).

## A2 — code defect fix in Stage 2 evenness check (recorded 2026-09-28 02:37 +10:00)

`s2_phase.py` evaluated `(tangent_rel_diff or 1) < 1e-9`; Python treats an exactly-zero difference (κ = 0.10 gave 0.0) as false and substituted 1, failing M2A.3 spuriously. Fixed to substitute 1 only when the value is missing. Metric definition and thresholds unchanged. Stage 2 summary recomputed from the unchanged checkpointed scans/analyses.

## A3 — Stage 7 self-test harness fixes (recorded 2026-09-28 03:30 +10:00)

The Stage-7 directed-rounding self-test converted mpmath.iv endpoints through `mp.mpf` at the default 15-digit context and used the unsigned mantissa, so every exact-rational containment check reported a spurious failure (1600/1600). The harness now reads the raw endpoint tuples (sign, mantissa, exponent) exactly; with that fix all 1600 containment checks pass. This concerns the test harness only (no certificate, metric or threshold changed). Stage 7 summary regenerated from the unchanged checkpointed variant runs.

## A4 — Stage 3B implementation parameters (disclosure recorded 2026-09-28 05:25 +10:00, AFTER the first 3B block outcomes)

The parameters below were fixed in code before the 3B launch at 03:41:37 (file mtimes: `bcert_primary.py` 02:41:57, `bcert_block.py` 03:26:31, `s3b_bcert.py` 03:29:17), but this disclosure is written only after the first block results arrived; it changes nothing retroactively.
* Tiling differs from the pre-registered "base tile width 0.005, b-bisection depth ≤ 4": single-tile Krawczyk tests (03:05–03:20) showed parametric inclusion needs sub-tiles of width ≈ 1e-4 (contraction 0.27–0.55 there; failure at ≥ 2e-4 with the enclosure widths obtained). The implementation therefore uses 0.01-wide blocks, each split into 100 sub-tiles of width 1e-4, with adaptive (v, b)-bisection over the sub-tile index range (depth ≤ 7) for the full-line exclusion. This is finer, not coarser, than pre-registered.
* Strict-convexity windows ±0.1 around the hull of each block's root boxes, used to prove uniqueness of the stationary point near each root.
* Per-block far-field visit budget 2.5 × 10⁶ cells (unresolved cells ⇒ block not certified).
* Internal split of the pre-registered 150-min cap: 72% primary (deadline 05:29), 28% Arb replay.
* `bcert_arb.py` (replay only) was edited at 03:53 to check Krawczyk images against the primary's root boxes and to allow bounded bisection of replayed leaves. The replay never ran, because no candidate was primary-certified.
Outcome under these parameters: every completed block passed all local predicates (Krawczyk inclusion, curvature, coefficient signs, slope, convexity windows, coalescent and tail margins) but exhausted the far-field visit budget with ≈ 410 unresolved cells. The result is recorded as a budget failure, not as a counterexample.
