# Replay commands

`replay/replay_and_compare.py` runs all of the commands below. Use this page to rerun a single step by hand. Work inside a prepared workspace:

```
python replay/replay_and_compare.py prepare
cd _replay_work
```

The frozen scripts read and write files in the current directory. **Never run them inside `certificates/`**: that would overwrite the archived records.

Run the commands of each chain in the order shown. R4 and R5 read outputs of R3, so run R3 first in the same workspace.

| Chain | Commands | Establishes | Implementation |
|---|---|---|---|
| R1 | `python stage12_full_arb_replay.py 21`<br>`python stage12_full_arb_replay.py 31`<br>`python stage12_cover_audit.py` | Re-evaluation of every cell, inner and root obligation of the N = 21 and N = 31 z = 2 certificates: 37,753 + 27,584 outer cells, 2,801 inner cells and 86 root or crossing obligations. Exact-decimal audits of the outer and inner covers. | Independent Arb (python-flint 0.9.0, 320-bit) |
| R2 | `python stage12_generalN.py 11`<br>`python stage12_generalN.py 15`<br>`python stage12_generalN.py 21`<br>`python stage12_generalN.py 31`<br>`python stage12_generalN.py 41` | Criterion C1–C7 at each length: a joint Krawczyk root, complete periodic isolation, the coalescent bound, and λ_N and c_N enclosures | Primary (mpmath.iv) |
| R3 | `python stage14_continuum.py`<br>`python stage14_continuum_interval.py`<br>`python stage14_continuum_global_hybrid.py`<br>`python stage14_replay_global.py`<br>`python stage14_continuum_endpoint_interval.py` | Continuum crossing: local interval root, 27,708-cell partition of [−100, 100], replay, endpoint | Primary |
| R4 | `python stage16_independent_arb.py` | The same partition re-checked uniformly for 0 ≤ N⁻² ≤ 35,377⁻² (includes the continuum) | Independent Arb |
| R5 | `python stage17_uniform_trial.py 10001`<br>`python stage17_refine_uniform.py 10001`<br>`python stage17_largeN_margin_summary.py`<br>`python stage17_write_largeN.py`<br>`python stage17_annotate_condition_thresholds.py` | Parameter-interval certificate for every odd N ≥ 10,001 at fixed normalized geometry. The first step deliberately reports four failed predicates on the unrefined partition; the refinement replaces them with six passing leaves. | Primary |
| R6 | `python stage17_epsilon_refine.py 0.035 14287`<br>`python stage17_epsilon_root_tiling.py 0.035 4096`<br>`python stage17_epsilon_inner_tiled.py _0p035 3000000`<br>`python stage17_epsilon_joints.py _0p035`<br>`python stage17_epsilon_coverage_audit.py 0.035`<br>`python stage17_write_epsilon_certificate.py 0.035`<br>`python stage17_epsilon_coeff_crosscheck.py` | N = 21 amplitude interval [0.213190630172278, 0.283190630172277] | Primary |
| R7 | `python stage13_build.py` | N = 21 lower-endpoint constants | Primary |
| T1 | `python paper/transition_theory/code/transition_theory.py`<br>`python paper/transition_theory/code/second_order_diagnostic.py` | Deterministic part of the sealed first-order transition theory: the branch quantities; the slope g, checked by envelope and finite differences in float64 and 50-digit mpmath; and every predicted P(A) and global MSE at the 39 Monte Carlo settings | Numerical approximation, not a certificate |

## Notes

- **R6 refinement argument.** The second argument of `stage17_epsilon_refine.py` is the number of outer parent cells to refine; all 14,287 are required. The coverage audit and the certificate writer read `stage17_epsilon_outer_refinement_0p035_14287.json` and assert `parents == 14287`. A partial invocation such as `0.035 40` writes `..._0p035_40.json`, and the chain cannot complete.
- **Arb scripts and `stage12_deps/`.** The Arb scripts insert `stage12_deps/` at the front of the import path. Keep that directory absent, so that python-flint comes from the environment, or install python-flint 0.9.0 into it for the current platform: `pip install python-flint==0.9.0 --target stage12_deps`.
- **T1 and the validation.** T1 writes `results/transition_theory/deterministic/` inside the workspace. Its roundoff-level diagnostic fields (finite-difference residuals, synthetic-sample moments) are reported but not held to the 1e-10 tolerance. The one-shot comparison with the Monte Carlo (`validate.py`) is not a chain: it refuses to run again, and its sealed outputs are recorded in `paper/transition_theory/RESULT_RECORD.sha256`.
- **Figures.** Regenerate them from archived results with `python paper/figure_scripts/make_revision_figures.py` in the workspace, or run `python replay/replay_and_compare.py figures`.
