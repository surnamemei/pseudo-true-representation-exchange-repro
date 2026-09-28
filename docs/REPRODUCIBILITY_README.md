# Reproducibility guide — current manuscript

## Stage 24 final TSP revision after the adversarial validation (current)

**Active package.** Main source `paper/main_tsp_reviewfriendly.tex` (body in `paper/sections/body_reviewfriendly.tex`, with `continuum_stage14.tex`, `evidence_table.tex` = Table I, `crossN_table.tex` = Table II, `robustness_table.tex` = Table III, `noise_table.tex` = Table IV), PDF `paper/main_tsp_reviewfriendly.pdf` (10 IEEE two-column pages). Supplement `supplement/supplement.tex` (S1, S1b, S2–S8), PDF `supplement/supplement.pdf` (10 pages). Title: *Global Pseudo-True Representation Exchange in Under-Modelled Spectral Fitting*. Claim changes, evidence levels and file list: `paper/REVISION_NOTES.md`; skeptical assessment: `paper/FINAL_TSP_SELF_REVIEW.md`. Earlier sections of this guide remain valid for the certificate chains they describe, except where corrected below.

**Evidence sources.** Every new robustness, phase and noise number in the revision is read from the completed outputs in `validation/adversarial_overnight/` (`MORNING_REPORT.md`, `CLAIM_AUDIT.md`, `HOSTILE_REVIEW.md`, `RESULT_MATRIX.csv`, `AMENDMENTS.md`, `FREEZE.md`, `PREREGISTRATION.md`, `stage*/DONE.json`, `stage*/STAGE*_REPORT.md`). No scientific experiment was run during the revision. The only computations were figure regeneration from archived results and one closed-form consistency check (see `paper/REVISION_NOTES.md`, "Tests run").

**Rebuild the manuscript and supplement** (TeX Live; the project used WSL):

```
cd paper && latexmk -pdf -interaction=nonstopmode main_tsp_reviewfriendly.tex
cd supplement && latexmk -pdf -interaction=nonstopmode supplement.tex
python paper/validate_reviewfriendly.py
```

The validator is a static check: every input file exists, labels are unique, every `\ref`/`\eqref` and citation resolves, environments balance, no manual equation tags are used, the seven main figures exist and are referenced, Table I is included, and the main text points to each of S1–S8. It does not check wording; the evidence-label audit of the prose is recorded in `paper/REVISION_NOTES.md`. The final run reports `errors 0`.

**Regenerate all figures** from archived results only (no search or solver runs; the only evaluation is the closed-form tangent loss drawn in Fig. 2):

```
python paper/figure_scripts/make_revision_figures.py
```

| Printed figure | File | Source data |
|---|---|---|
| Fig. 1 | `paper/figures/fig1_continuation.pdf` | `branch_continuation.csv`, `stage2_reference.json` (certified crossing) |
| Fig. 2 | `paper/figures/fig2_tangent_landscape.pdf` | closed-form N=21 tangent loss at lambda_21 and lambda_21 +/- 0.002 (`stage3_tangent_reference.json`) |
| Fig. 3 | `paper/figures/fig3_largeN_fixed_geometry.pdf` | `lambdaN_asymptotics.csv`, `stage14_asymptotic_coefficients.json` |
| Fig. 4 | `paper/figures/fig4_z_law.pdf` | `z_continuation_vs_asymptotics.csv` |
| Fig. 5 | `paper/figures/fig5_robustness.pdf` | validation `stage1_window`, `stage3_bwidth`, `stage4_imbalance` |
| Fig. 6 | `paper/figures/fig6_phase_regimes.pdf` | validation `stage2_phase` |
| Fig. 7 | `paper/figures/fig7_noise_local_global.pdf` | validation `stage6_noise/summary.csv` |
| Fig. S1 | `supplement/figures/figS1_weak_phase_resolution.pdf` | copy of archived `paper/figures/fig7_phase_resolution.pdf` |
| Fig. S2 | `supplement/figures/figS2_window_branches.pdf` | copy of validation `stage1_window/window_branches.pdf` |
| Fig. S3 | `supplement/figures/figS3_posthoc_fixed_delta.pdf` | validation `stage2_phase/POSTHOC_fixed_delta_all_switches.json` (post hoc) |
| Fig. S4 | `supplement/figures/figS4_noise_distributions.pdf` | validation `stage6_noise/trials.csv`, `noiseless_roots.json`, `summary.csv` |

Older files in `paper/figures/` (`fig1_geometry`, `fig2_objective`, `fig3_critical_law`, `fig4_objective`, `fig4_persistence`, `fig5_noise_probability`, `fig6_bimodality`, `fig7_phase_resolution`, `fig_stage13_*`, `fig_stage17_*`) belong to earlier manuscript versions and frozen archives; the current manuscript does not include them.

**Adversarial validation (completed; not rerun for the revision).** From `validation/adversarial_overnight/`:

```
python -u run_all.py                    # stages 0,1,2,3,6,4,5,7,8; resumable; skips stages with DONE.json
python -u run_all.py --stages 0,1       # subset
bash run_all.sh                         # same, unbuffered, teeing to overnight.log
python tools/replay_sandbox.py prepare  # Stage-0 sandbox copy of the frozen project
python tools/replay_sandbox.py run      # rerun the seven certificate chains inside the sandbox
python tools/replay_sandbox.py compare  # field-by-field comparison with the frozen outputs
```

The preregistration (`PREREGISTRATION.md`, hash in `prereg_hash.json`), the four amendments A1–A4 with their timing (`AMENDMENTS.md`; A4 was documented only after the first results of the b-interval attempt), and the freeze (`FREEZE.md`, `freeze_manifest.json`, 766 SHA-256-hashed files) precede all reported numbers. Provenance limits are stated in Supplementary Sec. S8: the Stage-1 window protocol was frozen after earlier exploratory window results had been viewed, so it is confirmatory rather than fully blind; the preregistered Stage-2 phase classification is `INCONCLUSIVE`, and the fixed-phase scaling fits are post-hoc sensitivity analyses.

**Replay of every theorem-level certificate chain** (Stage 0; frozen scripts, unmodified, run with the sandbox as working directory; all exit codes 0; every output byte-identical or differing only in elapsed-time fields or JSON key order; `stage0_baseline/replay_comparison.json`):

| Chain | Commands | Time (s) |
|---|---|---:|
| R1 finite-spacing Arb replay | `stage12_full_arb_replay.py 21`, `stage12_full_arb_replay.py 31`, `stage12_cover_audit.py` | 119 |
| R2 five tangent instances | `stage12_generalN.py 11`, `15`, `21`, `31`, `41` | 2089 |
| R3 continuum certificate | `stage14_continuum.py`, `stage14_continuum_interval.py`, `stage14_continuum_global_hybrid.py`, `stage14_replay_global.py`, `stage14_continuum_endpoint_interval.py` | 154 |
| R4 independent Arb large-N | `stage16_independent_arb.py` | 45 |
| R5 N >= 10,001 floor | `stage17_uniform_trial.py 10001`, `stage17_refine_uniform.py 10001`, `stage17_largeN_margin_summary.py`, `stage17_write_largeN.py`, `stage17_annotate_condition_thresholds.py` | 697 |
| R6 finite-width epsilon interval | `stage17_epsilon_refine.py 0.035 14287`, `stage17_epsilon_root_tiling.py 0.035 4096`, `stage17_epsilon_inner_tiled.py _0p035 3000000`, `stage17_epsilon_joints.py _0p035`, `stage17_epsilon_coverage_audit.py 0.035`, `stage17_write_epsilon_certificate.py 0.035`, `stage17_epsilon_coeff_crosscheck.py` | 4391 |
| R7 lower-transition constants | `stage13_build.py` | 10 |

**Correction to the Stage 17 section below.** The finite-width chain consumes the complete outer refinement `..._0p035_14287.json`; the first command must be `python stage17_epsilon_refine.py 0.035 14287`, not `... 0.035 40` (a partial-limit invocation).

Runtimes are single-process on an AMD Ryzen 9 9950X3D (16 cores, 61 GB RAM), Windows 11, Python 3.13.9, mpmath 1.3.0, python-flint 0.9.0 (FLINT 3.6.0). An independent 50-digit joint solve of the N=21 stationarity and equal-cost equations reproduces the 110-digit crossing to 5.1e-51.

**Scope of the independent (Arb) replay.** A separate python-flint/Arb implementation that imports no primary evaluator (static import audit) replays every cell, inner and root obligation of the N=21 and N=31 crossing certificates (R1) and the 27,708-cell continuum partition for 0 <= N^-2 <= 35,377^-2 (R4). The sharper N >= 10,001 floor (R5) and the N=21 finite-width amplitude interval (R6) have been rerun with the primary directed evaluator only; they have no second-backend replay.

**Integrity.**
- *End of validation (2026-09-28 05:30).* `tools/integrity_check.py` recorded all 766 frozen files unchanged. That record is preserved byte-for-byte as `paper/freeze_records/integrity_check_end_of_validation.json`.
- *Final pre-submission check.* Command: `cd validation/adversarial_overnight && python tools/integrity_check.py` (exit 0). It rewrites `integrity_check.json`; copies are in `paper/freeze_records/`. Result: 766 checked, 732 unchanged, 34 changed, 0 missing, no new top-level files.
- *What changed.* All 34 changes are manuscript-side: main and supplement sources, figures regenerated from archived results, the manuscript checker, the bibliography, LaTeX build outputs, this guide, `paper/CURRENT_VERSION.md` and `final_submission_readiness.md`. No certificate script, certificate input or certificate output changed.
- *One correction outside the manifest.* During the revision the Stage-3B margin summary in `validation/adversarial_overnight/tools/report_stages346.py` was corrected from two-block to all-24-block minima: contraction <= 0.67, R_vv,A >= 4.7e-5, coalescent margin >= 8.3e-4, tail margin >= 3.59e-3. These values were checked against `stage3_bwidth/bcert/blocks/*.json`. The derived `STAGE3_REPORT.md`, `STAGE4_REPORT.md` and `STAGE6_REPORT.md` were regenerated by that deterministic script; only the Stage-3B text changed. No raw validation output was modified.
- *Records.* The complete inventory is `paper/FINAL_CHANGED_FILES_AUDIT.md`, and the freeze record is `paper/FINAL_SUBMISSION_FREEZE.md`. Per-file hashes are in `paper/stage24_freeze_diff.json` (frozen files) and `paper/stage24_revision_manifest.csv` (deliverables).

**Attempted explicit b-interval certificate.** A parameterized interval certificate of the continuum crossing uniform in b was attempted for b in [9.75, 10.25] and nested subintervals within a 150-minute cap (`stage3_bwidth`, `s3b_bcert.py`, `bcert_*.py`). In all 24 completed blocks (b in [9.88, 10.12]) every local obligation passed (contraction <= 0.67, R_vv,A >= 4.7e-5, coalescent margin >= 8.3e-4, tail margin >= 3.59e-3; minima over all blocks); each block then exhausted its 2.5-million-cell full-line exclusion budget in the flat region between branch A and the central chart. This is an efficiency limit, not a counterexample; no failed predicate was observed. The theorem-level statement remains the existential local b-neighbourhood; b in [9.5, 10.5] is numerical only.

## Stage 22 continuum b-dependence and final claim freeze

The active scientific manuscript is `paper/main_tsp_reviewfriendly.tex`; the active mathematical supplement is `supplement/supplement.tex`. The earlier `TSP_SUBMISSION_FREEZE_v1` ZIP remains an unchanged prior submission snapshot. Stage 22 adds one continuum weak-location diagnostic and a source audit, without new N families, finite-spacing bridge claims, or applications.

Run `python stage22_b_map.py` to regenerate `lambda_inf_vs_b.csv` and `lambda_inf_vs_b.pdf`. The script uses the existing phase-pi continuum tangent loss, with `b` external, and 65-decimal-digit mpmath stationary/equal-cost solves at integer `b=6,...,14`. Each row includes `lambda_inf(b)`, A/B satellite coordinates, common cost, positive curvatures, nonzero satellite amplitudes, crossing slope, coalescent gap, nearest sampled regular competitor margin, a conservative tail lower margin, and an explicit evidence status. It scans `v` on `[-100,100]` with 8,001 samples plus near-zero probes and local refinements; simple sinc derivative inequalities bound `|v|>=100`. **These non-b=10 globality checks are numerical, not directed interval certificates.** `b=10` retains the existing independently certified full-line result. Its strict margins and the implicit-function theorem give an existential open b neighborhood with the same global noncoalescent exchange, but no numerical radius is certified.

At `b≈13.9758640010`, a separate 110-digit numerical solve locates the A regular satellite's contact with the central/confluent chart (`v_A=0`). This is not an A/B coalescence. The `b=14` row is therefore flagged near coalescent, even though its sampled stationary root and coalescent cost gap are positive. The standalone figure connects only integer b=6 through 13 and shows b=14 separately. Interpretation and all limitations are in `b_dependence_analysis.md`.

`quartic_term_source_audit.md` verifies the Hermitian `2 Re sum(conj(r0)*r1)` term against active LaTeX, direct high-precision code, and the real/imaginary directed-interval implementation. The audit passes. `final_b_novelty_audit.md` records the targeted primary-literature collision search and conservative positioning. The active PDFs and `final_submission_readiness.md` give the final submission state. Scientific work stops here; only formatting, references, cover letter, and metadata remain.

## Stage 21 final package and claim boundary

The active manuscript is `paper/main_tsp_reviewfriendly.tex` and the mathematical supplement is `supplement/supplement.tex`; their PDFs are built in place. The earlier `TSP_SUBMISSION_FREEZE_v1` archive remains an unchanged historical submission state. The current novelty claim is a signal-specific global best lower-order representation switch, its small-spacing and finite-record laws, phase-scaling dependence, and certified finite-spacing examples. It does not claim novelty for Maxwell crossings, generic nonuniqueness, local minima, pseudo-true parameters, or noise-induced switching.

**Phase scope.** The strong pair is `2 cos(z s + delta)`, with `phi` reserved for the weak-tone phase. The existing full-domain finite-spacing certificates are phase aligned (`delta=0`). A strict-margin/implicit-function corollary extends the small-spacing global noncoalescent exchange to an *unspecified open interval* of `delta=kappa z` around `kappa=0` for each covered finite `N`, and to the continuum tangent problem. There is no explicit uniform `kappa` radius. For fixed nonzero `delta`, the first-order tangent class gives the natural `epsilon=Theta(z)` competition scale; its N21 crossing coefficient and displayed phase map are numerical only. See `strong_pair_phase_analysis.md`, `fixed_delta_scaling.md`, `kappa_scaled_phase_theorem.md`, and `kappa_crossing_map.csv`.

**Spacing continuation.** Run `python stage21_z_figure.py` to read the archived Stage-18 `p_bridge_numerical.csv` and regenerate `z_continuation_vs_asymptotics.csv` and `.pdf`. It performs no new continuation solve. The 44 sampled points span `z=0.001` through `z=2`; only the terminal `z=2` crossing and its N21 amplitude interval have independent full-domain certificates. Relative to the archived *numerical equal-cost continuation*, the maximum sampled absolute relative errors are 6.3383% for `lambda_21*z^2` and 0.2228% for `lambda_21*z^2+c_21*z^4`, both at `z=2`. This comparison supplies no connected global certificate.

**Finite-SNR model-order sanity check.** Run `python stage21_order_selection.py --trials 120`. It writes trial-level `order_selection_sanity.csv`, grouped `order_selection_sanity_summary.csv`, and a CDF figure `order_selection_sanity.pdf`. The fixed design is N21, z=2, `eta=(epsilon-epsilon_c)/epsilon_c` in {-0.1,0,0.1}, SNR 20/30/40 dB, and 120 independent records per setting (seed 20260928). The exact same sample-noise convention as the archived branch experiment is used: independent circular complex Gaussian noise with `sigma2=(||x||_2^2/N)*10^(-SNR_dB/10)`, so real and imaginary parts have variance `sigma2/2`. Both orders use continuous-frequency variable-projection nonlinear least squares with unrestricted complex amplitudes and several fixed initializations; every twentieth order-two fit is audited with a full-torus grid and refinements. These numerical optimizations are not interval-global certificates.

For order K, the criterion is `BIC_K=2*N*log(RSS_K/N)+3*K*log(2*N)` up to a common constant. It counts 2N real observations and 3K real mean parameters (two per complex amplitude plus one frequency); the one common unknown variance parameter cancels in differences. `DeltaBIC=BIC_3-BIC_2<0` selects order three. The experiment selected order three in all 1,080 tested records, including all 360 at 20 dB. Thus the present baseline is best presented as a deliberate/fixed-order approximation, not a typical BIC order miss. This check is illustrative, not a new scientific claim. See `order_selection_sanity.md` for the distribution and search limitations.

**Noisy A/B width calculation moved out of the mathematical supplement.** In the archived 400-trial-per-amplitude branch experiment, define `eta=(epsilon-epsilon_c)/epsilon_c`. At each SNR, 13 binomial counts `(k_j,n_j=400)` are fitted by maximum likelihood to `p_j=Phi(a+b*eta_j)`. The reported 10--90% width is `W=2*Phi^{-1}(0.9)/|b|`. The approximate SE uses the inverse expected binomial-probit information `I=sum_j n_j*varphi(a+b*eta_j)^2/[p_j*(1-p_j)]*(1,eta_j)^T*(1,eta_j)`, followed by the delta method for W. The model-conditional 95% intervals at 20, 30, and 40 dB are [0.36591,0.40368], [0.12310,0.13551], and [0.03684,0.04185]; widths are 0.38480, 0.12930, and 0.03935. These intervals exclude probit-model discrepancy and unobserved search failures. Recompute from archived `branch_probability.csv` and `transition_width.csv` with `python stage20_noise_width_uncertainty.py`; the standalone `supplement/S7_noise_uncertainty.tex` is retained as a source note but is not included in the six-page mathematical supplement.

**Bibliographic positioning.** The prior-work distinctions in the active main text use the verified records for Pakrooh--Scharf--Pezeshki (IEEE TSP 64(9):2345--2354, 2016, DOI `10.1109/TSP.2016.2521617`), James--Anderson--Williamson (IEEE TSP 43(4):817--821, 1995, DOI `10.1109/78.376834`), and Fortunati--Gini--Greco--Richmond (IEEE Signal Processing Magazine 34(6):142--157, 2017, DOI `10.1109/MSP.2017.2738017`). They address stochastic threshold/subspace-swap or general misspecified estimation, not the present deterministic globally optimal 3-to-2 switch.

## Stage 17 active certificate and submission set

The active theorem uses the conservative, margin-driven sufficient floor **every odd N >= 10,001** for b=10 and phase pi, in addition to the five separately certified finite lengths. The former floor 35,377 and its independent Arb replay remain archived below. The new floor is validated by a directed parameter-interval proof over 0 <= N^-2 <= 10,001^-2, including the four adapted full-line competitor cells. The exact kernel tail majorant is below 4.846e-46 at the floor; no preset 1e-30 target selects it. See `tightened_N0_derivation.md`, `tightened_N0_certificate.json`, and `largeN_condition_thresholds.csv`. The new 10,001 certificate uses the original directed evaluator with a different exact rational tail bound; an independent Arb replay of this new floor has not been performed.

The N21, z=2 finite-width amplitude switch is recorded in `finite_width_epsilon_certificate.md`, `.json`, and `epsilon_interval_obligations.csv`. The canonical certified closed interval is [0.213190630172278, 0.283190630172277], of width 0.069999999999999. For every amplitude in that interval, A and B are the only globally competitive fits, their roots continue without collision, and their cost difference has a unique transverse zero. The proof retains exactly quadratic amplitude dependence in each fixed-projector loss, then certifies outer exclusion, adaptive inner cells, root slabs and shared slab boundaries. The exact-decimal coverage audit checks every original parent and adjacent slab. This interval does not connect the small-z theorem to z=2.

`inner_cell_accounting.md` and `.csv` reconcile the 3,284 archived crossing inner visits with all 2,801 unique terminal Arb replay obligations over the N21/N31 cases. `structural_stability_corollary.md` proves only an existential local parameter neighborhood. `final_hostile_review.md` and `final_submission_blockers.md` give the final claim audit. The main manuscript and supplement are `paper/main_tsp_reviewfriendly.tex` and `supplement/supplement.tex` with PDFs in those directories.

`stage17_artifact_manifest.csv` records SHA-256 hashes and byte counts for the active proof logs, sources, reports, and two PDFs. Regenerate it after any edit with `python stage17_build_manifest.py`.

To replay the new large-N floor, run:

```powershell
python stage17_uniform_trial.py 10001
python stage17_refine_uniform.py 10001
python stage17_largeN_margin_summary.py
python stage17_write_largeN.py
python stage17_annotate_condition_thresholds.py
```

The first run intentionally reports four failed predicates on the unrevised frozen partition; the targeted refinement replaces them with six passing terminal leaves. The final writer requires the refined cover and positive predicate margins.

To replay the selected finite-width certificate at radius 0.035 with 4,096 root slabs, run:

```powershell
python stage17_epsilon_refine.py 0.035 14287   # corrected at Stage 24; "0.035 40" was a partial-limit invocation
python stage17_epsilon_root_tiling.py 0.035 4096
python stage17_epsilon_inner_tiled.py _0p035 3000000
python stage17_epsilon_joints.py _0p035
python stage17_epsilon_coverage_audit.py 0.035
python stage17_write_epsilon_certificate.py 0.035
python stage17_epsilon_coeff_crosscheck.py
```

The output `stage17_epsilon_coeff_crosscheck.json` compares independent coefficient and direct-evaluator paths. The interval logs can be large; the concise certificate and per-parent obligations CSV are the review entry points. Smaller completed radii 0.0025, 0.005, 0.01, and 0.02 remain separate records. The chosen endpoints are sufficient, not proved maximal.

The old width-0.005 generic-amplitude slab failure is an operational history item. Stage 17 replaced that enclosure with coefficient-wise quadratics and root continuation; it is not an unresolved scientific contradiction. The 5 direct anchors and all older experiments retain their stated scope. No new noise, phase, or record-length sweeps were run.

## Archived Stage 16 independent audit and submission set

At Stage 16 the anonymous package was a 9-page main manuscript and 5-page supplement. The independent python-flint/Arb replay in stage16_independent_arb.py verifies the frozen 27,708-cell continuum/all-large-N partition uniformly in t=1/N^2 from zero through 35377^-2: 27,708 passed, ten original cells refined with 26 subdivision nodes, zero failed. It checks the parameterized Krawczyk root box, both curvature and satellite-coefficient signs, transverse slope, coalescent gap, and full finite-period tail. See independent_N0_replay.md, .json, and .csv. **35,377 was a conservative sufficient Stage-16 certification threshold**; separate direct certificates already cover N=11,15,21,31,41.

stage16_symbolic_asymptotics.py independently expands the centered midpoint kernel and perturbs the joint stationary/equal-cost equations with SymPy. It confirms a2=-0.23156010864980659834... and a4=-0.28314327490967740462.... It also writes stage16_residual_scaling.csv. During this audit, two Python integer-division ingress paths in stage14_asymptotics.py were repaired before the long coefficient decimals and four numerical large-N validation rows were regenerated. The displayed 11-digit manuscript coefficients and five previously certified finite-N anchors are unchanged. See lambda_asymptotic_independent_check.md.

The frozen submission index is certificate_archive_index.md; the selected artifact SHA-256 list is artifact_manifest.csv. CODE_README.md gives the shortest independent replay. final_hostile_review.md and final_submission_blockers.md record the claim-strength review and verdict. The main manuscript now includes a short, explicitly numerical noise consequence drawn from existing archived trials; no new noise experiment was run.

For that archived noise experiment, each complex sample has circular Gaussian variance sigma^2=(||x||_2^2/N)10^(-SNR_dB/10); the real and imaginary noise parts each have variance sigma^2/2. The reported transition width is the 10%-to-90% A-selection interval in eta=(epsilon-epsilon_c)/epsilon_c. There are 400 trials per tested amplitude point. The noisy global optimizations are empirical audits, not interval proofs.

For a fresh WSL environment, create a virtual environment and install python-flint==0.9.0 as described in CODE_README.md. The original Stage 15 checker remains a separate backend. The Arb script takes the frozen candidate center and cell partition as inputs; it does not import the original evaluator.

## Active package and claim scope

- Main: `paper/main_tsp_reviewfriendly.tex`, `paper/sections/body_reviewfriendly.tex`, `paper/main_tsp_reviewfriendly.pdf`.
- Mathematical supplement: `supplement/supplement.tex`, included S1, S1b, S2–S5 sources, and `supplement/supplement.pdf`.
- Review response: `response_to_major_review_II.md`.
- Stage-16 scientific claim: the conditional finite-N criterion had five individually certified instances and a uniform interval certificate for every odd N>=35377. Stage 17 sharpens that sufficient floor to 10,001 and establishes a finite N21 amplitude interval. The small-spacing radius remains N-dependent and existential; there is no certified connected z-bridge.
- Pre-Stage-15 active sources/PDFs are preserved under `stage15_before_revision`. Other `main*.tex` variants and historical stage reports are not the active manuscript.

## Archived Stage 15 effective threshold and finite-N expansion

At Stage 15, the effective sufficient threshold was `N0_interval=35377`, verified uniformly for `t=1/N^2` in `[0,35377^-2]`. `N0_analytic` from stand-alone closed-form global-loss bounds was not obtained; the first odd N satisfying the selected analytic kernel remainder inequality was 35377, and its complete interval replay passed. This historical proof threshold was superseded by Stage 17. No complete direct interval sweep below 10,001 was attempted; the five frozen finite-N anchors remain separate.

`stage15_uniform_interval.py 35377` replays all 27,708 frozen competitor cells with a parameter interval, checks a joint Krawczyk inclusion, and verifies coalescent and expanding-period tail margins. `stage15_finalize_certificate.py` checks the result, safety-factor margins, frozen partition provenance, and writes `explicit_N0_certificate.json`. The derived kernel/moment/loss expansions are in `largeN_expansion_derivation.md`; analytic root-shift formulas are in `lambdaN_asymptotic_theorem.md`. `stage15_asymptotic_data.py` creates `lambdaN_largeN_data.csv`, `lambdaN_residual_scaling.csv`, and `largeN_theory_vs_data.pdf`. Rows N=101,201,501,1001 are validation only.

```powershell
python stage15_uniform_interval.py 35377
python stage15_finalize_certificate.py
python stage15_asymptotic_data.py
```

The complete proof and its scope limits are in `explicit_N0_derivation.md`. The N^-2 and N^-4 coefficient decimals are high-precision evaluations of analytic derivatives, not separately interval-enclosed. Main and supplement remain 9 and 5 pages.

## Stage 14 continuum baseline

Original finite-N and finite-z certificates remain frozen. Stage 14 added a normalized continuum tangent loss at b=10, phase pi, a directed interval joint crossing, a replayed 27,708-cell partition on [-100,100], a coalescent gap, and an analytic |v|>=100 tail bound. These prove full-line continuum globality at the specified crossing. Its initial transfer was existential; Stage 15 makes the threshold effective. The existing C8 exhaustion applies to the covered odd records. No uniform small-z radius is certified.

The proof objects are `stage14_continuum_local_interval.json`, `continuum_global_certificate.json`, `stage14_continuum_global_partition.csv`, `stage14_continuum_global_replay.json`, and `stage14_continuum_endpoint_interval.json`. The candidate solver output `stage14_continuum_candidate.json` is numerical, not a proof premise. Reproduce the proof objects in this order:

```powershell
python stage14_continuum.py
python stage14_continuum_interval.py
python stage14_continuum_global_hybrid.py
python stage14_replay_global.py
python stage14_continuum_endpoint_interval.py
```

`stage14_asymptotics.py` derives the N^-2 and N^-4 coefficients by analytic implicit differentiation and writes `stage14_asymptotic_coefficients.json`, `lambdaN_asymptotics.csv`, and `continuum_vs_finiteN_figure.pdf`. Its N=101,201,501,1001 rows are numerical validations only; the five earlier anchors retain individual interval certification. The displayed large-N coefficients are high-precision numerical evaluations of analytic formulas, not separate interval enclosures. Derivation and limits are in `lambdaN_asymptotic_fit.md` and `largeN_persistence_theorem.md`.

The active main and supplement compile with the existing WSL TeX environment:

```powershell
wsl.exe -e sh -lc 'cd /mnt/f/Research/paper && latexmk -pdf -interaction=nonstopmode -halt-on-error main_tsp_reviewfriendly.tex'
wsl.exe -e sh -lc 'cd /mnt/f/Research/supplement && latexmk -pdf -interaction=nonstopmode -halt-on-error supplement.tex'
```

## Stage 13 endpoint audit

The current manuscript additionally proves that for N=21, b=10, phase pi, the coalescent tangent representation is strictly improvable for every positive scaled omitted amplitude. The directed interval constants and threshold comparisons are in lower_transition_certificate.json; the derivation is in lower_transition_certificate.md. Regenerate the 27 high-precision stationary-continuation rows, categorical data and cost figure with:

    python stage13_build.py

The output tangent_cost_vs_lambda.csv is numerical away from the original certified A/B crossing and the analytic near-zero neighborhood. The diagram and response files state this boundary. A separate full-period interval root-isolation attempt at lambda=0.02 was stopped after more than 60,000 cells without a completed certificate; its partial log is not evidence of global ordering. No original N21/N31 finite-z certificate or Stage-12 tangent proof file was changed. The Stage-13 edit backup is in stage13_before_revision.

## 1. General-N tangent certificates

`stage12_generalN.py` recomputes a high-precision tangent crossing for its explicit N, then performs joint interval Krawczyk validation and complete scalar derivative isolation over a cover of the frequency circle. It uses finite trigonometric sums and deflated moment expansions near v=0. It never substitutes a sampled grid for an interval cover.

```powershell
python stage12_generalN.py 11
python stage12_generalN.py 15
python stage12_generalN.py 21
python stage12_generalN.py 31
python stage12_generalN.py 41
```

The files `stage12_generalN/N*_reference.json` are high-precision candidates; `N*_certificate.json` contains the interval joint-root check, all scalar stationary root/cost/curvature enclosures and the coalescent bound. `N*_partition.csv` records the full terminal frequency cover, and `N*_roots.csv` the isolated stationary points. `generalN_certified_instances.csv` collects the checked conditions, including full interval strings and conservative displayed margin bounds.

The implementation constructs a parameterized version of the existing tangent interval equations in a separate Python module. It sets N explicitly in all finite sums, moment formulas and sample loops. Reference roots and quartic coefficients are recomputed for that N. The imported two-variable derivative-jet class is generic arithmetic; its original module's N21 objective is not called by this tangent verifier.

The analytical premises and proofs are in:

- `generalN_conditions.md`: C1–C8, general-phase finite-sum formula, coalescent projection.
- `generalN_theorem.md`: conditional theorem and coefficient proof.
- `lemma_compactified_subspaces.md`: periodic quotient, projector continuity and closure attainment.
- `lemma_one_central_scaling.md`: explicit rank-five derivative, uniform singular-value/Taylor bound.
- `lemma_both_central_limit.md`: recurrence proof allowing arbitrary rates and divergent amplitudes.

C1–C7 are finite-dimensional interval checks. C8 follows analytically for this family with odd N>=3; it is not inferred from absence of numerical competitors. v=0 is excluded from the “other regular minimum” count and handled by the coalescent bound. Stage 14 extends the verified criterion to all sufficiently large odd N, with an existential threshold.

## 2. Complete second-backend replay

Original interval backend: mpmath.iv 1.3.0. N21 originally used 110-digit reference roots and 90/65/45 decimal digits for local/outer/inner work; N31 used 100-digit root refinement and 90/70/45 digits. The second backend uses **python-flint 0.9.0, FLINT 3.6.0, Arb at 320-bit midpoint precision**. Arb's radius arithmetic is native outward ball arithmetic; increasing midpoint precision alone does not remove interval dependency overestimation.

A workspace-local python-flint installation is in `stage12_deps`; the replay scripts prepend that directory to their import path. To reproduce in a new environment, install python-flint 0.9.0 (or reproduce with an explicitly reported later version). Python 3.13, NumPy, SciPy and mpmath are also used by the candidate-generation scripts.

```powershell
python -m pip install python-flint==0.9.0 --target stage12_deps
python stage12_full_arb_replay.py 21
python stage12_full_arb_replay.py 31
python stage12_cover_audit.py
python stage12_reports.py
```

`stage12_arb_backend.py` implements finite-sum objectives and analytic derivative jets over Arb. It does not call mpmath.iv. `stage12_full_arb_replay.py` reads archived **coordinates and parameter brackets**, then freshly computes energies, feasible incumbents, projector displacement, objective/gradient/Hessian enclosures, Krawczyk preconditioners and crossing slopes. Stored old inequality values are not validity inputs. The original 1e-95 amplitude padding is included, as is the 1e-12 physical-radian outer radius inflation.

Outputs:

- `stage12_replay/N21_predicates.csv`, `N31_predicates.csv`: all archived cell and local/crossing obligations.
- `N21_refinement_leaves.csv`, `N31_refinement_leaves.csv`: child coordinates and successful predicates for locally refined parents.
- `N21_summary.json`, `N31_summary.json`: backend versions, counts, margins and runtime.
- `coverage_audit.json`: exact-decimal sweep-line rectangle-union proofs for all four outer partitions and eight expected inner neighborhoods.
- `refinement_coverage_audit.json`: exact union checks for all 172 refined parent cells.
- `full_second_backend_replay.csv`, `arb_full_replay_comparison.csv`: identical complete combined per-obligation tables, provided under both requested names.
- `full_second_backend_replay_summary.md`, `arb_full_replay_N21.md`, `arb_full_replay_N31.md`: interpretation and aggregate results.

There are **65,337 outer cells, 2,801 inner cells and 86 root/crossing obligations**, totaling **68,224**. Every obligation passes. Arb needs one bisection for 128 N21 and 44 N31 inner parents, giving 344 successful child leaves. Another 21 parents use a different valid exclusion predicate. This is full independent arithmetic replay with local refinement, not identical enclosures or bitwise output equality. No contradictory certificate result or unresolved parent remains. Printing uses Arb balls that include decimal-output rounding; decisions use the unprinted full-precision objects.

The original proof partitions are not regenerated by the second backend. They are treated as candidate covers, whose coverage and all proof obligations are independently verified. This distinction keeps replay independent of the original floating adaptive search. It is not proof-assistant formalization or third-party review.

## 3. Amplitude continuation and interpretation

```powershell
python stage12_branch_continuation.py
python stage12_constrained_check.py
python stage12_reports.py
```

The first script traces A from epsilon=0 and B backward from its crossing to zero, then forward slightly beyond it. Steps are at most 0.0005. It uses analytic gradients and Hessians, exact linear amplitude elimination, and independent 65-digit stationary refinement of the reported landmarks. Outputs are `branch_continuation.csv`, `branch_continuation_figure.pdf` (also copied into the main figures), and `stage12_continuation_checks.json`. The constrained script adds a full-period numerical fixed-satellite comparison. `branch_continuation_interpretation.md` explains the branch identities and evidence scope.

B exists at zero omitted amplitude as a higher-loss minimum; its satellite moves from about 14.6313 to 12.4634 at equality, while b=10. This is a finite-record pseudo-true representation. The plotted trajectories do not certify globality over their full amplitude range. Existing finite-spacing global certificates retain their exact bracket/sample scope.

The current figure generator for the older phase/resolution and noise figures is `stage11_revision_figures.py`; `paper/make_figures.py` alone would restore an obsolete version of the bimodality annotation. The new Fig. 1 is generated by the continuation script above.

## 4. Frozen finite-z artifacts and historical metadata correction

The original certificates remain unchanged:

|Artifact|SHA-256|
|---|---|
|`validated_globality_certificate.json`|`A3D6C4E381A21A7FA75CCD362D5A6A8478D1C14D24B680D1222DA747B54AABB3`|
|`N31_global_certificate.json`|`E4E32DC3595BCAF18EE48DF8C95152187C3944D6A61E6980B6EFCC3C87D57CAF`|

N21 data: `interval_boxes.csv`, `stage2_inner_cells.csv`, `stage2_reference.json`, `stage2_crossing_interval.json`, `stage2_crossing_global.json`, local/large Krawczyk checks and `stage2_cell_audit.json`. N31 data: `N31_interval_cells.csv`, `N31_inner_interval_cells.csv` and its certificate JSON. N21 includes the crossing plus the separately certified epsilon_c±0.035 amplitude samples. They do not establish the intervening interval.

The historical N31 JSON has eight superseded **comparison metadata** fields inside `high_precision_numerical_root`:

`epsilon_leading`, `epsilon_through_quartic`, `remainder`, `remainder_over_d4`, `relative_error`, `remainder_after_quartic`, `remainder_after_quartic_over_d6`, `relative_error_quartic`.

The shared high-precision solver evaluated the N31 equations correctly, then read N21 tangent coefficients to populate these comparisons. No crossing/root/cost/Hessian/Gram/cell predicate reads them. The corrected N31 quartic prediction at z=2 is 0.24857675463874848, with 0.227313% discrepancy from the N31 crossing. Use the new general-N coefficient file for theory comparisons; do not silently edit the frozen JSON.

The earlier dependency audit (`N31_leakage_audit.md/.csv`, `stage11_N31_poison_test.json`, `ARCHIVE_NOTE_N31_CERTIFICATE_REVISED.md`) records explicit N31 seeds, sample-count assignment, center overrides before inner predicates, and freshly checked radius/preconditioner choices. Poisoning N21 coefficients changed only those eight metadata fields. Poisoning default N21 inner centers reproduced all 825 N31 inner rows. Earlier clean reruns in `stage10_replay` matched N21 artifacts and N31 parsed fields (only JSON key order differed). These provenance checks supplement, and do not replace, the current full Arb replay.

## 5. Existing experiments and boundaries

No new noise or application experiment was added. Existing phase, resolution, asymmetry and noise data remain in the archive. The numerical phase cusp is at phi/pi≈0.33928149; the central-chart collision is near 0.46495169. Representative phase roots have interval checks, but no complete phase-region global certificate. The wider amplitude-slab attempt is recorded in `finite_epsilon_certificate.json`; its failed enclosures remain a limitation. The small-spacing theorem still has no numerical positive radius or certified bridge to z=2.

Historical noise files include `noise_branch_trials.csv`, `branch_probability.csv`, `mixture_statistics.csv`, `noise_model_fit.csv`, and `transition_width.csv`. Phase/resolution files include `phase_crossing_map.csv`, `stage11_phase_validation.json`, and `resolution_sweep.csv`. These are numerical evidence unless explicitly covered by a stated interval result.

## 6. Stage 18: attempted connected N21 \(p\)-bridge

The Stage-18 ledger is `p_bridge_certificate.md` and `p_bridge_certificate.json`. `p_bridge_numerical.csv` contains 44 high-precision A/B crossing samples from \(p=10^{-6}\) to 4. `p_bridge_slabs.csv` records 37 local pilot/archived slabs plus one archived full-domain global slab; its status column must be read before using any row as proof. `p_bridge_competitor_margins.csv` separates numerical third-candidate margins from the rigorous archived exclusion margin. `stage18_p_bridge_numerical.py`, `stage18_slab_pilot.py`, and `stage18_write_certificate.py` regenerate these ledgers; the frozen Stage-4 files remain untouched.

The small-\(p\) theorem has an existential radius with no explicit positive overlap point. The archived global finite-end interval is \(z\in[1.9999997,2]\), or \(p\in[3.99999880000009,4]\). The local Krawczyk interval is wider, \(z\in[1.9999,2]\), but lacks full-domain global exclusion outside the narrower interval. Isolated global nodes \(z=1,1.5\) and positive numerical third-candidate margins do not fill the gap. The single bounded subdivision repair first produced passing local slabs of widths about \(10^{-9}\) at \(z=0.1\), \(10^{-7}\) at 0.5, \(10^{-6}\) at 1, and \(10^{-5}\) at 1.5 and 1.9. This exposes ordinary-coordinate interval dependency and computational cost; no structural event has been certified over the unsolved gap. Consequently the paper keeps the small-spacing and finite-spacing theorems separate.

The quartic code/manuscript conjugation check is in `quartic_conjugation_audit.md`; the targeted literature result is in `final_novelty_audit.md`. The main manuscript's general-phase expression now uses \(G_0\) consistently, the all-large-odd-\(N\) evidence rows say analytical plus computer-assisted, and the continuum theorem follows the finite-\(N\) criterion on which it depends. Detailed replay cell arithmetic remains in this guide and machine-readable logs rather than the main text.

The Stage-14 continuum compact cover has 27,708 terminal cells: 25,274 cost-excluded, 1,420 derivative-sign excluded, 1,012 curvature/endpoint excluded, and two A/B isolators. In the uniform large-odd-\(N\) transfer, 27,704 pass their original predicates; two pass after one bisection each and two pass stronger no-root tests. The conservative cost-cell and transferred compact margins are \(6.28\times10^{-9}\) and \(4.12\times10^{-10}\), respectively. These totals describe the archived complete predicate logs; they do not replace them.

## 7. Build and review

The manuscript and supplement are multi-file IEEEtran projects. Compile with an existing LaTeX installation:

```text
cd paper
latexmk -pdf -interaction=nonstopmode -halt-on-error main_tsp_reviewfriendly.tex
cd ../supplement
latexmk -pdf -interaction=nonstopmode -halt-on-error supplement.tex
```

`python paper/validate_reviewfriendly.py` checks included sources, labels, citations and figure pointers. PDF layout is reviewed after rendering all pages with `pdftoppm`; previews are under `stage12_pdf_preview`. The reference audit is in `verified_computing_related_work.md`. All generated proof/result artifacts are indexed and hashed in `stage12_artifact_manifest.json` after final compilation. Hashes establish file identity, not mathematical validity; the supplied predicates establish the latter.
