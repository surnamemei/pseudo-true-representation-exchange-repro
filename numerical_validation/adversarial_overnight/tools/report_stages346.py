"""Write stage reports for Stages 3 (3A, and 3B if present), 4 and 6 (no new computation)."""
import csv
import json
from pathlib import Path

H = Path(__file__).resolve().parent.parent


def J(p):
    p = H / p
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else None


# ------------------------------------------------------------------ Stage 3
s3 = J("stage3_bwidth/DONE_3A.json")
s3b = J("stage3_bwidth/DONE_3B.json")
L = ["# Stage 3 — parameter width in b", "", "## 3A. Numerical scan (continuum tangent problem; diagnostic only, never 'certified')", "",
     f"Independent Gauss–Legendre (256-node) realified least-squares implementation, cross-checked at b = 10 against the frozen closed form: "
     f"max |ΔR| = {s3['crosscheck_b10']['max_abs_diff']:.1e}, |Δλ∞| = {s3['crosscheck_b10']['lambda_diff_vs_frozen']:.1e} (pass).", "",
     "| b | λ∞(b) | v_A | v_B | R_vv,A | R_vv,B | β_A | β_B | slope ∂λ(R_A−R_B) | coalescent gap | nearest regular competitor gap | tail margin | pass |",
     "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
for r in s3["continuum"]:
    L.append(f"| {r['b']:.1f} | {r['lam']:.8f} | {r['vA']:.4f} | {r['vB']:.4f} | {r['Rvv_A']:.3e} | {r['Rvv_B']:.3e} | {r['beta_A_re']:.4f} | "
             f"{r['beta_B_re']:.4f} | {r['slope']:.4f} | {r['coalescent_margin']:.3e} | {r['third_margin']:.3e} | {r['tail_margin']:.3e} | {r['pass']} |")
L += ["", f"**Largest symmetric interval supported by the frozen numerical criteria: b ∈ [{s3['diagnostic_interval'][0]:.1f}, {s3['diagnostic_interval'][1]:.1f}] "
      f"(radius {s3['diagnostic_symmetric_radius']}, the largest radius tested).** Evidence level GLOBAL_NUMERICAL (tangent problem); diagnostic, not certified.", "",
      "Secondary: finite N = 21, z = 2 record (full-torus crossing analysis at each b):", "",
      "| b | verdict | ε_c | J* | third margin (rel.) | min Hessian eig |", "|---|---|---|---|---|---|"]
for r in s3["finite_N21_z2"]:
    L.append(f"| {r['b']:.1f} | {r['verdict']} | {r.get('eps_c', float('nan')):.8f} | {r.get('J_star', float('nan')):.6f} | "
             f"{r.get('third_margin', float('nan')):.4f} ({r.get('third_margin_rel', float('nan')):.3f}) | {r.get('min_eig', float('nan')):.4f} |")
L += ["", "b = 10 is not a tuned point: across b ∈ [9.5, 10.5] the continuum crossing coefficient changes by < 1.6% and every margin stays positive; "
      "the certified finite z = 2 exchange persists at every b with ε_c ∈ [0.2476, 0.2503] and third-competitor margins ≥ 1.03·J*.", ""]
L += ["## 3B. Attempted explicit interval certificate (continuum problem, uniform in b)", ""]
if s3b:
    L += [f"Verdict: **{s3b['verdict']}** (runtime {s3b['runtime_minutes']:.1f} min of the {s3b['cap_minutes']:.0f}-min pre-registered cap; "
          f"{s3b['blocks_ok']}/{s3b['blocks_attempted']} attempted 0.01-wide blocks passed).", "",
          "| candidate b interval | blocks | primary (mpmath.iv) status | Arb replay |", "|---|---|---|---|"]
    for c in s3b["candidates"]:
        L.append(f"| [{c['interval'][0]}, {c['interval'][1]}] | {c['blocks']} | {c['status_primary']} | {c.get('arb_replay', '—')} |")
    if s3b["verdict"] != "CERTIFIED_B_INTERVAL":
        L += ["", "**Why no interval was certified (diagnosis).** Every one of the 24 completed blocks (b ∈ [9.88, 10.12]) passed all *local* obligations at all 100 sub-tiles: "
              "parametric Krawczyk inclusion (max contraction ≤ 0.67), positive reduced curvatures (min R_vv,A ≥ 4.7e-5), nonzero satellite coefficients (β_A ≥ 0.120, |β_B| ≥ 0.063), "
              "crossing slope ≥ 0.106, strict convexity of both ±0.1 windows, coalescent margin ≥ 8.3e-4 and tail margin ≥ 3.59e-3 (minima over all 24 blocks). "
              "Each block then exhausted the 2.5-million-visit far-field budget while sweeping v upward from −100: the sweep stopped at v ≈ −1.32. "
              "The flat region between well A and the central chart consumed the budget (≈ 5.7e5 leaves in v ∈ [−2, −1.32], 2.5e5 in [−3, −2]). "
              "The near-zero region and the whole B neighbourhood were never reached. "
              "The cause is dependency overestimation of naive mpmath.iv enclosures in b over flat regions: many cells had to be split down to single 1e-4 sub-tiles. "
              "It is not a certified counterexample: at every sampled b the numerical margins are positive (3A). "
              "The second-wave blocks (b beyond ±0.12) were abandoned at the 108-min primary deadline, and their workers were terminated after the verdict was written. "
              "A feasible future attempt would use Arb as the primary evaluator (measured ≈ 40× faster per cell here) and centred (mean-value) forms in b for the cost and gradient predicates."]
    if s3b.get("certificate_record"):
        cr = s3b["certificate_record"]
        L += ["", "Certificate record (widest certified candidate):", "", "```", json.dumps(cr, indent=1), "```"]
else:
    L += ["3B not yet complete when this report was generated."]
L += ["", "Method: b tiled into 0.01 blocks × 100 sub-tiles of width 1e-4; per sub-tile a parametric Krawczyk inclusion for (v_A, v_B, λ) with a mean-value form in b; "
      "strict convexity of R on ±0.1 windows around the root boxes (unique stationary point per window); full-line exclusion on [−100, 100] by cost/gradient/concavity/monotone predicates "
      "with adaptive (v, b)-bisection; per sub-tile coalescent and |v| ≥ 100 tail margins; exact rational coverage of every sub-tile's partition; independent python-flint/Arb (256-bit) replay "
      "recomputing every enclosure, preconditioner, λ-enclosure and incumbent from the stored geometry only. The earlier Stage-23 single-box attempt failed because its boxes were ~50× too wide for the Hessian enclosure."]
(H / "stage3_bwidth/STAGE3_REPORT.md").write_text("\n".join(L), encoding="utf-8")

# ------------------------------------------------------------------ Stage 4
s4 = J("stage4_imbalance/DONE.json")
L = ["# Stage 4 — strong-pair amplitude imbalance", "", f"Classification: **{s4['classification']}** (evidence GLOBAL_NUMERICAL).", "",
     "A₋ = 1 − η multiplies the tone at −z, A₊ = 1 + η the tone at +z (the one nearer the weak tone at b = 10); N = 21, z = 2, b = 10, φ = π, rectangular; ε ∈ [0, 1].", "",
     "| η | verdict | ε_c | J* | A pair | B pair | min Hessian eig (A, B) | ε_c·slope/J* | third margin (rel.) | min within-fit separation |",
     "|---|---|---|---|---|---|---|---|---|---|"]
for r in sorted(s4["rows"], key=lambda r: r["eta"]):
    L.append(f"| {r['eta']:+.2f} | {r['verdict']} | {r['eps_c']:.8f} | {r['J_star']:.6f} | ({r['A_u1']:.4f}, {r['A_u2']:.4f}) | ({r['B_u1']:.4f}, {r['B_u2']:.4f}) | "
             f"{r['A_eig']:.4f}, {r['B_eig']:.4f} | {r['eps_c'] * r['slope'] / r['J_star']:.3f} | {r['third_margin']:.4f} ({r['third_margin_rel']:.3f}) | {r['min_sep']:.3f} |")
L += ["", "Equal strong amplitudes are not special: across ±20% imbalance the separated exchange persists with ε_c ∈ [0.2369, 0.2482] and third-competitor margins ≥ 1.01·J*. "
      "The response is smooth and nearly symmetric in η (ε_c(−η) ≈ ε_c(η)). This is numerical evidence; the manuscript's existential local-stability corollary is the only proof-level statement."]
(H / "stage4_imbalance/STAGE4_REPORT.md").write_text("\n".join(L), encoding="utf-8")

# ------------------------------------------------------------------ Stage 6
s6 = J("stage6_noise/DONE.json")
ph6 = J("stage6_noise/POSTHOC_labels_and_widths.json")
summ = list(csv.DictReader((H / "stage6_noise/summary.csv").open()))
sa = s6["solver_invariance"]
L = ["# Stage 6 — noisy estimator: local vs global behaviour", "",
     f"N = 21, z = 2, b = 10, φ = π; ε = ε_c(1 + η), ε_c = 0.2481906301722774; SNR 20/30/40 dB with the frozen noise convention; **{s6['n_trials_per_setting']} records per (η, SNR)** "
     "(39,000 records; the pre-declared reduction rule was not triggered). Evidence level: GLOBAL_NUMERICAL (numerical global search per record; not certified).", "",
     "## 6A. Solver invariance (before trusting the Monte Carlo)", "",
     f"600 frozen records (η ∈ {{−0.1, −0.025, 0, 0.025, 0.1}} × 3 SNR × 40). Production (blind global search, no oracle seeds) vs reference (4× denser grid, 4× more random starts, oracle A/B seeds): "
     f"**{sa['prod_disagreements']}/{sa['records']} disagreements** in objective (> 1e-8) or label. Archived two-oracle-seed pipeline vs reference: **{sa['legacy_disagreements']}/{sa['records']}**. "
     f"Reference labels: A {sa['ref_label_counts']['A']}, B {sa['ref_label_counts']['B']}, OUT {sa['ref_label_counts']['OUT']}. No production upgrade was needed.", "",
     "## 6B. Selection probabilities and errors", "",
     "| SNR | η | P(A) | P(B) | P(out) | global MSE to noiseless winner | A-conditioned MSE (own root) | B-conditioned MSE (own root) | local sandwich trace A | local sandwich trace B |",
     "|---|---|---|---|---|---|---|---|---|---|"]
for r in summ:
    f = lambda k: (f"{float(r[k]):.4g}" if r.get(k) not in (None, "") else "—")
    L.append(f"| {r['snr']} | {float(r['eta']):+.3f} | {float(r['P_A']):.3f} | {float(r['P_B']):.3f} | {float(r['P_OUT']):.3f} | {f('MSE_to_winner')} | "
             f"{f('cond_MSE_own_A')} | {f('cond_MSE_own_B')} | {f('local_trace_A')} | {f('local_trace_B')} |")
h = s6["hypotheses"]
L += ["", "## 6C. Local theory vs Monte Carlo (pre-registered hypotheses)", "",
      "Local theory = frequency block of the misspecified sandwich covariance (σ²/2)H⁻¹Re(JᴴJ)H⁻¹ at each noiseless branch root (H includes residual curvature; FD-validated to "
      f"{s6['local_hessian_fd_max_rel_err']:.1e}). It is called a *local Hessian (sandwich) covariance approximation*, not a CRB/MCRB: the MCRB's regularity assumption (unique interior pseudo-true parameter) fails at η = 0 by construction.", "",
      f"* **H1 (branch-conditioned agreement away from the exchange; |η| ≥ 0.10, 30/40 dB, ≥ 100 selections): {h['H1']}** — trace ratios "
      + ", ".join(f"{x['branch']}@{x['snr']}dB,η={x['eta']:+.2f}: {x['ratio']:.3f}" for x in s6["H1_detail"]) + ".",
      f"* **H2 (global failure at the exchange; MSE to A ≥ 10 × local trace of A at η = 0): {h['H2']}** — ratios "
      + ", ".join(f"{x['snr']} dB: {x['ratio']:.0f}×" for x in s6["H2_detail"]) + ".",
      f"* **H3 (bimodality at η = 0; min(P_A, P_B) ≥ 0.2): {h['H3']}** — " + ", ".join(f"{x['snr']} dB: P_A={x['P_A']:.3f}, P_B={x['P_B']:.3f}" for x in s6["H3_detail"]) + ".",
      f"* Pre-registered 'meaningful TSP consequence' (H1 ∧ H2): **{h['TSP_consequence']}**.", "",
      "## Transition widths and comparison with the archived experiment", "",
      "| SNR | 10–90% probit width, pre-registered labels (OUT counts as not-A) | POST_HOC: nearest-root labels | archived (manuscript) |", "|---|---|---|---|"]
arch = {"20": "0.385 [0.366, 0.404]", "30": "0.129 [0.123, 0.136]", "40": "0.039 [0.0368, 0.0419]"}
for snr in ("20", "30", "40"):
    L.append(f"| {snr} dB | {s6['probit_widths'][snr]['width_10_90']:.4f} | {ph6['probit_width_nearest_root'][snr]['width_nearest_root_labels']:.4f} | {arch[snr]} |")
L += ["", f"POST_HOC_SENSITIVITY: all {ph6['n_OUT']} OUT records occur at 20 dB; {ph6['OUT_nearest_A_fraction'] * 100:.1f}% are nearest to A at distance "
      f"{ph6['OUT_dmin_quantiles'][0]:.2f}–{ph6['OUT_dmin_quantiles'][-1]:.2f} (median {ph6['OUT_dmin_quantiles'][2]:.2f}) — displaced A-family fits, not a third representation. "
      "The archived label rule (lower of two seeded refinements) assigns these to A.", "",
      "Findings: (i) the 30 and 40 dB widths reproduce the archived values; (ii) the 20 dB width is **not robust**: 0.479 (pre-registered labels), 0.418 (nearest-root labels) vs archived 0.385 "
      "with a model-conditional CI of ±0.02 that excludes both — the archived interval understates label-rule and η-design sensitivity at 20 dB; "
      "(iii) at 20 dB the local sandwich trace of A over-predicts the A-conditioned MSE (0.23–0.25 vs 0.33–0.40) because conditioning on the 1.0-radius label truncates the A-family tail — "
      "local theory is only reliable here at ≥ 30 dB.", "",
      "Figures: `stage6_noise/selection_probability.png|pdf`, `stage6_noise/bimodality.png`, `stage6_noise/global_vs_conditional_error.png|pdf`, `stage6_noise/local_theory_comparison.png|pdf`. "
      "Trial-level data: `stage6_noise/trials.csv` (39,000 rows)."]
(H / "stage6_noise/STAGE6_REPORT.md").write_text("\n".join(L), encoding="utf-8")
print("reports written")
