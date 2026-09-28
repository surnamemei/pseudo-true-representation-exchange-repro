# FREEZE — adversarial overnight invariance validation

Freeze time: **2026-09-28 02:17:56 +10:00** (local, AUSEST). Written before any new experiment in this run.
Companion files: `freeze_manifest.json` (SHA-256 of 766 project files, 1.407 GB, hashed in 4.3 s) and
`PREREGISTRATION.md` (all stage protocols, metrics and classification rules, frozen at the same time).
Any later protocol change must be logged in `AMENDMENTS.md` **before** the affected rerun, with its reason.

## 1. Version-control state

* `F:/Research` is **not a git work tree** (`git rev-parse` fails; the only `.git` directory under it belongs to the unrelated `GaussSt_TopoID/` project). There is therefore no HEAD commit.
  The SHA-256 manifest `freeze_manifest.json` is the provenance substitute: every top-level project file plus `paper/`, `supplement/`, `manuscript/`, `stage12_generalN/`, `stage12_replay/` is hashed (path, bytes, mtime, sha256).
* No `git init`, commit, push, upload or submission is performed by this run.
* "Dirty/untracked" analogue: the last packaged snapshot is `TSP_SUBMISSION_FREEZE_v1.zip` (2026-09-27 20:29, sha256 `edf1ad64…9e31c04`). **66 hashed files were modified after it** (Stages 20–23). They fall in two groups:
  * Recorded in `REPRODUCIBILITY_README.md` / `final_submission_readiness.md` (00:11–01:12): Stage 20 phase regimes (`stage20_*`, `kappa_crossing_map.csv`, `fixed_delta_*`, `kappa_scaled_phase_theorem.md`, `strong_pair_phase_analysis.md`), Stage 21 order selection and z-figure, Stage 22 b-map (`stage22_b_map.py`, `lambda_inf_vs_b.*`), audits, and rebuilt manuscript/supplement PDFs and sources.
  * **Unrecorded work in progress (01:45–01:53), not referenced by the README**: `stage23_b_interval.py` + `b_interval_certificate.json` (a *failed* directed b-interval Krawczyk attempt: every width 0.25…0.01 failed inclusion), `stage23_b_arb_replay.py` + `b_interval_arb_replay.json` (Arb replay of the same failed boxes), `stage23_estimator_analysis.py` + `estimator_mode_statistics.csv`, `local_covariance_validation.csv`, `estimator_analysis_validation.json`, `estimator_mode_mixture.pdf` (post-processing of archived noise trials), and `stage23_window_study.py` + `window_robustness_*` (a 72-second exploratory rectangular/Hann/Hamming window scan). A backup directory `stage23_before_revision/` suggests a manuscript revision was being prepared. **These files are treated as frozen prior exploratory evidence, not as results of this run.**
* **Disclosure (pre-registration integrity):** before freezing I read the Stage-23 exploratory outputs, including `window_robustness_summary.csv` (numerical crossings reported at ε≈0.1471 for Hann and ≈0.1836 for Hamming with weight normalization Σw=N), the Stage-20 κ map, and the Stage-22 b map. The protocols in `PREREGISTRATION.md` were therefore not written blind to those exploratory numbers. To limit the effect, the Stage 1 ε interval is fixed from the rectangular result plus a broad range (not from the Hann/Hamming numbers), and classification thresholds are fixed numerically in advance.

## 2. Current manuscript result claims (as frozen)

Active sources: `paper/main_tsp_reviewfriendly.tex` + `paper/sections/body_reviewfriendly.tex` (9 pages) and `supplement/supplement.tex` (6 pages). Claim ledger: `FINAL_CLAIM_LEDGER.md`. The manuscript is **not edited** by this run.

| ID | Claim (paraphrase of frozen text) | Frozen status |
|---|---|---|
| M1 | Theorem 1: any odd N satisfying C1–C8 has a global noncoalescent A/B exchange for small z with ε_c(z)=λ_N z²+c_N z⁴+O(z⁶); radius existential. | Analytical, conditional |
| M2 | Corollary: C1–C8 hold for N∈{11,15,21,31,41} at b=10, φ=π; N=21: λ₂₁=0.0659804111269914, c₂₁=−0.00101774929746980, v_A=−4.15344048664877, v_B=12.4680948443009. | Analytical + computer-assisted |
| M3 | Theorem 2: continuum crossing λ_∞=0.0665069703923955…, v_A=−4.06003021618…, v_B=12.42390445040…; all odd N≥10,001 satisfy C1–C8; λ_N=λ_∞−0.23156010865/N²−0.28314327491/N⁴+O(N⁻⁶). | Analytical + computer-assisted |
| M4 | Proposition 1: N=21 and N=31, z=2, b=10, φ=π: unique global A/B equal-cost point in the certified bracket; ε_c(21)=0.2481906301722774…, ε_c(31)=0.2491430898525018…. | Interval-certified global |
| M5 | Proposition 2: N=21, z=2: A uniquely global on [0.213190630172278, ε_c), tie at ε_c, B uniquely global on (ε_c, 0.283190630172277]. | Interval-certified global |
| M6 | Corollary (local parameter stability): sufficiently small changes in weak location, phase and strong-amplitude ratio preserve a separated transverse exchange; neighbourhood existential. | Deduction from strict margins |
| M7 | Scaled strong-pair phase δ=κz: z² law persists on an unspecified open κ-neighbourhood of 0; R_{N,κ}=R_{N,0}+4κ²P_N(v); fixed nonzero δ is a different first-order (ε=Θ(z)) regime; N21 fixed-δ coefficient ε/z≈0.6274474|sin δ| numerical. | Analytical + numerical |
| M8 | b-dependence: strict margins give an existential open b-neighbourhood of 10; integer b=6…13 numerical crossings with λ_∞(b) in ≈[0.0664,0.0870]; b≈13.9759 A-chart contact with central chart. | Existential + numerical |
| M9 | Resolution sweep z∈{2,2.5,π,3.5,4,5,2π}: numerical crossings, A/B lowest two minima, third-minimum margins 0.88296–10.5029. | Numerical |
| M10 | Noise consequence (N=21, z=2): 400 trials/point; P(A)=0.52 at crossing (20 dB); probit 10–90% widths 0.385/0.129/0.039 at 20/30/40 dB; "near the exchange, one local pseudo-true parameter does not summarize global estimator output". | Numerical, model-specific |
| M11 | Discussion: "the pseudo-true parameter can be nonunique at the exchange and switch discontinuously as an omitted component varies"; "local covariance about one fitted mode does not describe the other competitive mode". | Interpretive |
| M12 | BIC sanity check selected order three in all 1,080 records; the paper positions order two as deliberate fixed-order approximation. | Illustrative |

## 3. Exact baseline parameters (frozen)

* Record: odd N, t = −(N−1)/2,…,(N−1)/2, s_t = t/N. Baseline **N = 21, z = 2, b = 10, φ = π**, equal unit strong amplitudes, strong pair phase aligned at the record centre (δ = 0).
* Signal: x_t = 2cos(z s_t) + ε e^{iφ} e^{i b s_t} = e^{−izs_t} + e^{+izs_t} − ε e^{ibs_t} (φ=π).
* Physical frequencies: ω₁,₂ = ∓z/N, ω₃ = b/N. Normalised fitted coordinates u_j = N ν_j, torus u ∈ [−πN, πN).
* Certified crossing (110-digit reference, `stage2_reference.json`): ε_c = 0.248190630172277398203609851355…; A = (−4.772447353918614298…, 0.605170320276861964…); B = (−0.057653561926679808…, 12.463419356801659806…); J_A = J_B = 0.811807871917534443556…; dΔJ/dε = 8.259051045764865….
* Certified finite-width interval: [0.213190630172278, 0.283190630172277].
* Continuum (N→∞) anchor: λ_∞ = 0.0665069703923955354559…, v_A = −4.060030216180106…, v_B = 12.42390445039620…, R* = 0.00288787652883089…; curvatures R_vv,A ∈ [1.33681896603e-4, …], R_vv,B ≈ 9.6036910554e-4; slope ≈ 0.107817513650; coalescent margin ≥ 8.5116817e-4; tail margin ≥ 3.7492892e-3; regular cost-cell margin ≥ 6.28e-9.

## 4. Objective definitions (frozen)

* Fit class: x̂_t = B₁e^{iν₁t} + B₂e^{iν₂t}, B_j ∈ ℂ unrestricted, ν_j ∈ ℝ/2πℤ, unordered, no grid, no penalty.
* Objective (rectangular, existing): J(u₁,u₂;x) = min_B ‖x − VB‖² = ‖(I−P_V)x‖², V = [e^{iu₁s}, e^{iu₂s}]; nonconfluent closed form Eq. (review-J) with the centred Dirichlet kernel D_N(q) = sin(q/2)/sin(q/(2N)).
* Coalescent closure: the limit subspace span{e^{ims}cos(hs), e^{ims} s·sinc(hs)} (mean m, half-separation h); at h = 0 this is span{e^{ims}, s e^{ims}}.
* Tangent (small-z) loss: R_N(v,λ) = min_{α,β∈ℂ, κ∈ℝ} Σ_t |q_t − α − 2iκs_t − βe^{ivs_t}|², q_t = −s_t² + λe^{iφ}e^{ibs_t}; coalescent class {complex 1, complex s, real s²}.
* Continuum loss R_∞: replace D_N/N by d(q) = 2 sin(q/2)/q and S_k/N by μ₂ = 1/12, μ₄ = 1/80, μ₆ = 1/448.
* New for this run (pre-registered in `PREREGISTRATION.md`): weighted objective J_w = min_B Σ_t w_t |x_t − (VB)_t|², Σ_t w_t = N.

## 5. Existing numerical search settings, seeds and precision (frozen, for reference)

* Stage-5 core (`three_to_two_tone_stage5_core.py`): BFGS on (u₁,u₂), gtol 2e-10, maxiter 250; confluent substitution when |u₂−u₁| < 1e-6; grid 128×128 (grid_minima), top 16; FD Hessian step 2e-4; dedup radius 0.02.
* Resolution sweep (manuscript): 512×512 full-torus grid, ≤200 grid-minimum refinements, 32 random starts, strong-pair initialisation, multiprecision stationary solves.
* Archived noise experiment (`three_to_two_tone_stage5_noise.py`): RNG `default_rng(20260926)`, 400 trials per η, SNR 20/30/40 dB; **each noisy record refined only from the two noiseless branch roots A₀, B₀**; every 20th trial audited with a 192-grid search (top 24). Branch label = lower of the two seeded refinements, overridden by the audit when it finds a lower minimum.
* Stage-21 order selection: seed 20260928, 120 records per setting.
* Stage-23 window study: grid 384 (crossing grid 768), 12 random starts (seed 1207), 45 grid seeds, ε ∈ [0.05, 0.8] with 31 points.
* Interval precision: N21 original mpmath.iv 110-digit roots, 90/65/45 digits local/outer/inner; N31 100-digit roots, 90/70/45; continuum iv.dps 65–85; Arb replays python-flint 0.9.0 / FLINT 3.6.0 at 320-bit midpoint precision.

## 6. Existing branch-assignment rules (frozen)

* A: the branch continuing from the exact strong-pair fit (−2, 2) at ε = 0; B: the branch continuing backward from the crossing to (−0.0200738, 14.6313006) at ε = 0.
* Archived noise labels: A if the A₀-seeded refinement has J ≤ the B₀-seeded one, else B; audit overrides to "other" if the audit minimum is > 0.05 from both.
* Stage-23 estimator analysis: pseudo-true reference = A for η ≤ 0, B for η > 0.

## 7. Existing noise model (frozen; reused exactly)

Independent circular complex Gaussian noise w_t with σ² = (‖x‖₂²/N)·10^{−SNR_dB/10}; Re w_t and Im w_t each ~ N(0, σ²/2), independent across t; x is the noiseless record at the tested ε. SNR ∈ {20, 30, 40} dB. η = (ε − ε_c)/ε_c with ε_c = 0.2481906301722774.

## 8. Existing fitting / search pipeline (frozen)

Variable projection over the unordered frequency pair; complex amplitudes eliminated by linear least squares (numpy `lstsq`, MKL); local BFGS refinement from grid minima + structured and random starts; high-precision (mpmath) stationary/crossing solves for landmarks; interval certificates (mpmath.iv primary, Arb replay) for theorem-level statements. Grid-based global statements are numerical only.

## 9. Software, precision and hardware at freeze

Python 3.13.9 (Anaconda, MSC v.1929, win64); numpy 2.3.5 (BLAS/LAPACK: MKL `mkl-sdl`); scipy 1.16.3; mpmath 1.3.0; matplotlib 3.10.6; sympy 1.14.0; numba 0.62.1; tqdm 4.67.1; python-flint 0.9.0 (workspace-local `stage12_deps/`, FLINT 3.6.0). Windows 11 Pro 10.0.26100; AMD Ryzen 9 9950X3D (16 cores / 32 threads); 61.4 GB RAM (27.6 GB free at freeze); NVIDIA RTX 5090 present but already busy with desktop workloads (≈20 GB used) — **this run uses CPU only** (GPU nondeterminism and contention avoided). Worker processes run with `OMP_NUM_THREADS=MKL_NUM_THREADS=OPENBLAS_NUM_THREADS=1` and `MKL_CBWR=COMPATIBLE` for reproducibility. WSL Ubuntu-24.04 is available but not required.

## 10. Rules adopted for this run

1. Frozen baseline files are never modified. Certificate replays run in a sandbox copy (`stage0_baseline/replay_sandbox/`) and are compared to the frozen outputs; a post-run hash check against `freeze_manifest.json` verifies that no frozen file changed.
2. Evidence levels: `CERTIFIED`, `GLOBAL_NUMERICAL`, `LOCAL_NUMERICAL`, `ILLUSTRATIVE`, `FAILED`, `INCONCLUSIVE` — exactly one per result.
3. Metrics and classification thresholds are fixed in `PREREGISTRATION.md`; anything examined later is tagged `POST_HOC_SENSITIVITY` and never replaces a failed primary metric.
4. Negative results are reported as found; no re-searching for parameters that restore the phenomenon.
5. Existing certificate code is not edited; new code lives under `validation/adversarial_overnight/`.
6. No manuscript edits, commits, pushes, uploads or submissions.
