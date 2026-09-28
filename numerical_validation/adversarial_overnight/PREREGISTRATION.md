# PREREGISTRATION — stage protocols, primary metrics and classification rules

Frozen together with `FREEZE.md` (2026-09-28 ~02:25 +10:00), before any new experiment. Hash recorded in `STATUS.json` at the first runner start. Changes after that go to `AMENDMENTS.md` with a reason, **before** the affected rerun. Alternative metrics examined after results are tagged `POST_HOC_SENSITIVITY`.

## P0. Common numerical machinery (all finite-record stages)

**Record/signal.** t = −(N−1)/2…(N−1)/2, s_t = t/N,
x_t = A₋ e^{−i(z s_t+δ)} + A₊ e^{+i(z s_t+δ)} + ε e^{iφ} e^{i b s_t}. Baseline A± = 1, δ = 0, φ = π, b = 10, N = 21, z = 2.

**Objective.** J_w(pair) = min_{B∈ℂ²} Σ_t w_t |y_t − Σ_j B_j e^{iu_j s_t}|², symmetric weights w, Σ_t w_t = N (rectangular w ≡ 1 reproduces the frozen objective exactly).

**Global coordinates.** Unordered pair ↔ (m, h): u₁ = m − h, u₂ = m + h. The fitted subspace is spanned by c = e^{ims}cos(hs), d = e^{ims} s·sinc(hs); for symmetric w these are w-orthogonal, so
J(m,h) = Σw|y|² − |⟨c,y⟩_w|²/‖c‖²_w − |⟨d,y⟩_w|²/‖d‖²_w,
which is smooth on the whole (m,h) plane, even in h, 2πN-periodic in m, and equals the coalescent closure at h = 0. Analytic gradient and Hessian are used.

**Global search protocol GS(level).** For a record y:
1. (m,h) grid: m uniform on [−πN, πN) with M_m points, h uniform on [0, πN/2] with M_h points (h = 0 included). Grid local minima (3×3, periodic in m) sorted by J; keep K_grid lowest.
2. Structured seeds: strong pair (m,h) = (0, z) and the generating pairs {−z,z}, {−z,b}, {z,b}; central+satellite pairs {0, v_j} for n_sat satellites v_j uniform on the torus; continuation seeds (≤ 8 minima from neighbouring scan points, when scanning).
3. n_rand deterministic pseudo-random (m,h) starts, RNG = `numpy.random.default_rng(SeedSequence([20260928, stage, config_id, point_id]))`.
4. Batch damped-Newton refinement in (m,h) (eigenvalue-shifted Hessian, backtracking), stop at ‖∇J‖ ≤ tol_g or 400 iterations.
5. Keep converged points whose Hessian has min eigenvalue ≥ −1e-9·max(1,J); deduplicate by unordered-pair torus distance < 1e-3; classify **regular** if within-fit separation |u₂−u₁| (torus) ≥ 1.0, otherwise **coalescent-class**.

| level | M_m | M_h | K_grid | n_sat | n_rand | tol_g |
|---|---|---|---|---|---|---|
| production (P) | 1024 | 257 | 64 | 48 | 32 | 1e-10 |
| reference (R) | 2048 | 513 | 200 | 96 | 128 | 1e-11 |

(M_m, M_h scale with N/21 so the grid spacing ≈ 0.129 in u-units is kept for other N.)

**Crossing analysis CA (per configuration).** ε-scan on a frozen grid of 101 equispaced points on [0, ε_max] (stage-specific scaling below); GS(P) at each point. A *switch* is a pair of consecutive points whose best minima are > 1.0 apart. For each switch, the two branches P (below) and Q (above) are locally continued (warm-started Newton) and the root of ΔJ = J_P − J_Q is solved by Brent (xtol 1e-13); a branch whose continuation jumps > 0.5 is a **fold** (not an exchange). At ε_c: GS(R) validation, Hessians, amplitudes, slope by central difference (±1e-6).

**Qualifying separated regular crossing (used everywhere):** C1 P,Q distinct with pair distance ≥ 1.0; C2 both regular (separation ≥ 1.0, Hessian min eigenvalue > 1e-8, max|B_j| < 1e3); C3 |dΔJ/dε| ≥ 1e-3; C4 third-competitor margin (lowest other distinct minimum, coalescent class included, minus J*) ≥ max(1e-4, 0.02·J*); C5 GS(R) at ε_c finds no minimum below J* − 1e-9·max(1,J*) other than P,Q and recovers both.
* C1–C5 hold → crossing **present** (evidence `GLOBAL_NUMERICAL`).
* C1–C3 hold but C4 or C5 fails → **INCONCLUSIVE** (if C5 finds a strictly lower minimum, the switch is recorded as *not global*).
* No switch satisfying C1–C3 on the range → **absent**.
**Boundary-hit rule:** if no qualifying crossing is found on [0, ε_max], extend once to [ε_max, 3ε_max] (101 points); never further. "Approaches coalescence": within-fit separation of a winner < 1.0, or max|B| ≥ 1e3, or coalescent-class margin < 0.01·J*.

## Stage 0 — baseline reproduction (priority 1)

Numerical (N=21, z=2, b=10, φ=π, rectangular):
* B0.1 float64 CA crossing within |Δε_c| ≤ 1e-10 of 0.2481906301722774.
* B0.2 independent mpmath (50 digits) joint solve of ∇J_A = ∇J_B = 0, J_A = J_B: |Δε_c| ≤ 1e-30 vs the 110-digit reference; A, B roots within 1e-30; J* within 1e-30 of 0.8118078719175344435….
* B0.3 Hessians of A and B positive definite.
* B0.4 GS(R) at ε_c: A and B are the two lowest minima; third margin within 1e-6 of the archived 0.88296 (rounding tolerance 1e-5).
* B0.5 ordering: J_A < J_B at ε_c − 1e-3 and J_A > J_B at ε_c + 1e-3 (also at the certified interval ends 0.213190630172278 and 0.283190630172277).
Certificate replays (sandbox copies, original scripts unmodified; frozen outputs compared field-by-field): R1 `stage12_full_arb_replay.py 21/31` + `stage12_cover_audit.py`; R2 `stage12_generalN.py 11,15,21,31,41`; R3 Stage-14 continuum chain (`stage14_continuum.py`, `_interval.py`, `_global_hybrid.py`, `stage14_replay_global.py`, `_endpoint_interval.py`); R4 `stage16_independent_arb.py`; R5 Stage-17 N≥10,001 floor chain; R6 Stage-17 finite-width ε chain at radius 0.035; R7 `stage13_build.py` lower-transition constants. Pass = same verified/pass flags and all printed margins identical to the stated digits.
**Rule:** if B0.1–B0.5 fail, or any replay contradicts its frozen certificate, the run is `NOT_READY`; new scientific stages then run only as diagnostics.

## Stage 1 — window / metric invariance (priority 2)

Fixed: N=21, z=2, b=10, φ=π, A±=1, δ=0, order two. Primary windows (symmetric, Σw=N as weights in J_w):
rectangular; Hann = `numpy.hanning(21)` (zero end samples); Hamming = `numpy.hamming(21)`; DPSS = `scipy.signal.windows.dpss(21, NW=2.5)` first taper (sign made positive).
ε-range: **[0, 1.0]** (= 4× the rectangular ε_c; fixed from the rectangular result, not from exploratory window numbers), 101 points; boundary rule P0.
Per window: CA; record best/second/third minima, frequencies, amplitudes, Hessian eigenvalues, gaps at every scan point; crossing ε_c, third margin, coalescence proximity.
Per-window classification: `PERSISTS` (qualifying crossing present), `DOES_NOT_PERSIST` (absent after boundary rule), `INCONCLUSIVE`.
Secondary (pre-registered, not used for the primary verdict): (a) **Hann²** weights (taper applied to data and model; w = hann², Σw=N); (b) **rectangular→Hann homotopy** w_θ ∝ (1−θ)·1 + θ·hann, θ ∈ {0,0.1,…,1}: local continuation of the A/B crossing from θ=0 with GS(P) validation at each θ, to test whether the Hann exchange is the continuation of the rectangular one.
Outputs: `window_summary.csv`, `window_minima.npz`, `window_report.md`, plots.

## Stage 2 — record-origin / phase invariance (priority 3)

Strong pair 2cos(zs+δ); weak tone fixed at φ = π unless stated.
**2A (δ = κz).** κ ∈ {0, ±0.05, ±0.10, ±0.20, ±0.30, ±0.40, ±0.50}; diagnostic point ±0.3882278594 (frozen central-collision value, not a grid redefinition).
(i) Tangent level (N=21), independent implementation by direct realified least squares (not the closed forms): λ₂₁(κ), v_A, v_B, curvatures, coalescent margin, third regular-minimum margin (full-period scan, 8192 points + refinement).
(ii) Finite records N=21 at z ∈ {0.25, 0.5, 1.0, 2.0}: CA with ε = λz², λ-grid [0, 1.0] (101 points), boundary rule P0.
Metrics: M2A.1 crossing present at each (κ, z); M2A.2 relative deviation |ε_c/z² − λ₂₁(κ)|/λ₂₁(κ) at z = 0.25, 0.5, 1 (must decrease as z↓ and be < 5% at z = 0.25); M2A.3 evenness |ε_c(κ) − ε_c(−κ)|/ε_c(κ) < 1e-6 (finite) and < 1e-9 (tangent); M2A.4 quadratic shift: λ₂₁(κ) − λ₂₁(0) vs the frozen prediction 1.92838866512018·κ², relative error < 10% at κ = 0.05.
**2B (fixed δ).** δ/π ∈ {0.05, 0.10, 0.20, 0.30, 0.40}; N=21; fitting window z ∈ {0.05, 0.1, 0.2, 0.4}; out-of-window point z = 0.8 (reported only). CA with ε = μz, μ-grid [0, 1.0] (101 points), boundary rule P0. The analysed crossing is the qualifying crossing whose satellites are closest to the tangent-predicted pair (v≈1.618, v≈10.305); others are reported.
Fits on the 4 in-window points: log–log OLS exponent α; one-parameter least squares ε_c = C₁z and ε_c = C₂z² (relative residual RMS); secondary: drop smallest z, drop largest z.
**2C.** Table κ ↔ record-centre displacement τ = κN samples for N ∈ {21, 31, 101, 1001}; fraction of record = κ. Secondary (pre-registered): a *true* record shift also rotates the weak tone by bκ; tangent-level N=21 check of δ = κz **with** φ = π + bκ on the same κ grid.
**Stage-2 classification (z² law only):**
* `QUADRATIC_LAW_LOCALLY_STABLE_BUT_ALIGNMENT_SPECIFIC` if M2A.1–M2A.4 hold for all |κ| ≤ 0.20 **and** for every fixed δ the exponent α ∈ [0.8, 1.2] with the C₁z fit residual below the C₂z² residual.
* `QUADRATIC_LAW_EXTREMELY_FRAGILE` if the finite-record crossing or the z² consistency (M2A.1/M2A.2) already fails at |κ| = 0.05.
* `INCONCLUSIVE` otherwise.

## Stage 3 — parameter width (priority 4)

**3A numerical (continuum tangent, independent Gauss–Legendre implementation, cross-checked at b=10 against the frozen closed form to < 1e-10).** b ∈ {9.5, 9.6, …, 10.5}. Continuation of the A/B crossing from b = 10; per b: λ_∞(b), v_A, v_B, R_vv,A, R_vv,B, β_A, β_B, slope, coalescent gap, nearest regular competitor gap (scan v ∈ [−100,100], 20001 points + refinement), tail margin (analytic |v| ≥ 100 bound as in Stage 22). Pass per b: curvatures > 0, |β| ≥ 1e-3, slope ≥ 1e-3, competitor gap > 0, coalescent gap > 0, tail margin > 0, |v_A|,|v_B| ≥ 0.1. **Diagnostic interval** = largest r ∈ {0.1,…,0.5} with all grid b in [10−r, 10+r] passing (evidence `GLOBAL_NUMERICAL` for the tangent problem; never "certified").
Secondary: N=21, z=2 finite-record CA (Stage-1 protocol, rectangular) at each b.
**3B interval certificate attempt (after 3A).** Candidates in the given order [9.75,10.25], [9.90,10.10], [9.95,10.05], [9.98,10.02], [9.99,10.01]. Implementation: b-tiles of base width 0.005 processed centre-out; a candidate is certified iff every tile inside it passes (equivalent to testing candidates in descending order and stopping at the widest success). Per tile, uniformly in b: parametric Krawczyk inclusion for (v_A, v_B, λ) with contraction < 1; R_vv > 0; β_A > 0 > β_B; ∂_λ(R_A − R_B) > 0; Gram > 0; full-line exclusion on [−100,100] by cost/gradient/monotone/root-box predicates with zero unresolved cells; coalescent gap > 0; analytic tail gap > 0; exact-decimal coverage of v-cells and b-tiles. Outward interval arithmetic (mpmath.iv primary, same architecture as Stage 14), independent Arb (python-flint) replay of all final predicates for any success. Limits: tile b-bisection depth ≤ 4, v-cell width ≥ 1e-9, **total 3B wall-clock cap 150 min**; unfinished tiles count as not certified. If nothing succeeds: `NO_EXPLICIT_B_INTERVAL_CERTIFIED`.

## Stage 4 — strong-pair amplitude imbalance (priority 6)

A₋ = 1 − η (tone at −z), A₊ = 1 + η (tone at +z); η ∈ {0, ±0.02, ±0.05, ±0.10, ±0.20}; baseline otherwise; rectangular; CA on ε ∈ [0, 1.0] + boundary rule.
Classification: `PERSISTS_OVER_TESTED_IMBALANCE` if a qualifying crossing is present at all 9 values; `NARROW_IMBALANCE_ROBUSTNESS` if present for all |η| ≤ 0.05 but absent at some |η| ∈ {0.10, 0.20}; `DOES_NOT_PERSIST` if absent at some |η| ≤ 0.02; `INCONCLUSIVE` otherwise (any INCONCLUSIVE cell that decides the class).

## Stage 5 — matching convention (priority 7)

N ∈ {11, 21, 31, 41, 61}. Comparison A: z = 2, b = 10 fixed (ω's scale as 1/N). Comparison B: Δ = 2/21, ω₃ = 10/21 rad/sample fixed (z = 2N/21, b = 10N/21). CA (rectangular, ε ∈ [0,1.0], boundary rule). One table, one figure, one paragraph; `matching_language_audit.md` lists manuscript sentences that could be read as fixed-physical-frequency statements (no edits).

## Stage 6 — noisy estimator / local vs global (priority 5)

N=21, z=2, b=10, φ=π, ε = ε_c(1+η), ε_c = 0.2481906301722774; η ∈ {−0.20, −0.15, −0.10, −0.075, −0.05, −0.025, 0, 0.025, 0.05, 0.075, 0.10, 0.15, 0.20}; SNR ∈ {20, 30, 40} dB; noise exactly as FREEZE §7; seeds `SeedSequence([20260928, 6, snr, η-index, block])`, blocks of 50 trials.
**Estimators.** Production = *blind* GS(P) (grid + satellites + random starts; no oracle seeds). Reference = GS(R) + oracle seeds (noiseless A/B roots). Legacy = archived two-oracle-seed refinement (diagnostic only).
**6A.** Frozen subset: η ∈ {−0.10, −0.025, 0, 0.025, 0.10} × 3 SNR × first 40 records = 600 records. Disagreement = |J_prod − J_ref| > 1e-8·max(1,J_ref) or different branch label. Tolerance: ≤ 1% (≤ 6/600). If exceeded: double M_m, M_h, n_rand and K_grid of production once, rerun 6A, record in `AMENDMENTS.md`. Legacy disagreement rate is reported as a diagnostic of the archived experiment.
**Labels.** A if unordered-pair distance to the noiseless A root at that ε < 1.0; B likewise; else OUT.
**6B.** 1000 trials per (η, SNR) (39,000 records). Reduction rule (declared now): after the first two blocks of every setting, if projected wall-clock for the remainder > 150 min, reduce to 500 trials for all settings.
Metrics: P_A, P_B, P_OUT (Wilson 95%); empirical frequency distributions; global MSE to the noiseless global winner (A for η<0, B for η>0; at η=0 both references reported); MSE to the true strong pair (−2, 2); branch-conditioned mean, MSE (to own noiseless root) and 2×2 covariance.
**6C.** Local theory = frequency block of the misspecified local sandwich covariance (σ²/2)·H⁻¹Re(JᴴJ)H⁻¹ (6 real parameters, H includes the residual-curvature term) at each noiseless branch root — called **"local Hessian (sandwich) covariance approximation"**, not a CRB/MCRB.
Pre-registered hypotheses: **H1** for |η| ≥ 0.10 and SNR ∈ {30, 40}, for every branch with ≥ 100 selections, 0.8 ≤ trace(conditional MC covariance)/trace(local prediction) ≤ 1.25; **H2** at η = 0 for each SNR, global MSE (to A) ≥ 10 × local trace of A; **H3** at η = 0, min(P_A, P_B) ≥ 0.2 at each SNR. "Meaningful TSP consequence" = H1 ∧ H2 (H3 supporting). Failures are reported.

## Stage 7 — numerical invariance (priority 8)

Points: baseline crossing; weakest-third-margin qualifying crossing from Stages 1–4; Hann crossing; κ = 0.20 at z = 0.5 and fixed δ/π = 0.20 at z = 0.2; 50 noisy records at η = 0, 30 dB (Stage-6 seeds).
Variations: grid (M_m,M_h) ∈ {(512,129),(1024,257),(2048,513)}; n_rand ∈ {8,32,128}; K_grid ∈ {16,64,200}; tol_g ∈ {1e-7,1e-9,1e-11}; precision float64 vs mpmath-50 polishing of the crossing. Question: does any classification or noisy branch label change? (ε_c changes > 1e-6 relative are also flagged.)
Certificate path: directed-rounding containment tests (exact rational checks), transcendental enclosure tests vs 120-digit mpmath for mpmath.iv and Arb, static import audit of replay independence, and bitwise/parsed reproducibility of one certificate run repeated twice.

## Stage 8 — optional real sinusoids (only with ample time left)

Real record x_t = cos((ω₀−Δ)t) + cos((ω₀+Δ)t) − ε cos((ω₀+ω₃)t), Δ = 2/21, ω₃ = 10/21, N = 21, ω₀ ∈ {π/2, 6π/21}; fit two real sinusoids (four real amplitudes, two frequencies in (0,π)); grid + refinement; ε ∈ [0, 1.0]. Tag `OPTIONAL_EXPLORATORY`.

## Final decision logic

Exactly the user's hierarchy (STRONG / NEEDS ONE FOCUSED VALIDATION PASS / MAJOR RESCOPE RISK / NOT READY). Operationalisation fixed now:
* "exchange persists under at least one/two reasonable nonrectangular windows" = ≥ 2 of {Hann, Hamming, DPSS} `PERSISTS`.
* "parameter robustness is nontrivial" = Stage 3A diagnostic r ≥ 0.2 **and** Stage 4 `PERSISTS_OVER_TESTED_IMBALANCE` or `NARROW_IMBALANCE_ROBUSTNESS` with |η| ≤ 0.10 passing.
* "clean distinction between general exchange and symmetry-specific z² law" = Stage 2 class `QUADRATIC_LAW_LOCALLY_STABLE_BUT_ALIGNMENT_SPECIFIC` with qualifying crossings also present in the fixed-δ family.
* "meaningful TSP consequence" = H1 ∧ H2.
* "major numerical confound" = any Stage-7 classification change, or Stage 6A disagreement > 1% after the one allowed production upgrade.
