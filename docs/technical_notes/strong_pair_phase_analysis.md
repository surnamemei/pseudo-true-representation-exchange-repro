# Stage 20: strong-pair relative phase audit

**Verdict: PHASE-GENERALIZED.** The phase-aligned \(z^2\) law is not an isolated asymptotic point: a global A/B small-spacing exchange persists for an existential open neighborhood of scaled phase offsets \(\delta=\kappa z\). A genuinely different fixed-nonzero-\(\delta\) competition scale is \(\epsilon=\Theta(z)\), and its N21 tangent landscape has a strong numerical A/B crossing. That fixed-phase *global crossing coefficient* is not interval-certified, so the stronger label `TWO-REGIME-THEORY` is premature.

## What changed mathematically

For fixed \(\delta\ne0\), the pair has first-order term \(-2z s\sin\delta\). A single central complex-amplitude/real-frequency exponential produces an imaginary, not real, linear frequency tangent. The correct leading problem must include a two-central confluent chart as well as a central-plus-satellite chart. The exact regular scalar loss is \(\lambda^2H_N(v)+4\sin^2\delta P_N(v)\) at \(\epsilon=\lambda z\). At N21, the numerical normalized crossing gives \(\epsilon_c/z\approx0.6274473964|\sin\delta|\) for the five tested phases, with positive sampled coalescent and third-regular gaps. See `fixed_delta_scaling.md` and `fixed_delta_phase_map.csv`.

For \(\delta=\kappa z\), the target is \(-(s+\kappa)^2-\lambda e^{ibs}\) at order \(z^2\). Exact parity reduction gives \(R_{N,\kappa}=R_{N,0}+4\kappa^2P_N\). Certified aligned-case strict margins and the implicit-function theorem prove persistence for some unspecified \(|\kappa|<\kappa_0(N)\). The continuum crossing is also locally stable. The crossing coefficient is even in \(\kappa\); its finite N21 and continuum quadratic terms are \(1.9283886651\kappa^2\) and \(1.9441711059\kappa^2\), respectively. See `kappa_scaled_phase_theorem.md` and `kappa_crossing_map.csv`.

The N21 numerical continuation encounters a vanishing coalescent margin near \(|\kappa|=0.3882278594\), so no connected noncoalescent theorem is asserted across that value. Values at \(\pm0.50\) are separately sampled regular crossings, not an extension of the certified neighborhood.

## Physical timing interpretation

The strong pair is \(e^{i\delta}e^{izs}+e^{-i\delta}e^{-izs}=2\cos(zs+\delta)\), with relative phase \(2\delta\). Shifting the observation-window center by \(\tau\) samples from a beat-envelope center gives \(\delta=\Delta\tau=(\tau/N)z\). Thus a fixed fractional offset \(\tau/N=\kappa\) belongs to the quadratic scaled-phase regime. Holding \(\delta\) fixed as \(z\to0\) instead moves the envelope center a growing number of record lengths away. Both parameter families are mathematically legitimate; neither is universally the physical one.

## Related-work boundary

The deterministic representation crossing is distinct from the **noise-driven** subspace-swap threshold analyzed by [Pakrooh, Scharf, and Pezeshki, IEEE TSP 2016](https://doi.org/10.1109/TSP.2016.2521617) and the single-tone ML threshold/outlier behavior of [James, Anderson, and Williamson, IEEE TSP 1995](https://doi.org/10.1109/78.376834). The pseudo-true-target framing is established by [Fortunati et al., IEEE SPM 2017](https://doi.org/10.1109/MSP.2017.2738017). These are context and terminology, not exact collisions or evidence for the present phase theorem.

## Proof and numerical boundaries

- Proven: first-order versus second-order asymptotic charts; exact scalar tangent identities; even \(\kappa\) dependence; existential open-\(\kappa\) global persistence from the existing aligned certificates; analytic formula for crossing curvature.
- Numerical: five fixed-phase global tangent scans and their critical coefficient; all displayed nonzero \(\kappa\) rows; the approximate coalescent contact; displayed high-precision derivative values.
- Unproved: an explicit certified phase radius; interval-global fixed-phase crossing; finite-\(z\) connected phase map; continuation across the coalescent contact.

The independent projector audit in `theta_factor_audit.md` found the correct \(\sqrt2 Lr_m\) implementation. The revised manuscript restricts its original finite-\(z\) certificates to record-center phase alignment. The active Stage-19 ZIP remains the prior frozen submission state and is not silently overwritten by this Stage-20 analysis.
