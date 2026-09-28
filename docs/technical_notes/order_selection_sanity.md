# Stage 21: model-order sanity check at the certified baseline

**Result.** BIC selected order three in all 1,080 tested records. This does not support accidental BIC underfitting at the studied 20--40-dB baseline. The deterministic 3-to-2 switch remains relevant to a deliberately imposed two-tone representation or a downstream method with fixed order. This experiment is an interpretation check, not a paper contribution or a theorem.

## Fixed design and criterion

The signal is the phase-aligned N21, z=2, b=10, weak-phase-pi family. We tested `eta=(epsilon-epsilon_c)/epsilon_c` equal to -0.10, 0, and 0.10, with `epsilon_c=0.2481906301722774`, at SNR 20, 30, and 40 dB. Each of the nine cells contains 120 independent records, with RNG seed 20260928. As in the paper's noise experiment, `sigma2=(||x||_2^2/N)*10^(-SNR/10)` is the **complex** sample variance, and each real/imaginary component has variance `sigma2/2`.

Both candidate orders use continuous-frequency variable-projection nonlinear least squares: at every trial frequency vector, unrestricted complex amplitudes are eliminated by least squares. The two-tone search starts from the deterministic A, B, and strong-pair seeds; every twentieth record is additionally audited by a 128-by-128 full-torus grid and 16 local refinements. None of the 54 audited records improved the best seeded two-tone fit by more than `1e-7`. The three-tone search uses the generating frequencies, augmentations of A/B, and a residual-periodogram seed; its five local optimizations retain the smallest residual. All selected three-tone gradient norms are below `9e-6`. These are careful numerical fits, not certified global optima.

With 2N independent real observations and 3K real mean parameters at order K (two per complex amplitude and one frequency), the common unknown noise-variance parameter cancels in

`DeltaBIC = BIC_3-BIC_2 = 2N log(RSS_3/RSS_2) + 3 log(2N)`.

Negative `DeltaBIC` selects order three. Using `N=21`, the penalty difference is `3 log(42)`. We did not adjust it after seeing the data.
Counting N complex observations instead would use the smaller penalty `3 log(21)` and therefore also select order three in every tested record.

## Results

| SNR dB | eta | order 2 / 120 | order 3 / 120 | DeltaBIC median | DeltaBIC 5--95% sample quantiles |
|---:|---:|---:|---:|---:|---:|
| 20 | -0.10 | 0 | 120 | -24.68 | [-40.55, -10.96] |
| 20 | 0 | 0 | 120 | -27.56 | [-44.08, -14.25] |
| 20 | +0.10 | 0 | 120 | -28.55 | [-44.94, -13.25] |
| 30 | -0.10 | 0 | 120 | -98.08 | [-114.72, -84.62] |
| 30 | 0 | 0 | 120 | -103.93 | [-125.65, -89.26] |
| 30 | +0.10 | 0 | 120 | -109.13 | [-121.60, -90.07] |
| 40 | -0.10 | 0 | 120 | -194.32 | [-207.15, -178.80] |
| 40 | 0 | 0 | 120 | -203.22 | [-217.97, -186.57] |
| 40 | +0.10 | 0 | 120 | -203.67 | [-222.18, -187.04] |

The largest observed criterion difference is -0.144, so one 20-dB record was close to a BIC tie; the remaining 1,079 still favored three tones. Zero order-two selections among 360 trials at 20 dB gives an approximate one-sided 95% binomial upper bound of 0.83% for the aggregate tested mixture of three eta values. It does not bound other SNRs, epsilon values, or other selectors.

The [trial table](order_selection_sanity.csv) contains every residual, selection, criterion difference, fitted frequency, and audit flag. [Grouped quantiles](order_selection_sanity_summary.csv) and the [empirical CDF plot](order_selection_sanity.pdf) are derived from it. Regenerate all three with `python stage21_order_selection.py --trials 120`.
