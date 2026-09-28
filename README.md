# Global Pseudo-True Representation Exchange in Under-Modelled Spectral Fitting: certificates and reproducibility code

[![Zenodo DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23013068.svg)](https://doi.org/10.5281/zenodo.23013068)

When a spectral model with a fixed number of sinusoids is fitted by least squares to a record that contains more components, it estimates a *pseudo-true* lower-order representation of the finite record. For a specified family of finite three-tone records fitted with two tones, the globally optimal two-tone representation switches discontinuously between two distinct, separated, locally regular fits as an omitted weak component grows. For a centred, phase-aligned close pair, the critical omitted amplitude follows ε_c = λ_N z² + c_N z⁴ + O(z⁶) in the normalized spacing z. The global statements rest on analytical arguments whose computational hypotheses are verified by directed-rounding interval arithmetic.

This repository is a code-and-data companion, not a copy of the manuscript. It contains:
- **Code.** The interval-arithmetic certificate scripts (mpmath.iv), and an independent Arb/python-flint replay implementation.
- **Records.** The machine-readable certificate records, the preregistered numerical validation behind the manuscript's robustness figures, and frozen SHA-256 manifests.
- **Tooling.** A one-command replay that reruns the frozen scripts and compares their outputs with the archived records.
- **First-order transition theory.** The sealed analysis behind the manuscript's noisy-estimation section: a first-order law predicting noisy branch selection and the global error near the exchange, with its one-shot validation (`transition_theory/`).

## What is certified

The certified results use directed interval arithmetic together with analytical arguments. The configuration is a centred, phase-aligned strong pair with the weak tone at normalized location b = 10 and phase π.

- **Small spacing.** A checkable global criterion (C1–C8) for the exchange at small spacing, verified for N = 11, 15, 21, 31 and 41. The small-spacing radius is existential.
- **Continuum.** The full-line continuum tangent crossing.
- **Large N.** Transfer to **every odd N ≥ 10,001 at fixed normalized geometry**: z = NΔ and b = Nω₃ are held fixed. This is not a statement about fixed physical frequencies.
- **Finite spacing.** Full-domain global exchanges at spacing z = 2 for N = 21 and N = 31.
- **Finite width.** For N = 21 at z = 2, a certified amplitude interval [0.213190630172278, 0.283190630172277]. On it, A and B are the only globally competitive fits, and their cost difference has a unique transverse zero.
- **Lower endpoint.** The N = 21 lower-endpoint constants: the coalescent class is never globally optimal for positive omitted amplitude.

**The independent Arb replay does not cover every certified result.**
- **Covered.**
  - Every obligation of the two z = 2 crossing certificates (68,224 checks).
  - The continuum partition, uniformly for 0 ≤ N⁻² ≤ 35,377⁻².
- **Not covered.** The N ≥ 10,001 floor, the finite-width interval, the five tangent instances and the lower-endpoint constants have been checked only by the primary implementation, through original-code replays.

## What is numerical only

- **Global numerical robustness** (full-torus searches under a preregistered protocol):
  - weighting by Hann, Hamming, DPSS and Hann² windows;
  - weak-tone location b ∈ [9.5, 10.5]. This is numerical, **not a certified b interval**: the certified statement is an existential neighbourhood of b = 10, and an explicit b-interval certificate attempt did not complete;
  - strong-amplitude imbalance up to 20%;
  - strong-pair phase and record-origin offsets;
  - matching conventions across N;
  - a Monte Carlo study of local covariance against global error;
  - solver and precision invariance.
- **Numerical, outside the protocol:**
  - the continuum b map at integer 6 ≤ b ≤ 13;
  - a BIC order check;
  - the numerical continuation from small spacing to z = 2. No certified connected branch is claimed.
- **Exploratory:** real-valued sinusoids.
- **Post-hoc sensitivity:** fixed-phase exponent fits and the 20 dB nearest-root labelling.

**Validated first-order approximation** (not a theorem, certificate or bound):
- **What it predicts.** The law P(A) ≈ Φ(−ΔJ₀/(√2σ‖r_A − r_B‖₂)) ≈ Φ(−Kη) and a two-branch mixture predict, from noiseless branch quantities and the prescribed noise level alone:
  - the branch-selection probability near the N = 21, z = 2 exchange;
  - its 10–90% width, which scales as σ;
  - the global mean squared error.
- **Validation at 30 and 40 dB.**
  - Predicted widths are 0.1285 and 0.0406, against observed 0.1295 and 0.0419.
  - The mean |ΔP(A)| is 0.005 and 0.002.
  - The global MSE is within a factor 0.74–1.35 for |η| ≤ 0.1.
- **Provenance.** The theory was sealed before the per-setting Monte Carlo outcomes were read.
  - The empirical widths were already known, so the width comparison is not blind.
  - No new Monte Carlo was run.
- **Stress case.** At 20 dB, displaced A-family fits fall outside the approximation.

**The critical law is not universal.** ε_c = λ_N z² + c_N z⁴ + O(z⁶) is proved for the centred phase-aligned family, where symmetry removes the first-order term. A fixed nonzero strong-pair phase gives an O(z) regime.

## Quick start

```
git clone https://github.com/surnamemei/pseudo-true-representation-exchange.git
cd pseudo-true-representation-exchange
python3 -m venv .venv && . .venv/bin/activate            # Windows: .venv\Scripts\activate
pip install -r requirements.txt                           # or: conda env create -f environment.yml
sha256sum -c manifests/SHA256SUMS                         # optional integrity check
python replay/replay_and_compare.py all R1 R3 R4 R7 T1    # short check, about 3 minutes
python replay/replay_and_compare.py all                   # all eight chains in parallel, just under an hour
python replay/replay_and_compare.py figures               # redraw the manuscript figures from archived results
```

The replay tool works in a scratch copy, so archived records are never overwritten:
1. It assembles the flat working directory in which the scripts were frozen (`_replay_work/`, from `replay/workspace_map.csv`).
2. It deletes every declared output.
3. It runs the frozen scripts unmodified, one Python process per command.
4. It compares each regenerated output with the archived record. Verdicts, labels and counts must match exactly; numeric enclosure values may differ only in far digits (relative difference ≤ 1e-10).

The per-chain commands and what each establishes are in [`docs/REPLAY_COMMANDS.md`](docs/REPLAY_COMMANDS.md). The finite-width chain starts with `python stage17_epsilon_refine.py 0.035 14287`: all 14,287 outer parent cells must be refined.

## Expected runtimes

Times are wall-clock seconds per chain, measured on an AMD Ryzen 9 9950X3D. Each command is a single-threaded process.

| Chain | Establishes | Windows 11, Python 3.13 | Linux (WSL2), Python 3.12 |
|---|---|---:|---:|
| R1 | Independent Arb replay of both z = 2 crossing certificates | 119 | 89 |
| R2 | Five tangent instances (criterion C1–C7) | 2,089 | 1,433 |
| R3 | Continuum crossing certificate | 155 | 116 |
| R4 | Independent Arb replay of the continuum partition, N ≥ 35,377 | 45 | 30 |
| R5 | Every odd N ≥ 10,001, fixed normalized geometry | 697 | 450 |
| R6 | N = 21 finite-width amplitude interval | 4,391 | 3,202 |
| R7 | N = 21 lower-endpoint constants | 10 | 11 |
| T1 | Deterministic predictions of the first-order transition theory (not a certificate) | — | 9 |
| | **Total, single process** | **7,506 (2.1 h)** | **5,331 (1.5 h)** |

The Linux replay of 2026-09-28 is recorded in [`replay/records/linux_wsl_2026-09-28/`](replay/records/linux_wsl_2026-09-28/): exit codes, runtimes, logs, output-by-output comparison, environment and file-access trace. All seven chains exited 0, and every declared output matched the archived record. Each output was byte-identical, or identical up to CRLF/LF line endings, JSON key order or elapsed-time fields. The one exception is floating-point-preconditioned enclosure values in the N ≥ 10,001 floor and the finite-width root tiling, which differ only in far digits (relative difference ≤ 3 × 10⁻¹⁴). No verdict, label or count changed.

## Hardware and software notes

- **Python and packages.** Python ≥ 3.11 with the pinned packages in `requirements.txt`: NumPy, SciPy, mpmath 1.3.0, SymPy, Matplotlib, python-flint 0.9.0 (FLINT 3.6.0) and tqdm. The frozen runs used CPython 3.13.9 on Windows 11 (NumPy/MKL). The seven chains and this replay tool were rerun with CPython 3.12.3 on Ubuntu 24.04 under WSL2 (NumPy/OpenBLAS). The tool is written to be portable, but it has been tested on Linux only.
- **Arithmetic.**
  - The primary certificates use `mpmath.iv` with directed rounding at 45–110 decimal digits.
  - The independent replay uses Arb balls at 320-bit midpoint precision.
  - Floating point proposes candidates, partitions and preconditioners only.
- **Platform differences.** Preconditioners can differ between BLAS builds, which moves some enclosure bounds in far digits (relative 1e-14 or smaller in the recorded replay). No verdict changes. Files written in text mode on Windows have CRLF line endings, so a byte comparison needs line-ending normalization. The repository's `.gitattributes` turns off Git's line-ending conversion, so every clone, including one on Windows, keeps the archived bytes and passes the hash check.
- **The `stage12_deps/` path.** The frozen Arb scripts prepend a directory `stage12_deps/` to the import path. It is absent here, and python-flint is taken from the environment. Do not put a Windows build of python-flint there on Linux.
- **Resources.** A full replay needs about 1 GB of free disk and at least 8 GB of RAM. The largest command, R6's inner tiling, peaked at about 3.8 GB; the replay tool records peak memory per command on Linux and macOS. The chains are independent, except that R4 and R5 read outputs of R3; the tool orders them.

## Version correspondence

This repository corresponds to the manuscript *Global Pseudo-True Representation Exchange in Under-Modelled Spectral Fitting* by Jinghang Mei, submitted to *IEEE Transactions on Signal Processing* (2026).
- **Frozen scripts and records.** They are those of the reviewer archive `TSP_REVIEWER_ARCHIVE_2026-09-28_v2`, whose ZIP has SHA-256 `0c00688c760e135e4c508404042389b175958e9bb29381d90ef145bc80d80567`.
  - That archive accompanies the Stage 25 manuscript, which integrates the approved first-order transition theory.
  - It supersedes the first archive, `TSP_REVIEWER_ARCHIVE_2026-09-28`.
  - No certificate script, certificate record or raw validation output differs between the two.
- **Hash checks.** Every script and record here that is also listed in the project's 766-file freeze manifest (`numerical_validation/adversarial_overnight/freeze_manifest.json`) keeps its frozen SHA-256; `manifests/files.csv` records the check. The only exception is `docs/REPRODUCIBILITY_README.md`, documentation edited after the validation freeze. No certificate code or data differ.
- **What is not here.** The manuscript and its figures are not part of this repository.

## Evidence levels at a glance

| Level | Results |
|---|---|
| **Certified** | Small-spacing criterion (N = 11, 15, 21, 31, 41); continuum crossing; every odd N ≥ 10,001 at fixed normalized geometry; z = 2 exchanges (N = 21, 31); N = 21 finite-width amplitude interval; N = 21 lower-endpoint constants |
| **Independent replay** (Arb; separate code, same archived partitions) | Both z = 2 crossing certificates; the continuum partition for 0 ≤ N⁻² ≤ 35,377⁻² |
| **Original-code replay** (the frozen scripts rerun unmodified, on Windows and on Linux, 2026-09-28) | All seven certificate chains, and T1 for the transition theory. This is the only replay of the N ≥ 10,001 floor, the finite-width interval, the tangent instances and the lower-endpoint constants. |
| **Global numerical robustness** (preregistered) | Windows; b ∈ [9.5, 10.5]; imbalance ≤ 20%; phase and record-origin offsets; matching conventions; noise Monte Carlo; solver invariance |
| **Validated first-order approximation** (sealed before the per-setting Monte Carlo outcomes were read; compared once) | Noisy branch-selection probabilities, their σ-proportional width, and the global error near the N = 21, z = 2 exchange (30 and 40 dB; 20 dB is a stress case) |
| **Exploratory** | Real-valued sinusoids |
| **Post-hoc sensitivity** | Fixed-phase exponent fits; 20 dB nearest-root labelling |

The full table, with scopes, qualifications and record files, is in [`docs/EVIDENCE_LEVELS.md`](docs/EVIDENCE_LEVELS.md). It also lists what is explicitly *not* claimed:
- a certified b interval;
- a certified small-to-finite-spacing branch;
- fixed-physical-frequency large-N behaviour;
- a universal z² law;
- a CRB or MCRB failure;
- an order-selection failure.

## Repository layout

| Path | Contents |
|---|---|
| `src/primary_interval/` | Primary directed-interval (mpmath.iv) evaluators and the drivers of chains R2, R3, R5, R6 and R7 |
| `src/independent_arb/` | Independent python-flint/Arb implementation (R1, R4) and the exact-decimal cover audit. A static audit shows it imports no primary evaluator. |
| `src/analysis/` | Numerical, non-certificate analysis scripts behind archived data and figures: continuation, asymptotic coefficients, b map, order check, noise-width uncertainty |
| `src/provenance/` | Earlier-stage scripts, including the generators of the archived candidate partitions. Kept for provenance; the replay does not rerun them, and their historical intermediate data are not included. |
| `certificates/` | Machine-readable certificate records: inputs, candidate partitions and outputs of R1–R7; earlier certificate records; `stage12_generalN/`; `stage12_replay/` |
| `numerical_validation/adversarial_overnight/` | Preregistered adversarial validation: protocol, preregistration hash, freeze manifest, amendments, stage scripts, scans, analyses, reports, and the Stage-0 replay record |
| `numerical_validation/analysis/` | Numerical results used by the manuscript's figures and text |
| `numerical_validation/figures/` | The figure script, and the archived inputs it copies |
| `transition_theory/` | The sealed first-order transition theory: plan, spec, seals, derivation, code, validation report, and its predictions and comparison (`results/`) |
| `replay/` | `replay_and_compare.py`, `workspace_map.csv`, frozen hashes of outputs too large to include, and replay records |
| `manifests/` | `SHA256SUMS` and `files.csv` for this repository, the Stage 24 and Stage 25 revision manifests, integrity-check records, and historical artifact manifests |
| `docs/` | Evidence levels, replay commands, the full project reproducibility guide (`REPRODUCIBILITY_README.md`, with historical sections), and technical notes: derivations, lemmas, certificate summaries, replay reports and formula audits |
| `examples/` | Short command sequences for common checks |

Scripts keep their research-stage names and bytes, so their SHA-256 matches the frozen manifests. They expect the flat working directory in which they were frozen; run them through the replay tool, which rebuilds that layout.

## Not included

| Item | Reason |
|---|---|
| The 539 MB finite-width inner-tiling intermediate (R6, step 3) | Regenerated by the replay. Its frozen SHA-256, both as frozen and with LF line endings, is in `replay/frozen_output_hashes.json`. |
| Records of the superseded smaller finite-width radii | No claim uses them |
| Raw per-leaf logs of the incomplete b-interval attempt (97 MB) | Per-block summaries are included |
| Replay sandboxes, byte-code caches, the Windows python-flint build | Temporary or platform-specific |
| A duplicate copy of the Arb per-obligation table | `certificates/arb_full_replay_comparison.csv` is the only copy |
| Manuscript sources, PDFs and figures; cover letter; submission metadata; author-internal reviews, audits and end-of-run validation summaries; historical backups; the obsolete first submission snapshot | Not code or evidence, or not public. The preregistered validation results are all included. |

## Licence and citation

- **Code:** MIT (`LICENSE`).
- **Certificate records, validation data and documentation:** CC BY 4.0 (`LICENSE-DATA.md`).
- **Citation:** see `CITATION.cff`.
- **Contact:** Jinghang Mei, ORCID [0009-0007-2901-3285](https://orcid.org/0009-0007-2901-3285).
