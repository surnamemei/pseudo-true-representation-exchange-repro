# Numerical validation and figure data

Nothing in this directory is an interval certificate.

- **`adversarial_overnight/`.** The preregistered adversarial validation:
  - `PREREGISTRATION.md` and `prereg_hash.json`, recorded before any new experiment;
  - `FREEZE.md` and `freeze_manifest.json`, the SHA-256 of 766 project files;
  - `AMENDMENTS.md` (A1–A4, with timing);
  - stage scripts (`s1_window.py` … `s8_real.py`, `advcore.py`, `run_all.py`), and the scans, analyses and reports of each stage;
  - `stage0_baseline/replay_comparison.json`: the Windows replay of the certificate chains.

  Three items are omitted: the Stage-0 replay sandbox, which is a byte copy of the project; the raw per-leaf logs of the incomplete b-interval attempt (per-block summaries are kept in `stage3_bwidth/bcert/blocks/*.json`); and the end-of-run summaries written for the author (`MORNING_REPORT.md`, `CLAIM_AUDIT.md`, `HOSTILE_REVIEW.md`, `FINAL_SUMMARY.txt`). The preregistered results, including the overall status defined in the preregistration, are in `STATUS.json`, `RESULT_MATRIX.csv`, `stage*/DONE.json` and the stage reports.
- **`analysis/`.** Numerical results used by the manuscript's figures and text:
  - amplitude continuation;
  - large-N coefficient table;
  - z-continuation comparison;
  - continuum b map;
  - order-selection trials;
  - archived phase, resolution and noise maps;
  - the attempted connected p-bridge ledger. The bridge is **not** certified.
- **`figures/`.** `make_revision_figures.py` redraws every manuscript figure from these archived results; run it with `python replay/replay_and_compare.py figures`. `inputs/` holds the archived figure file that the script copies.

Evidence levels: robustness results are global numerical searches under the preregistered protocol. The real-valued study is exploratory. The fixed-phase exponent fits are post-hoc sensitivity analyses, and the preregistered Stage-2 classification is inconclusive. See `../docs/EVIDENCE_LEVELS.md`.
