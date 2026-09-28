# Replay

- `replay_and_compare.py`: builds the workspace, runs chains R1–R7 and compares their outputs with the archived records. It also redraws the figures. Run `python replay/replay_and_compare.py -h` for usage.
- `workspace_map.csv`: maps each workspace path, the flat layout in which the scripts were frozen, to the file in this repository.
- `frozen_output_hashes.json`: SHA-256 of chain outputs too large to include.
- `records/linux_wsl_2026-09-28/`: the Linux replay of 2026-09-28. It contains exit codes and runtimes (`replay_state.json`), per-command logs, the output-by-output comparison (`output_comparison.json`), the environment (`environment.json`), the file-access trace of every chain (`trace_summary.json`) and the figure regeneration check.
- The earlier Windows replay (Stage 0 of the validation) is in `../numerical_validation/adversarial_overnight/stage0_baseline/replay_comparison.json`.
