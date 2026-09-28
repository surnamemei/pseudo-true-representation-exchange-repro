# Linux (WSL2) replay of the seven certificate chains, 2026-09-28

After the project moved from Windows (`F:/Research`) to WSL2, the seven theorem-level certificate chains R1–R7 were rerun from the frozen scripts, unmodified, in a scratch copy of the project. Every output was then compared with the frozen record.

The commands and output lists are those of `validation/adversarial_overnight/tools/replay_sandbox.py`, the Stage-0 tool behind Supplementary Table S5. The finite-width chain starts with `stage17_epsilon_refine.py 0.035 14287`.

**Result: all seven chains exit 0, and every declared output passes the comparison. No verdict, label or integer count differs from the frozen record.**

| Chain | Linux time (s) | Outcome |
|---|---:|---|
| R1 finite-spacing Arb replay | 89.3 | 4 outputs byte-identical; 2 identical after CRLF→LF; the 2 summaries differ only in `elapsed_seconds` |
| R2 five tangent instances | 1432.7 | 10 CSVs byte-identical; 10 JSONs identical after CRLF→LF |
| R3 continuum certificate | 115.7 | 3 byte-identical; 4 identical after CRLF→LF; 1 same content with a different JSON key order |
| R4 independent Arb, large N | 29.9 | both outputs byte-identical |
| R5 N ≥ 10,001 floor | 449.7 | numeric enclosure values differ at relative ≤ 4.2 × 10⁻²²; 3 SHA-256 fields hash regenerated inputs; all flags and counts identical |
| R6 finite-width ε interval | 3202.1 | 6 identical after CRLF→LF, including the 539 MB inner-tiling intermediate; 1 byte-identical; root tiling differs in 2,791 Krawczyk contraction bounds at relative ≤ 2.9 × 10⁻¹⁴; certificate JSON differs only in 5 SHA-256 fields of regenerated inputs |
| R7 lower-transition constants | 11.2 | 2 byte-identical; 1 identical after CRLF→LF (plus two diagnostic figure PDFs, whose embedded metadata differ) |

## Interpretation

- **Line endings.** Windows text-mode writes produced CRLF line endings; Linux writes LF. The CSV writer uses CRLF on both platforms, so CSVs match byte for byte.
- **Numeric differences.** The only numeric differences occur where a floating-point quantity (a Krawczyk preconditioner or a candidate centre) enters an interval computation. NumPy's BLAS differs between the Windows runs (MKL) and this environment (OpenBLAS). Such inputs are allowed to vary: every decision is made in interval arithmetic, and all predicates, verdicts and counts are identical.
- **python-flint.** The Arb replays (R1, R4) used the PyPI wheel of python-flint 0.9.0 (FLINT 3.6.0). The project's `stage12_deps/` holds a Windows build that cannot load on Linux. It was not copied into the sandbox and was not modified. Run in place under WSL, `stage12_arb_backend.py` fails to import flint because of that directory.

## Provenance of the run

- **Setup.** Sandbox: all top-level project files of 50 MB or less, plus `stage12_generalN/` and `stage12_replay/`. Every declared output was deleted first (`sandbox_info.json`).
- **Execution.** Chains ran concurrently, one single-threaded process per command; `environment.json` gives the Python, package and BLAS details.
- **Harness defect in R5.** The thread that started R5 crashed before running anything. It read the shared state file while another thread was writing it; this was a race in the harness, not a script failure. The harness was fixed to write atomically, and R5 was then run by itself in the same sandbox after R3 had finished. `replay_state.json` shows the resulting exit codes and times.
- **Tracing.** A read-only Python audit hook (`sitecustomize.py`, loaded via `PYTHONPATH`) recorded every file each chain opened and every local module it imported (`trace_summary.json`). This trace defines the replay dependencies packaged in the reviewer archive and the public repository.

## Figures

`paper/figure_scripts/make_revision_figures.py` was also rerun (exit 0) in a copy of the project. `figure_trace.json` lists its inputs. `figure_regeneration.json` compares the output with the frozen figures:
- the two copied figures are byte-identical;
- the others differ only in font rendering and PDF metadata, with the same plotted data.

## Files

`replay_state.json`, `output_comparison.json`, `environment.json`, `sandbox_info.json`, `trace_summary.json`, `figure_trace.json`, `figure_regeneration.json`, and `logs/` (one log per command). Absolute scratch paths are replaced by `<sandbox>`, `<venv>` or `<scratch>`.
