# Examples

Run everything from the repository root after `pip install -r requirements.txt`.

| Script | What it does | Time |
|---|---|---|
| `verify_integrity.sh` | Checks every file against `manifests/SHA256SUMS` | seconds |
| `replay_independent_arb.sh` | Arb replay of both z = 2 crossing certificates (R1), and of the continuum partition for large N (R3, then R4) | ~3 min |
| `replay_finite_width_interval.sh` | The N = 21 finite-width amplitude chain (R6), starting with `stage17_epsilon_refine.py 0.035 14287` | ~55 min |
| `replay_all.sh` | All seven chains in parallel, then the comparison with the archived records | ~1 h |
| `redraw_figures.sh` | Redraws the manuscript figures from the archived numerical results | ~1 min |
| `inspect_certificate.py` | Prints the key fields of the main certificate records | seconds |
