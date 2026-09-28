# Source code

Every script keeps its frozen research-stage name and bytes. The scripts expect the flat working directory they were frozen in, so run them through `../replay/replay_and_compare.py`, which rebuilds that layout in `_replay_work/`.

| Directory | Contents |
|---|---|
| `primary_interval/` | Directed-rounding interval evaluators (`mpmath.iv`) and the drivers of chains R2, R3, R5, R6 and R7, with the earlier-stage modules they import. Also here: the superseded N ≥ 35,377 floor certificate (`stage15_uniform_interval.py`, `stage15_finalize_certificate.py`). Floating point proposes candidates, partitions and preconditioners; every decision is made in interval arithmetic. |
| `independent_arb/` | The second implementation, in python-flint/Arb: `stage12_arb_backend.py`, `stage12_full_arb_replay.py` (R1), `stage16_independent_arb.py` (R4), and `stage23_b_arb_replay.py` (Arb check of the early Stage-23 b-interval boxes, all of which failed; not a claim). A static import audit shows that none imports a primary evaluator. `stage12_cover_audit.py` checks exact-decimal coverage of the archived partitions, and `stage12_reports.py` writes the replay summaries. |
| `analysis/` | Numerical, non-certificate analysis scripts: amplitude continuation, asymptotic coefficients, the continuum b map, the order-selection check, noise-width uncertainty, the z-continuation comparison, and a closed-form coefficient check (`kcheck.py`). Their outputs are numerical evidence only. |
| `provenance/` | Earlier-stage scripts, including the floating-point search code that proposed the archived candidate partitions. They are kept so that the history can be read. The replay does not rerun them, and their historical intermediate data are not included. |
