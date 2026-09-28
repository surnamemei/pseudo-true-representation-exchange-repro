# Full second-backend replay summary

**PASS for both complete finite-spacing certificates.**

|N|Archived obligations|PASS|FAIL|Arb-refined inner parents|Validated child leaves|
|---:|---:|---:|---:|---:|---:|
|21|39792|39792|0|128|256|
|31|28432|28432|0|44|88|

Total: 65,337 outer cells, 2,801 inner cells and 86 local/crossing obligations = **68,224 checks**. All archived parents are validated. Arb uses one bisection for 172 inner parents; another 21 parents use a different exclusion predicate. No contradictory result was found. The full arithmetic check is accompanied by an exact-decimal union-cover audit for all four outer case partitions and eight inner neighborhoods.

The second backend is python-flint 0.9.0 / FLINT 3.6.0 / Arb at 320 bits, separate from the original mpmath.iv endpoint arithmetic. Its fresh energy, objectives, analytic derivative jets, projector displacement and Krawczyk evaluations do not read the old inequality values as proof inputs. Candidate centers, parameter brackets and cell coordinates are shared certificate data. This is an independent arithmetic replay, not an independently discovered theorem or formal proof-assistant verification.

See `arb_full_replay_N21.md`, `arb_full_replay_N31.md`, the two complete comparison CSVs, and the per-parent refinement leaves. Historical 44-check MPFR results remain archived but are no longer the extent of independent verification.
