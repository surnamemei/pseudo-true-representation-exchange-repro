# Complete independent Arb replay: N=31

**PASS — all 28,432 archived cell/local/crossing obligations validated.**

- Backend: python-flint 0.9.0, FLINT 3.6.0, Arb, 320-bit midpoint precision.
- Outer categories: excluded=27481, near_A=69, near_B=34.
- Inner terminal cells: 825.
- Full root inclusion, contraction, positive Hessian principal minors, positive Gram, crossing endpoint signs, and positive cost-gap slope recomputed.
- Minimum outer exclusion margins by case: {'crossing': 0.00010809124995015376}.
- Original inner predicates needing local Arb refinement: 44; one bisection each, 88 validated children, 132 refinement tree nodes including parents.
- Total different inner proof paths: 49 (including refinements); remaining differences use another valid exclusion predicate.
- Contradictory certificate results: **0**. Unresolved parent cells: **0**.
- Runtime: 45.33 s for this replay; timing excludes coverage audit.

## Meaning of agreement

This is complete validation of the archived mathematical obligations, with local refinement where Arb ball propagation is wider. It is not bitwise agreement between endpoint and ball arithmetic. An inconclusive unrefined Arb bound is an enclosure limitation; both child cells prove the same parent's exclusion. The archived files are unchanged. The replay does not use stored lower bounds as premises: it freshly evaluates the signal energy, incumbent, projector change, gradients, Hessians and preconditioners.

All outer exclusions, accepted-cell containments and inner terminal cells are checked, including the confluent edge of the outer domain. N21 includes below/crossing/above amplitude records, with the entire narrow crossing bracket used in the crossing case; N31 includes its full crossing bracket. Exact-decimal rectangle-union audits verify no gap in the outer or inner covers. Radius inflation covers the difference between decimal pi endpoints and the true frequency torus.

## Files

- `stage12_replay/N31_predicates.csv`: one result for every archived obligation; Arb ball values include output-rounding radii.
- `stage12_replay/N31_refinement_leaves.csv`: all refined child coordinates and exclusion predicates.
- `stage12_replay/N31_summary.json`: machine-readable counts and versions.
- `stage12_replay/coverage_audit.json`: independent union-cover verification.
- `stage12_full_arb_replay.py`, `stage12_arb_backend.py`, `stage12_cover_audit.py`: replay implementation.

This supports “Both complete finite-spacing certificates were independently replayed using two interval-arithmetic implementations,” with the local-refinement qualification stated in the paper.
