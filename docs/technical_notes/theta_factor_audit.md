# Centered-projector factor audit

**Result: PASS.** The mean-frequency term is `sqrt(2) * L * r_m`, where `L=(N-1)/2`. No active certificate implementation uses `sqrt(2*L) * r_m`.

| Location | Expression | Interpretation |
|---|---|---|
| `paper/sections/body_reviewfriendly.tex`, full-domain exclusion paragraph | `\sqrt2 Lr_m` | Main-text mean-frequency term, with `L=(N-1)/2`; complete formula is in the supplement. |
| `supplement/S2_interval_certificate_protocol.tex`, Eq. `app-theta` | `\sqrt2 Lr_m` | Derived from `||\partial_m Q||_F <= \sqrt2 L`. |
| `three_to_two_tone_stage2_outer_global.py`, original N21 backend | `iv.sqrt(2)*10*rmi` | `10=(21-1)/2`. |
| `stage7_n31_certificate.py`, N31 backend | `iv.sqrt(2)*tmax*rmi` | `tmax=(N-1)//2`. |
| `stage17_epsilon_quadratic_trial.py`, finite-width N21 backend | `iv.sqrt(2)*10*rm` | Same N21 factor over amplitude slabs. |
| `stage12_full_arb_replay.py`, independent Arb replay | `I(2).sqrt()*tm*rm` | `tm=(N-1)//2`; separately recomputed projector bound. |

The distinction matters: replacing `sqrt(2)*L` with `sqrt(2*L)` would make the mean-frequency contribution too small for `L>1` and could invalidate outer exclusions. The inspected original and replay paths use the larger, correct factor. This audit checks the source expression and its parameter definitions; the existing interval replay checks the resulting predicates.
