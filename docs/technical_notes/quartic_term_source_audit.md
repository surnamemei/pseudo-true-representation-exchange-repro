# Quartic Hermitian cross-term source audit

**PASS.** The active manuscript, active mathematical supplement, high-precision candidate code, and directed-interval code use the same Hermitian cross term. No implementation disagreement was found.

| Source | Inspected expression | Result |
|---|---|---|
| `paper/sections/body_reviewfriendly.tex`, Eq. `review-r1` | `T_j=2\Re\sum_t\overline{r_{0,j,t}}r_{1,j,t}` | Unambiguously conjugates the leading residual. |
| `supplement/S1_complete_local_proof.tex`, global theorem/quartic subsection | `T_j=2\Re\sum\overline{r_{0,j}}r_{1,j}` | Same Hermitian ordering. |
| `three_to_two_tone_stage3_tangent.py`, lines 44--45 | `r1=[u**4/12-1j*kappa*u*alpha+kappa*kappa*u*u ...]`; `2*mp.re(sum(mp.conj(r0[k])*r1[k] ...))` | Direct complex implementation. |
| `three_to_two_tone_stage3_tangent_interval.py`, `cubic_envelope` | `acc+=2*(rre*wre+rim*wim)` | Equals `2 Re(conj(r0)*r1)` after real/imaginary splitting. |

At weak phase `phi=pi`, the code uses real `alpha`, `beta`, and `kappa`. It sets `r0 = rre+i*rim` and `r1=wre+i*wim`, with `wre=s^4/12+kappa^2*s^2` and `wim=-kappa*alpha*s`. Then `Re(conj(r0)*r1)=rre*wre+rim*wim`; the plus sign on the imaginary product is correct. The expression for `r1` follows by expanding `(2+alpha*p)exp(i*kappa*p*s)+beta*p*exp(i*v*s)` through `p^2`, where `p=z^2`.

The code's candidate and interval formulas therefore implement the same quartic coefficient used in the paper. This is a source-equivalence audit; it does not add a new theorem or re-run the already recorded interval certificate.
