# Stage 18: quartic Hermitian cross-term audit

For a complex residual expansion (r(p)=p r_0+p^2r_1+O(p^3)),

\[
\|r(p)\|_2^2=p^2\|r_0\|_2^2+2p^3\operatorname{Re}\sum_t\overline{r_{0,t}}r_{1,t}+O(p^4).
\]

The active manuscript, `paper/sections/body_reviewfriendly.tex` at equation `eq:review-r1`, and `supplement/S1_complete_local_proof.tex` both write (T_j=2\operatorname{Re}\sum_t\overline{r_{0,j,t}}r_{1,j,t}). The coefficient-producing code, `three_to_two_tone_stage3_tangent.py`, function `cubic_cost`, evaluates `2*mp.re(sum(mp.conj(r0[k])*r1[k] for k in range(N)))`. Thus the manuscript and code use the same Hermitian convention. No quartic coefficient needs recomputation on this ground.

This audit addresses conjugation only. It does not independently validate the entire local expansion, interval certificates, or the unquantified (O(z^6)) radius.
