# Stage 5 — fair-matching / normalization audit across N

Evidence level of the numerical table: `GLOBAL_NUMERICAL` (pre-registered crossing analysis, rectangular weights, ε ∈ [0, 1]; reference-level validation at each crossing). The language audit is a reading of the frozen sources; **no manuscript file was edited.**

## Two different asymptotic questions

* **Comparison A — fixed normalized geometry** (the manuscript's convention): z = NΔ = 2 and b = Nω₃ = 10 held fixed, so the physical half-spacing Δ = 2/N and the weak-tone frequency ω₃ = 10/N shrink as the record grows. The strong pair stays at 0.637 DFT bins and the weak tone at 1.59 bins for every N.
* **Comparison B — fixed physical frequencies**: Δ = 2/21 and ω₃ = 10/21 rad/sample held fixed (the N = 21 baseline signal observed over longer or shorter records), so z = 2N/21 and b = 10N/21 grow with N and the pair becomes resolved.

| N | A: z, b | A: ε_c | B: z, b (bins of pair) | B: ε_c | B: status / third-margin (rel.) |
|---|---|---|---|---|---|
| 11 | 2, 10 | 0.24347 | 1.048, 5.238 (0.33 bins) | 0.10790 | present / 0.70 |
| 21 | 2, 10 | 0.24819 | 2, 10 (0.64 bins) | 0.24819 | present / 1.09 |
| 31 | 2, 10 | 0.24914 | 2.952, 14.76 (0.94 bins) | 0.55823 | present / 1.70 |
| 41 | 2, 10 | 0.24949 | 3.905, 19.52 (1.24 bins) | 0.84086 | present / 0.87 |
| 61 | 2, 10 | 0.24974 | 5.810, 29.05 (1.85 bins) | 0.95972 | **inconclusive**: third competitor within 0.13% of J* |

Figure: `stage5_matching/matching.png|pdf`. Table: `stage5_matching/matching_table.csv`.

**Paragraph for the authors.** Under the manuscript's matching (Comparison A) the finite-spacing crossing converges smoothly as N grows (ε_c = 0.2435 → 0.2497 for N = 11 → 61 at z = 2), consistent with the continuum/large-N theorem, whose λ_∞ and "all odd N ≥ 10,001" statements are asymptotics **at fixed z = NΔ and b = Nω₃**, i.e. for strong tones whose physical separation shrinks like 1/N. Holding the *physical* frequencies fixed instead (Comparison B) is a different question: ε_c moves by a factor ≈ 9 between N = 11 and N = 61, and by N = 61 (pair separation 1.85 bins, weak tone at 9.2 bins) the global switch is between "both strong tones" (−5.66, 6.05) and "one strong tone + the weak tone" (−6.78, 28.88), with the mirror choice (the other strong tone + weak tone) nearly tied — the ordinary resolved-tone model-order competition, not the close-pair mechanism. None of the existing theorems addresses Comparison B, and none should be read as doing so.

## Manuscript sentences that could be misread as fixed-physical-frequency statements

Risk levels: **H** = likely misread by a signal-processing reader skimming the abstract/intro; **M** = clear in context but ambiguous in isolation; **L** = explicit enough.

| # | Location (frozen source) | Sentence (verbatim excerpt) | Risk | Why |
|---|---|---|---|---|
| 1 | Abstract (`body_reviewfriendly.tex`) | "…a continuum certificate transfers the result to all odd $N\ge10{,}001$." | **H** | "the result" reads like "the exchange for a given signal persists for all long records"; the transfer is at fixed z = NΔ, b = Nω₃ (and small z). |
| 2 | Abstract | "for the specified $b=10$ family, $\lambda_\infty\approx0.066507$ as $N\to\infty$" | M | b = 10 is normalized; the reader must recall ω₃ = b/N to see the weak tone moves toward DC as N grows. |
| 3 | Introduction | "For the specified $b=10$, weak-tone phase $\phi=\pi$ family and large odd records, $\epsilon_c\sim0.06650697039(N\Delta)^2$" | M | (NΔ)² makes the scaling explicit, but "large odd records" invites a fixed-signal reading. |
| 4 | Introduction | "…certify its continuum limit and transfer to all sufficiently large odd lengths." | **H** | "lengths" suggests record length of a fixed signal. |
| 5 | Related work | "…transfer to all sufficiently large odd records, and two separate finite-spacing certificates." | M | as #4. |
| 6 | Theorem 2 title/statement (`continuum_stage14.tex`) | "Continuum crossing and all large odd $N$" / "For every odd $N\ge10{,}001$, conditions C1--C8 hold…" | M | Precise (C1–C8 are defined in normalized coordinates) but the title alone reads as a record-length persistence statement. |
| 7 | Section heading (`continuum_stage14.tex`) | "Continuum crossing and large-record persistence" | **H** | "large-record persistence" is exactly the fixed-physical-frequency reading that Stage 5 shows to be false in general. |
| 8 | Supplement S1b | "Continuum certificate and large-record transfer"; "Effective large-record threshold" | M | as #7, lower visibility. |
| 9 | Discussion | "The continuum theorem supplies an explicit conservative sufficient threshold for all large odd records" | M | as #4. |
| 10 | Discussion | "The centered midpoint kernel fixes the first two large-record corrections." | L | technical; "corrections" to λ_N in N at fixed normalized geometry. |
| 11 | Limitations | "…the uniform interval certificate covers every odd length above its stated sufficient threshold for equal, record-center phase-aligned strong amplitudes, $b=10$, and weak phase $\pi$." | M | lists every other scope condition but not "at fixed z = NΔ, b = Nω₃". |
| 12 | Conclusion | "…extends the small-spacing theorem to all sufficiently large odd lengths." | **H** | last sentence a reviewer reads; same issue. |
| 13 | Formulation | "Let $\Delta>0$ be the true physical strong-tone half-spacing, with … $\omega_3=b/N$." | L | this is the one place the convention is explicit; it should be echoed wherever N → ∞ is invoked. |

**Recommendation (for a later revision; nothing edited now):** wherever "all odd N ≥ 10,001", "large records", "N → ∞" or "large-record persistence" appears, add "at fixed normalized geometry z = NΔ and b = Nω₃ (strong-pair separation and weak-tone location fixed in DFT bins)", and add one sentence stating that holding the physical frequencies fixed while N grows resolves the pair and is not covered (Stage 5, Comparison B).
