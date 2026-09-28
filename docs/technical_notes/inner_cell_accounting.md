# Exact inner-cell accounting

**Result: complete coverage; no missed terminal obligation.** The apparently conflicting numbers count different objects and different cases.

| Scope | Visited tree nodes | Internal nonterminal parents | Archived terminal obligations |
|---|---:|---:|---:|
| N=21 crossing, A/B | 1,636 | 817 | 819 |
| N=31 crossing, A/B | 1,648 | 823 | 825 |
| **Crossing only** | **3,284** | **1,640** | **1,644** |
| N=21 auxiliary epsilon_c-0.035, A/B | 1,676 | 837 | 839 |
| N=21 auxiliary epsilon_c+0.035, A/B | 634 | 316 | 318 |
| **All archived inner searches** | **5,594** | **2,793** | **2,801** |

Each of the eight original inner searches is a full binary subdivision tree. Its identity is **visited = 2 × terminal − 1**. Therefore the 3,284 “Inner A/B” count in the manuscript's main certificate table is the sum of visited nodes for the *crossing cases only*. It includes 1,640 internal parent nodes that were subdivided and are not terminal proof obligations. The Arb report's 2,801 is the count of *archived terminal rows* across N=21 below/crossing/above and N=31 crossing: 1,976 N21 + 825 N31 = 2,801. The N21 auxiliary cases contribute 1,157 of these terminal rows. Thus 3,284 and 2,801 have different denominators.

The archived terminal logs have **zero duplicate coordinate rows within a case/branch**. The 24 accepted root-box terminal rows (16 N21, 8 N31) are included in the 2,801, rather than counted separately. Every archived terminal row has exactly one matching inner Arb predicate row with the same index, case, branch, and original predicate, and every such row passes. The complete neighborhood-coverage audits have exactly the same eight case/branch cell counts and zero uncovered slabs.

Arb's wider balls required refinement of 172 archived terminal parents (128 N21, 44 N31). Each was replaced by exactly two passing child leaves; 344 child records have complete parent-domain coverage audits. Hence the *effective* Arb terminal partition has **2,801 − 172 + 344 = 2,973** leaves. The replay's “refinement nodes” count is 516 = 3 × 172: it counts each failed parent evaluation plus its two children. Arb sometimes used a different valid predicate on an archived terminal row (193 such rows); this changes no coverage count. The N21 auxiliary cases are in the Arb replay, but not in Table III's crossing-only visited count.

The machine-auditable eight-row breakdown is inner_cell_accounting.csv; stage17_inner_accounting.py checks all identities, one-to-one replay matches, duplicate absence, subdivision coverage, and accepted neighborhood coverage. **No fatal replay gap exists.**
