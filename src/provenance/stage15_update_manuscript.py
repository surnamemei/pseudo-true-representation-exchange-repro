"""Apply the effective-threshold wording without changing frozen certificates."""
from pathlib import Path

root=Path(__file__).resolve().parent
body=root/'paper/sections/body_reviewfriendly.tex'
s=body.read_text(encoding='utf-8')
replacements={
 'every sufficiently large odd $N$, with $\\lambda_N=': 'every odd $N\\ge35{,}377$, with $\\lambda_N=',
 'We do not infer a connected finite-spacing bridge or an explicit large-$N$ threshold.': 'The small-spacing radius remains $N$-dependent; no connected finite-spacing bridge is inferred.',
 'a certified continuum crossing and a uniform-convergence theorem for every sufficiently large odd $N$': 'a certified continuum crossing and an interval-assisted theorem for every odd $N\\ge35{,}377$',
 'The small-$z$ radius and large-$N$ threshold are existential;': 'The small-$z$ radius is existential;',
 'The continuum theorem shows that the mechanism persists for every sufficiently large odd record, although its proof does not give an explicit threshold.': 'The continuum theorem shows that the mechanism persists for every odd record with $N\\ge35{,}377$; this is a conservative sufficient threshold.',
 'The continuum argument proves these conditions for every sufficiently large odd $N$, for equal strong amplitudes, $b=10$, and phase $\\pi$, but gives neither an explicit $N_0$ nor a uniform small-$z$ radius.': 'The uniform interval certificate proves these conditions for every odd $N\\ge35{,}377$, for equal strong amplitudes, $b=10$, and phase $\\pi$; it gives no uniform small-$z$ radius.',
 'a full-line continuum certificate plus uniform-convergence argument extends the small-spacing theorem to every sufficiently large odd record.': 'a full-line continuum certificate plus uniform interval transfer extends the small-spacing theorem to every odd $N\\ge35{,}377$.',
 'An explicit small-spacing radius, large-$N$ threshold, and connected finite-spacing bridge remain open.': 'An explicit small-spacing radius and connected finite-spacing bridge remain open.'
}
for a,b in replacements.items():
    if a not in s:raise ValueError(f'Missing expected manuscript phrase: {a}')
    s=s.replace(a,b)
anchor='The centered midpoint kernel fixes the first two large-record corrections to the crossing coefficient.'
if anchor not in s:raise ValueError('Missing engineering-law anchor')
s=s.replace(anchor,anchor+' For this specified family and sufficiently small $z=N\\Delta$, the leading practical law is $\\epsilon_c\\approx0.06650697039(N\\Delta)^2$, with an $O(N^{-2})$ coefficient correction.')
body.write_text(s,encoding='utf-8')

table=root/'paper/sections/evidence_table.tex'
s=table.read_text(encoding='utf-8')
s=s.replace('every sufficiently large odd $N$', 'every odd $N\\ge35{,}377$')
s=s.replace('All sufficiently large odd $N$', 'Odd $N\\ge35{,}377$')
s=s.replace('uniform finite-record transfer', 'uniform directed interval transfer')
table.write_text(s,encoding='utf-8')
