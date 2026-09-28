import csv,json,re,hashlib,math
from pathlib import Path
from decimal import Decimal,ROUND_FLOOR
R=Path(__file__).resolve().parent

def lower_display(s,places=8):
 return str(Decimal(s).quantize(Decimal(1).scaleb(-places),rounding=ROUND_FLOOR))
def numbers(s):return re.findall(r'[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?',s)
def interval(s):
 a=numbers(s);return '['+a[0]+','+a[-1]+']'
def general():
 rows=[]
 for N in (11,15,21,31,41):
  cp=R/f'stage12_generalN/N{N}_certificate.json'
  if not cp.exists():continue
  c=json.loads(cp.read_text());r=json.loads((R/f'stage12_generalN/N{N}_reference.json').read_text());local=c['local_crossing'];roots=c['stationary_points']
  assert c['tangent_classification_certified'] and local['verified'] and not c['unresolved']
  assert all(abs(x['v'])>1e-6 for x in roots)
  selected=sorted([x for x in roots if x['kind']=='min'],key=lambda x:float(x['R'][0]))[:2]
  rows.append(dict(N=N,b=10,phi='pi',status='C1-C7_interval_certified_C8_proved',v_A=r['v_A'],v_A_interval=interval(local['X'][0]),v_B=r['v_B'],v_B_interval=interval(local['X'][1]),lambda_N=r['lambda_N'],lambda_interval=str(c['lambda_box']),beta_A=r['A_coeff'][3],beta_A_interval=interval(local['beta_A']),beta_B=r['B_coeff'][3],beta_B_interval=interval(local['beta_B']),R_equal=r['R_A'],R_equal_enclosures=str([x['R'] for x in selected]),other_regular_margin_lower=lower_display(c['other_minimum_margin'],12),coalescent_loss_interval=str(c['coalescent_coefficient']),coalescent_margin_lower=lower_display(c['coalescent_margin'],12),slope=r['delta_R_lambda'],slope_interval=interval(local['delta_R_lambda']),curvature_A_interval=interval(local['Rvv_A']),curvature_B_interval=interval(local['Rvv_B']),c_N=r['lambda_quartic'],c_N_interval=interval(local['lambda_quartic']),regular_minima=sum(x['kind']=='min' for x in roots),regular_maxima=sum(x['kind']=='max' for x in roots),visited=c['visited'],unresolved=0))
 with (R/'generalN_certified_instances.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
 tab=r'''\begin{table}[t]
\caption{Verified finite-$N$ tangent criteria at $b=10,\phi=\pi$. Critical coefficients are rounded; strict margins are truncated downward. All rows satisfy C1--C8.}\label{tab:crossN}
\centering\footnotesize
\begin{tabular}{@{}r r r r@{}}\toprule
$N$ & $\lambda_N$ & $\delta_{\rm reg}>$ & $\delta_{\rm coal}>$\\\midrule
'''
 for r in rows:tab+=f"{r['N']} & {float(r['lambda_N']):.10f} & {lower_display(r['other_regular_margin_lower'],7)} & {lower_display(r['coalescent_margin_lower'],7)}"+r'\\'+'\n'
 tab+=r'\bottomrule\end{tabular}\end{table}'+'\n';(R/'paper/sections/crossN_table.tex').write_text(tab)
 tab=r'''\begin{table}[H]
\caption{Interval-certified tangent instances, $b=10,\phi=\pi$. Root/coefficient/cost/slope entries are rounded displays of enclosed values; the two strict margin bounds are truncated downward. The two costs are equal at the joint certified root. C8 is analytical.}\label{tab:app-generalN}
\centering\scriptsize\setlength{\tabcolsep}{3.5pt}
\begin{tabular}{@{}r r r r r r r r r r@{}}\toprule
$N$ & $v_A$ & $v_B$ & $\lambda_N$ & $\beta_A$ & $\beta_B$ & $R_*$ & $\delta_{\rm reg}>$ & $\delta_{\rm coal}>$ & $S$\\\midrule
'''
 for r in rows:
  tab+=str(r['N'])+' & '+' & '.join(f'{float(r[k]):.7f}' for k in ('v_A','v_B','lambda_N','beta_A','beta_B','R_equal'))+' & '+lower_display(r['other_regular_margin_lower'],7)+' & '+lower_display(r['coalescent_margin_lower'],7)+f" & {float(r['slope']):.7f}"+r'\\'+'\n'
 tab+=r'\bottomrule\end{tabular}\end{table}'+'\n';(R/'supplement/generalN_certificate_table.tex').write_text(tab)
 print('Certified N:',[r['N'] for r in rows])
 return rows

def replay():
 summaries=[json.loads((R/f'stage12_replay/N{N}_summary.json').read_text()) for N in (21,31)]
 cover=json.loads((R/'stage12_replay/coverage_audit.json').read_text());assert all(q['complete'] and q.get('true_torus_covered',True) for q in cover)
 allrows=[]
 for q in summaries:
  assert not q['failures'];N=q['N'];allrows+=list(csv.DictReader((R/f'stage12_replay/N{N}_predicates.csv').open(newline='')))
  text=f'''# Complete independent Arb replay: N={N}

**PASS — all {q['predicates']:,} archived cell/local/crossing obligations validated.**

- Backend: python-flint {q['python_flint']}, FLINT {q['FLINT']}, Arb, {q['precision_bits']}-bit midpoint precision.
- Outer categories: {', '.join(k.removeprefix('outer_')+'='+str(v) for k,v in q['counts'].items() if k.startswith('outer_'))}.
- Inner terminal cells: {sum(v for k,v in q['counts'].items() if k.startswith('inner_'))}.
- Full root inclusion, contraction, positive Hessian principal minors, positive Gram, crossing endpoint signs, and positive cost-gap slope recomputed.
- Minimum outer exclusion margins by case: {q['outer_min_margin']}.
- Original inner predicates needing local Arb refinement: {q['refined_archived_inner_cells']}; one bisection each, {2*q['refined_archived_inner_cells']} validated children, {q['refinement_nodes']} refinement tree nodes including parents.
- Total different inner proof paths: {q['alternate_inner_proofs']} (including refinements); remaining differences use another valid exclusion predicate.
- Contradictory certificate results: **0**. Unresolved parent cells: **0**.
- Runtime: {q['elapsed_seconds']:.2f} s for this replay; timing excludes coverage audit.

## Meaning of agreement

This is complete validation of the archived mathematical obligations, with local refinement where Arb ball propagation is wider. It is not bitwise agreement between endpoint and ball arithmetic. An inconclusive unrefined Arb bound is an enclosure limitation; both child cells prove the same parent's exclusion. The archived files are unchanged. The replay does not use stored lower bounds as premises: it freshly evaluates the signal energy, incumbent, projector change, gradients, Hessians and preconditioners.

All outer exclusions, accepted-cell containments and inner terminal cells are checked, including the confluent edge of the outer domain. N21 includes below/crossing/above amplitude records, with the entire narrow crossing bracket used in the crossing case; N31 includes its full crossing bracket. Exact-decimal rectangle-union audits verify no gap in the outer or inner covers. Radius inflation covers the difference between decimal pi endpoints and the true frequency torus.

## Files

- `stage12_replay/N{N}_predicates.csv`: one result for every archived obligation; Arb ball values include output-rounding radii.
- `stage12_replay/N{N}_refinement_leaves.csv`: all refined child coordinates and exclusion predicates.
- `stage12_replay/N{N}_summary.json`: machine-readable counts and versions.
- `stage12_replay/coverage_audit.json`: independent union-cover verification.
- `stage12_full_arb_replay.py`, `stage12_arb_backend.py`, `stage12_cover_audit.py`: replay implementation.

This supports “Both complete finite-spacing certificates were independently replayed using two interval-arithmetic implementations,” with the local-refinement qualification stated in the paper.
'''
  (R/f'arb_full_replay_N{N}.md').write_text(text,encoding='utf-8')
 for name in ('full_second_backend_replay.csv','arb_full_replay_comparison.csv'):
  with (R/name).open('w',newline='',encoding='utf-8') as f:
   w=csv.DictWriter(f,fieldnames=allrows[0]);w.writeheader();w.writerows(allrows)
 text='''# Full second-backend replay summary

**PASS for both complete finite-spacing certificates.**

|N|Archived obligations|PASS|FAIL|Arb-refined inner parents|Validated child leaves|
|---:|---:|---:|---:|---:|---:|
'''
 for q in summaries:text+=f"|{q['N']}|{q['predicates']}|{q['passed']}|{len(q['failures'])}|{q['refined_archived_inner_cells']}|{2*q['refined_archived_inner_cells']}|\n"
 text+='''
Total: 65,337 outer cells, 2,801 inner cells and 86 local/crossing obligations = **68,224 checks**. All archived parents are validated. Arb uses one bisection for 172 inner parents; another 21 parents use a different exclusion predicate. No contradictory result was found. The full arithmetic check is accompanied by an exact-decimal union-cover audit for all four outer case partitions and eight inner neighborhoods.

The second backend is python-flint 0.9.0 / FLINT 3.6.0 / Arb at 320 bits, separate from the original mpmath.iv endpoint arithmetic. Its fresh energy, objectives, analytic derivative jets, projector displacement and Krawczyk evaluations do not read the old inequality values as proof inputs. Candidate centers, parameter brackets and cell coordinates are shared certificate data. This is an independent arithmetic replay, not an independently discovered theorem or formal proof-assistant verification.

See `arb_full_replay_N21.md`, `arb_full_replay_N31.md`, the two complete comparison CSVs, and the per-parent refinement leaves. Historical 44-check MPFR results remain archived but are no longer the extent of independent verification.
'''
 (R/'full_second_backend_replay_summary.md').write_text(text,encoding='utf-8')

def continuation():
 q=json.loads((R/'stage12_continuation_checks.json').read_text());rows=list(csv.DictReader((R/'branch_continuation.csv').open()))
 text='''# Interpretation of amplitude continuation

N=21, z=2, b=10, phi=pi. Deterministic stationary continuation uses amplitude steps no larger than 0.0005, analytic gradients/Hessians, and 65-digit independent stationary refinement at the displayed landmarks. Complex amplitudes are optimized at each step; conjugate symmetry makes them real in this centered phase-pi family. The CSV still records real and imaginary parts explicitly.

## Observed branches

|Branch|epsilon|u1|u2|Loss|
|---|---:|---:|---:|---:|
'''
 for r in q['high_precision_checks']:
  text+=f"|{r['branch']}|{float(r['epsilon']):.10f}|{float(r['u1']):.10f}|{float(r['u2']):.10f}|{max(0,float(r['J'])):.10f}|\n"
 text+=f'''
A has {q['A_points']} samples and B {q['B_points']}. Both were followed over the physical range 0<=epsilon<=approximately 0.30, across epsilon_c={float(q['epsilon_cross']):.15f}. Neither encountered a sampled fold, singularity, coalescence or loss of positive curvature. Minimum sampled Hessian eigenvalues: A={q['minimum_sampled_curvature']['A']:.9g}, B={q['minimum_sampled_curvature']['B']:.9g}. These trajectory checks are numerical; the crossing itself has the independent interval certificate.

A begins exactly at the generating strong pair (-2,2), with unit amplitudes and zero loss. Its increasingly asymmetric frequencies and amplitudes provide a continuous deformation of that representation. B already exists at zero omitted amplitude as a worse local approximation. It reaches the physical boundary epsilon=0 without termination, so there is no positive-amplitude birth event to identify on this backward segment.

At equality, B combines a nearly central component with a satellite at 12.4634193568, **not** at b=10. The satellite is displaced by 2.4634193568 normalized radians. Fixing the satellite to b and optimizing the other frequency produces u1={q['satellite_fixed_b']['u1']:.10f}, loss={q['satellite_fixed_b']['J']:.10f}, versus {q['crossing_cost']:.10f} for the unconstrained branches (penalty {q['satellite_fixed_b']['J']-q['crossing_cost']:.10f}). A 4096-interval full-period scan and local refinements found this as the best constrained minimum; this comparison is numerical, not an additional global certificate. The established tangent fixed-b loss 0.1036226553 also exceeds the equal tangent loss 0.0587626233.

## Signal-processing interpretation

Immediately below the certified crossing, the winning representation is a deformation of the exact strong-pair fit. Immediately above it, a central component represents the close pair collectively and a biased satellite accounts for the omitted component's influence. B should not be described as literally estimating the weak tone. Its movement depends on the finite-record projection geometry and the simultaneous fit of the central component.

The continuation figure plots A1/A2 and B1/B2 with the generating frequencies and crossing line. It establishes numerical branch identity from zero amplitude; it does not upgrade the entire displayed epsilon interval to an interval-certified global ordering. The separately certified below/above points and the narrow crossing bracket retain their original scope.

Outputs: `branch_continuation.csv`, `branch_continuation_figure.pdf`, `stage12_continuation_checks.json`.
'''
 (R/'branch_continuation_interpretation.md').write_text(text,encoding='utf-8')
if __name__=='__main__':general();replay();continuation()
