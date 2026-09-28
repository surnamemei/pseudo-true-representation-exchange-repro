"""Assemble the finite-width N=21 global switch certificate and obligations."""
import csv
import hashlib
import json
import sys
from collections import Counter,defaultdict
from decimal import Decimal,ROUND_FLOOR,ROUND_CEILING
from pathlib import Path

R=Path(__file__).resolve().parent
radius=Decimal(sys.argv[1]) if len(sys.argv)>1 else Decimal('0.0025')
suffix='' if radius==Decimal('0.0025') else '_'+str(radius).replace('.','p')
outer_suffix=str(radius).replace('.','p')
outer_file=f'stage17_epsilon_outer_refinement_{outer_suffix}_14287.json'
inner_file=f'stage17_epsilon_inner_tiled{suffix}.json'
root_file=f'stage17_epsilon_root_tiling{suffix}.json'
joint_file=f'stage17_epsilon_joints{suffix}.json'
cover_file=f'stage17_epsilon_coverage_audit{suffix}.json'
outer=json.loads((R/outer_file).read_text())
inner=json.loads((R/inner_file).read_text())
roots=json.loads((R/root_file).read_text())
joints=json.loads((R/joint_file).read_text())
cover=json.loads((R/cover_file).read_text())
cross=json.loads((R/'stage2_crossing_interval.json').read_text())
assert outer['complete_outer'] and inner['complete'] and roots['all_pass']
assert joints['all_joints_pass'] and cover['all_exact_cover_identities_pass']
assert cross['bracket_certified'] and cross['endpoints'][0]['sign']=='negative'
assert cross['endpoints'][1]['sign']=='positive'
assert Decimal(joints['minimum_gram_determinant_lower'])>0
def parse_interval_point(s):return Decimal(s.strip('[]').split(',')[0].strip())
assert parse_interval_point(roots['minimum_slope_lower'])>0
final=cover['final_closed_interval']
slope_display=str(parse_interval_point(roots['minimum_slope_lower']).quantize(Decimal('0.001'),rounding=ROUND_FLOOR))
gram_display=str(Decimal(joints['minimum_gram_determinant_lower']).quantize(Decimal('0.001'),rounding=ROUND_FLOOR))
cross_lo_display=str(Decimal(cross['epsilon_lower']).quantize(Decimal('0.000000000000001'),rounding=ROUND_FLOOR))
cross_hi_display=str(Decimal(cross['epsilon_upper']).quantize(Decimal('0.000000000000001'),rounding=ROUND_CEILING))
files=['stage17_epsilon_quadratic_trial.py','stage17_epsilon_refine.py',
       'stage17_parametric_jet.py','stage17_epsilon_root_tiling.py',
       'stage17_epsilon_inner_tiled.py','stage17_epsilon_joints.py',
       'stage17_epsilon_coverage_audit.py',
       outer_file,root_file,inner_file,joint_file,cover_file,
       'stage2_crossing_interval.json','interval_boxes.csv','stage2_inner_cells.csv']
def sha(p):return hashlib.sha256((R/p).read_bytes()).hexdigest()
cert=dict(N=21,z=2,b=10,phase='pi',working_radius=str(radius),
          final_closed_epsilon_interval=final,certified_width=cover['final_width'],
          crossing_bracket=[cross['epsilon_lower'],cross['epsilon_upper']],
          outer=dict(archived_parents=outer['parents'],refined_parents=cover['outer_refined_parent_count'],
                     terminal_leaves=outer['leaves'],unresolved=len(outer['unresolved'])),
          inner=dict(archived_parents=inner['archived_terminal_rows'],
                     visited_nodes=inner['visited'],terminal_leaves=inner['leaves'],
                     unresolved=len(inner['unresolved']),max_spatial_depth=inner['max_u_depth'],
                     predicate_counts=inner['predicate_counts']),
          roots=dict(epsilon_slabs=roots['slabs'],
                     branch_krawczyk_checks=2*roots['slabs'],
                     all_pass=roots['all_pass'],
                     minimum_slope_lower=roots['minimum_slope_lower'],
                     shared_boundary_checks=joints['checks'],
                     minimum_shared_box_containment=joints['minimum_box_containment'],
                     minimum_gram_determinant_lower=joints['minimum_gram_determinant_lower']),
          exact_cover_audit=cover,
          cost_difference_endpoint_signs=['negative','positive'],
          unique_transverse_crossing=True,
          A_uniquely_global_below_crossing=True,
          B_uniquely_global_above_crossing=True,
          exactly_A_and_B_global_at_crossing=True,
          no_competitive_third_or_coalescent_fit=True,
          full_global_interval_certified=True,
          scope='one N=21,z=2 interval only; no connected bridge to small-z theorem',
          sha256={p:sha(p) for p in files})
(R/f'finite_width_epsilon_certificate{suffix}.json').write_text(json.dumps(cert,indent=2))

obligations=[]
out_by=defaultdict(list)
for leaf in outer['refined_leaves']:out_by[leaf['parent_index']].append(leaf)
for k in range(outer['parents']):
    x=out_by[k]
    obligations.append(dict(scope='outer',index=k,branch='',epsilon_lo=outer['epsilon_interval'][0],
        epsilon_hi=outer['epsilon_interval'][1],archived_parents=1,terminal_leaves=len(x),
        visited_nodes=2*len(x)-1,unresolved=0,
        minimum_margin=min((v['margin_lower'] for v in x),key=parse_interval_point),
        predicate_summary=';'.join(sorted(set(v['branch_incumbent'] for v in x)))))
in_by=defaultdict(list)
for leaf in inner['records']:in_by[leaf['parent_index']].append(leaf)
for k in range(inner['archived_terminal_rows']):
    x=in_by[k];cnt=Counter(v['predicate'] for v in x)
    margin=[parse_interval_point(v['margin_lower']) for v in x if v['margin_lower']]
    obligations.append(dict(scope='inner',index=k,branch=x[0]['branch'],
        epsilon_lo=roots['interval'][0],epsilon_hi=roots['interval'][1],
        archived_parents=1,terminal_leaves=len(x),visited_nodes=2*len(x)-1,unresolved=0,
        minimum_margin=str(min(margin)) if margin else '',
        predicate_summary=';'.join(f'{p}:{n}' for p,n in sorted(cnt.items()))))
for row in roots['records']:
    for w in ('A','B'):
        b=row['branches'][w]
        obligations.append(dict(scope='root_slab',index=row['index'],branch=w,
            epsilon_lo=row['epsilon_interval'][0],epsilon_hi=row['epsilon_interval'][1],
            archived_parents=0,terminal_leaves=1,visited_nodes=1,unresolved=0,
            minimum_margin=row['slope_lower'],
            predicate_summary=f"Krawczyk:{b['inclusion']};Hessian:{b['positive_Hessian']}"))
for row in joints['records']:
    obligations.append(dict(scope='joint',index=row['left_slab'],branch=row['branch'],
        epsilon_lo=row['epsilon'],epsilon_hi=row['epsilon'],archived_parents=0,
        terminal_leaves=1,visited_nodes=1,unresolved=0,
        minimum_margin=row['min_box_containment'],predicate_summary='common_root_in_both_boxes'))
with (R/f'epsilon_interval_obligations{suffix}.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=list(obligations[0]));w.writeheader();w.writerows(obligations)

md=rf"""# Finite-width global epsilon certificate at N=21, z=2

**Certified closed interval:** \([{final[0]},\,{final[1]}]\), width \({cover['final_width']}\). It is contained strictly inside both the outer and root/inner working intervals, each approximately \(\epsilon_c\pm{radius}\). The unique crossing lies in the outward-rounded bracket \([{cross_lo_display},\,{cross_hi_display}]\); the exact directed endpoint-sign bracket is in the JSON certificate.

For every epsilon in the displayed interval, A and B are separated strict local minima with positive Gram determinants and finite fitted amplitudes. They are the only globally competitive two-tone fits. The optimized cost difference has positive derivative (uniform lower bound \(>{slope_display}\)), so A is uniquely global below the single crossing, B above it, and both tie there. No third regular or coalescent fit is competitive.

## Exact quadratic structure

For a fixed fitted subspace \(P\), write \(x(\epsilon)=x_0+\epsilon x_1\). Then
\[
J(P;\epsilon)=\|(I-P)x(\epsilon)\|^2
=a(P)\epsilon^2+b(P)\epsilon+c(P),
\]
where \(a=\|(I-P)x_1\|^2\), \(b=2\operatorname{{Re}}\langle(I-P)x_0,(I-P)x_1\rangle\), and \(c=\|(I-P)x_0\|^2\). The code encloses these coefficients with directed intervals. For each outer spatial cell, the centered confluent basis gives a quadratic \(J_c(\epsilon)\), while fixed A/B fitted subspaces give feasible quadratic incumbents \(U_A,U_B\). The spatial projector displacement \(\theta\) is independent of epsilon. On each epsilon interval, the code bounds the minimum of the quadratic difference \(J_c-U_w\) and the maximum of \(2\theta\sqrt{{J_c E}}-\theta^2E\), after verifying \(\sqrt{{J_c}}-\theta\sqrt E>0\). Their strict inequality excludes the entire cell.

Inside the two accepted neighborhoods, Taylor lower bounds, gradients, and Krawczyk exclusion are also evaluated as coefficient-wise epsilon quadratics. Epsilon is subdivided only where a spatial cell needs it. Each root is enclosed by a parameterized Krawczyk map on {roots['slabs']:,} adjacent epsilon slabs; {joints['checks']:,} endpoint checks certify that neighboring slabs share the same root. The Gram determinant stays above \({gram_display}\), and the slope stays above \({slope_display}\).

## Coverage accounting

The 14,287 archived outer exclusion parents become {outer['leaves']:,} terminal leaves after refinement of {cover['outer_refined_parent_count']:,} parents. The 819 archived crossing inner parents generate {inner['leaves']:,} terminal three-dimensional leaves from {inner['visited']:,} visited nodes. All leaves pass; none is unresolved. The exact-decimal audit checks containment and volume identity for every archived parent, and checks that all {roots['slabs']:,} epsilon slabs meet exactly. The original directed bracket has opposite certified endpoint signs. The per-parent and per-slab record is [epsilon_interval_obligations{suffix}.csv](epsilon_interval_obligations{suffix}.csv); the full inner leaves and root checks are in the Stage-17 machine logs.

## Limit

This certifies a finite epsilon interval for \(N=21,z=2\). It does not prove a connected \(z\)-continuation from the small-spacing theorem. {"The separately certified \\(\\epsilon_c\\pm0.035\\) samples lie outside this interval." if radius<Decimal('0.035') else "The separately certified outer amplitude samples are distinct records; the displayed inward-rounded closed interval is the continuous claim."} The endpoints are sufficient, not maximal.
"""
(R/f'finite_width_epsilon_certificate{suffix}.md').write_text(md,encoding='utf-8')
print(json.dumps({k:v for k,v in cert.items() if k not in ('sha256','exact_cover_audit')},indent=2))
