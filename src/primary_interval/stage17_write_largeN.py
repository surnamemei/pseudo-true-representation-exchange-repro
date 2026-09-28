"""Assemble the composite uniform N>=10001 certificate and margin table."""
import csv
import hashlib
import json
from pathlib import Path
from fractions import Fraction
from mpmath import mp
from stage17_kernel_bound import upper_error

R=Path(__file__).resolve().parent
old=json.loads((R/'explicit_N0_certificate.json').read_text())
trial=json.loads((R/'stage17_interval_trial_N10001.json').read_text())
ref=json.loads((R/'stage17_refinement_N10001.json').read_text())
margin=json.loads((R/'stage17_largeN_margins.json').read_text())
assert trial['failed_total']==4 and ref['complete']
assert margin['original_cells']==27708 and margin['replaced_parents']==4
assert all(mp.mpf(x['margin_lower'])>0 for x in margin['predicate_margins'].values())
assert trial['local']['inclusion'] and trial['local']['signs'] and trial['root_boxes_linked']
assert mp.mpf(trial['coalescent_margin_lower'])>0 and mp.mpf(trial['tail_margin_lower'])>0
majorant=upper_error(10001)
assert majorant<Fraction(1,10**45)
def sha(p):return hashlib.sha256((R/p).read_bytes()).hexdigest()
continuum={
 'C2 curvature': '0.000133682',
 'C3 satellite coefficient': '0.0644229455',
 'C5 transversality': '0.1078175136',
 'C6 regular cost': '0.000000006281428732256794875',
 'C6 derivative sign': '0.000000001339436331021652144',
 'C6 monotone cell': '0.000000003793579161565239847',
 'C7 coalescent': '0.0008511681673580126948230596234041332452',
 'C8 finite-period tail': '0.0037492892323564881493884541052232812113'
}
finite={
 'C2 curvature':trial['local']['curvature_A_lower'],
 'C3 satellite coefficient':trial['local']['beta_B_abs_lower'],
 'C5 transversality':trial['local']['slope_lower'],
 'C6 regular cost':margin['predicate_margins']['cost_excluded']['margin_lower'],
 'C6 derivative sign':margin['predicate_margins']['gradient_excluded']['margin_lower'],
 'C6 monotone cell':margin['predicate_margins']['monotone_derivative_no_root']['margin_lower'],
 'C7 coalescent':trial['coalescent_margin_lower'],
 'C8 finite-period tail':trial['tail_margin_lower']
}
rows=[]
for key in continuum:
    c=mp.mpf(continuum[key]);f=mp.mpf(finite[key])
    rows.append(dict(condition=key,continuum_margin=continuum[key],
                     finite_uniform_lower_at_N10001=finite[key],
                     erosion_upper=str(max(mp.mpf(0),c-f)),
                     retained_fraction=str(f/c),
                     analytic_kernel_tail_upper_at_N10001='4.846e-46',
                     sufficient_joint_odd_N_floor=10001,
                     condition_specific_minimum='not optimized'))
rows.insert(0,dict(condition='C1/C4 joint root inclusion and equality',
                   continuum_margin='strict Krawczyk inclusion',
                   finite_uniform_lower_at_N10001='inclusion true; contraction upper 0.142794',
                   erosion_upper='not scalar',retained_fraction='not scalar',
                   analytic_kernel_tail_upper_at_N10001='4.846e-46',
                   sufficient_joint_odd_N_floor=10001,
                   condition_specific_minimum='not optimized'))
rows.insert(1,dict(condition='C1 tangent Gram factors',
                   continuum_margin='strictly positive',
                   finite_uniform_lower_at_N10001='A > 0.23505228; B > 0.92073167',
                   erosion_upper='not scalar',retained_fraction='not scalar',
                   analytic_kernel_tail_upper_at_N10001='4.846e-46',
                   sufficient_joint_odd_N_floor=10001,
                   condition_specific_minimum='not optimized'))
with (R/'largeN_condition_thresholds.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
files=['stage17_uniform_trial.py','stage17_kernel_bound.py','stage17_refine_uniform.py',
       'stage17_largeN_margin_summary.py','stage14_continuum_global_partition.csv',
       'stage17_interval_trial_N10001.json','stage17_refinement_N10001.json',
       'stage17_largeN_margins.json']
cert=dict(old_sufficient_odd_N_floor=35377,new_sufficient_odd_N_floor=10001,
          improvement_factor='3.5373462653734627',
          scope='all odd N>=10001, b=10, phi=pi, small-z tangent criterion C1-C8',
          parameter_interval=['0','1/10001^2'],
          kernel_tail_exact_rational=str(majorant),
          kernel_tail_decimal_upper='4.846e-46',
          used_derivative_enclosure='1e-45 through derivative order 46 on |q|<=110',
          near_zero_Taylor_enclosure='1e-25 for deflated A/v^4,B/v^2 and first two derivatives on |v|<=1',
          frozen_partition_cells=27708,failed_frozen_predicates=4,
          refined_parent_cells=4,refined_terminal_leaves=6,
          root_linked=trial['root_boxes_linked'],local=trial['local'],
          predicate_margins=margin['predicate_margins'],
          coalescent_margin_lower=trial['coalescent_margin_lower'],
          finite_period_tail_margin_lower=trial['tail_margin_lower'],
          all_uniform_conditions_pass=True,
          lower_trials={'N1501':'local and at least 30 frozen cells fail; not a theorem counterexample',
                        'N5001':'19 frozen cells fail and original root box is not linked; not a theorem counterexample'},
          minimality_claim=False,
          sha256={p:sha(p) for p in files})
(R/'tightened_N0_certificate.json').write_text(json.dumps(cert,indent=2))
print(json.dumps({k:v for k,v in cert.items() if k not in ('local','predicate_margins','sha256')},indent=2))
