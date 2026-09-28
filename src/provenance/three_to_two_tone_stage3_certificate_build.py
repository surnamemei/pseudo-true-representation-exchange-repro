"""Machine-readable index of the analytical and interval Stage-3 claims."""
import csv,hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parent
tan=json.loads((ROOT/'stage3_tangent_reference.json').read_text())
iv=json.loads((ROOT/'stage3_tangent_interval.json').read_text())
stage2=json.loads((ROOT/'validated_globality_certificate.json').read_text())
rows=list(csv.DictReader((ROOT/'lambdaN_theory_vs_certified.csv').open(newline='')))
scripts=['three_to_two_tone_stage3_tangent.py','three_to_two_tone_stage3_tangent_interval.py',
         'three_to_two_tone_stage3_full_crossings.py','three_to_two_tone_stage2_interval.py']
hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in scripts}
out=dict(model=dict(N=21,strong_amplitudes=[1,1],strong_frequencies_u=['-d','d'],
                    weak_frequency_u=10,weak_phase='pi',d_definition='N*Delta',
                    fit='two free complex-amplitude undamped tones'),
         local_small_d_theorem=dict(status='proved_with_existence_constants',
             branch_A_satellite_u=tan['v_A'],branch_B_satellite_u=tan['v_B'],
             lambda_quadratic=tan['lambda_N'],lambda_quartic=tan['lambda_quartic'],
             expansion='epsilon_c=lambda_quadratic*d^2+lambda_quartic*d^4+O(d^6)',
             local_transversality=tan['delta_R_lambda'],
             noncoalescence='|nu_s-nu_c|>=min(|v_A|,|v_B|)/(2*N) for sufficiently small d',
             interval_conditions=iv,
             explicit_d0=None,explicit_remainder_C=None,
             explicit_radius_status='not_computed; existence follows from analytic IFT'),
         global_baseline=dict(status='interval_checked',d=2,
             crossing=stage2['crossing'],
             full_crossing_bracket_global=stage2['crossing_global_bracket']['full_bracket_global'],
             outer_and_inner_case_status=[z['global_result'] for z in stage2['cases']]),
         open_region=dict(local_asymptotic='nonquantified open neighborhood in weak b, phase, positive real strong amplitude ratio',
                          global_near_d2='nonquantified open neighborhood by strict Stage-2 inequalities and transversality',
                          numerical_bounds=None,
                          connectedness_from_d0_to_d2='not_certified'),
         point_comparisons=rows,source_sha256=hashes,
         theorem_verified=bool(iv['verified'] and stage2['verified']),
         quantitative_region_certified=False)
(ROOT/'theorem_region_certificate.json').write_text(json.dumps(out,indent=2))
print('theorem_verified',out['theorem_verified'],'quantitative_region_certified',out['quantitative_region_certified'])
