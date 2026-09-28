"""Assemble Stage-18 evidence ledger without upgrading numerical rows to proofs."""
import csv, hashlib, json
from decimal import Decimal
from pathlib import Path

R=Path(__file__).resolve().parent
def read(name):
    with (R/name).open(newline='',encoding='utf-8') as f:return list(csv.DictReader(f))
def write(name,rows):
    with (R/name).open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def h(name):return hashlib.sha256((R/name).read_bytes()).hexdigest()

numerical=read('p_bridge_numerical.csv')
pilot=read('stage18_local_slab_pilot.csv')
base=json.loads((R/'stage4_certificate.json').read_text())
slabs=[]
for q in pilot:
    slabs.append(dict(p_lo=q['p_lo'],p_hi=q['p_hi'],z_lo=str(Decimal(q['z_mid'])-Decimal(q['z_width'])/2),
                      z_hi=str(Decimal(q['z_mid'])+Decimal(q['z_width'])/2),
                      status='local_crossing_interval_validated' if q['verified']=='True' else 'local_crossing_interval_failed',
                      local_validated=q['verified'],global_validated=False,
                      contraction_upper=q['contraction_upper'],A_Hdet_lower=q['A_Hdet_lower'],
                      B_Hdet_lower=q['B_Hdet_lower'],slope_lower=q['slope_lower'],
                      competitor_margin_lower='',source='stage18_local_slab_pilot.csv'))
for q in json.loads((R/'stage4_interval_boxes.json').read_text()):
    slabs.append(dict(p_lo=str(Decimal(q['z_lo'])**2),p_hi=str(Decimal(q['z_hi'])**2),
                      z_lo=q['z_lo'],z_hi=q['z_hi'],status='archived_local_crossing_interval_validated',
                      local_validated=q['verified'],global_validated=False,
                      contraction_upper=q['contraction_upper'],A_Hdet_lower=q['A_Hdet_lower'],
                      B_Hdet_lower=q['B_Hdet_lower'],slope_lower=q['slope_lower'],
                      competitor_margin_lower='',source='stage4_interval_boxes.json'))
g=base['connected_global_interval']
slabs.append(dict(p_lo=str(Decimal(g['z_interval'][0])**2),p_hi='4',
                  z_lo=g['z_interval'][0],z_hi=g['z_interval'][1],
                  status='archived_full_domain_global_interval_validated',
                  local_validated=True,global_validated=True,contraction_upper='',
                  A_Hdet_lower='',B_Hdet_lower='',slope_lower='',
                  competitor_margin_lower=g['outer_min_margin'],source='stage4_global_interval_result.json'))
write('p_bridge_slabs.csv',slabs)

margin=[]
for name in ('globality_margins.csv','stage4_bridge_third_branch_probe.csv'):
    for q in read(name):
        z=q.get('z_NDelta') or q.get('z')
        if not z:continue
        margin.append(dict(p=str(Decimal(z)**2),z=z,
                           third_candidate_cost=q.get('third_candidate_cost',''),
                           third_candidate_margin=q.get('third_candidate_margin') or q.get('third_margin',''),
                           confluent_margin=q.get('confluent_margin',''),
                           status='numerical_search_only',source=name))
margin.sort(key=lambda q:Decimal(q['p']))
margin.append(dict(p='[3.99999880000009,4]',z='[1.9999997,2]',third_candidate_cost='',
                   third_candidate_margin='',confluent_margin='',status='full_domain_exclusion_margin_not_third_best',
                   source='stage4_global_interval_result.json'))
write('p_bridge_competitor_margins.csv',margin)

fields=('A_hessian_min','B_hessian_min','A_gram_determinant','B_gram_determinant',
        'A_within_fit_separation','B_within_fit_separation','AB_frequency_pair_distance',
        'slope_dJ_depsilon')
mins={k:str(min(Decimal(q[k]) for q in numerical)) for k in fields}
best={}
for z in ('0.1','0.5','1.0','1.5','1.9'):
    yes=[q for q in pilot if q['z_mid']==z and q['verified']=='True']
    best[z]=max(yes,key=lambda q:Decimal(q['z_width'])) if yes else None
cert=dict(model='N=21,b=10,phi=pi; p=z^2',verdict='TSP-READY-WITH-DISCONNECTED-RESULTS',
          bridge_from_tangent_to_p4_certified=False,
          small_p_global_interval='exists 0<p<p0 with p0>0 and a tangent limit at p=0, but p0 is not explicit; no certified overlap with a finite-p slab',
          largest_explicit_connected_global_interval_p=['3.99999880000009','4'],
          other_certified_isolated_p=[1,2.25],
          archived_global_interval_z=g['z_interval'],
          local_crossing_interval_z=['1.9999','2'],
          numerical_continuation=dict(rows=len(numerical),p_sampled=[numerical[0]['p'],numerical[-1]['p']],
                                      no_fallbacks=all(q['tangent_fallback']=='False' for q in numerical),
                                      sampled_minima=mins,status='numerical_only'),
          bounded_subdivision_repair=dict(local_certified_pilot_widths_z={z:(q['z_width'] if q else None) for z,q in best.items()},
                                          full_domain_global_exclusion_attempted_on_new_slabs=False,
                                          conclusion='raw ordinary-coordinate interval dependency requires impractically narrow slabs near p=0; no connected global cover obtained'),
          obstruction='interval dependency and computational cost; no sampled structural event, but absence of such an event is not certified over the gap',
          overlap_maps=dict(chart_S_to_F='u_A1=v_A; u_A2=p*kappa_A; u_B1=p*kappa_B; u_B2=v_B; epsilon=p*lambda, p>0',
                            inverse='v_A=u_A1; kappa_A=u_A2/p; kappa_B=u_B1/p; v_B=u_B2; lambda=epsilon/p, p>0',
                            consistency='algebraic identity on p>0; no Stage-18 interval chart-overlap certificate near p=0'),
          source_sha256={name:h(name) for name in ('p_bridge_numerical.csv','stage18_local_slab_pilot.csv',
                      'p_bridge_slabs.csv','p_bridge_competitor_margins.csv','stage4_certificate.json',
                      'stage4_global_interval_result.json','stage4_interval_boxes.json')})
(R/'p_bridge_certificate.json').write_text(json.dumps(cert,indent=2),encoding='utf-8')
print(json.dumps(dict(verdict=cert['verdict'],slab_rows=len(slabs),margin_rows=len(margin),pilot_best={z:(q['z_width'] if q else None) for z,q in best.items()}),indent=2))
