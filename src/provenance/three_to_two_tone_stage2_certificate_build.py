"""Assemble the interval certificate index from independently saved logs."""
import csv,json,re
from pathlib import Path
from mpmath import iv
import mpmath,numpy,scipy
from three_to_two_tone_stage2_interval import objective_jet

ROOT=Path(__file__).resolve().parent
NUM=re.compile(r'[-+]?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?')
def lo(s):return float(NUM.search(s).group())

ref=json.loads((ROOT/'stage2_reference.json').read_text())
local=json.loads((ROOT/'stage2_local_interval_checks.json').read_text())
large=json.loads((ROOT/'stage2_large_krawczyk_checks.json').read_text())
outer=json.loads((ROOT/'stage2_outer_global_summary.json').read_text())
inner=json.loads((ROOT/'stage2_inner_global_summary.json').read_text())
audit=json.loads((ROOT/'stage2_cell_audit.json').read_text())
cross=json.loads((ROOT/'stage2_crossing_interval.json').read_text())
crossglobal=json.loads((ROOT/'stage2_crossing_global.json').read_text())

iv.dps=90
hess=[]
for z in large:
    c=z['x0'];r=iv.mpf(z['radius_u'])
    X=[iv.mpf(c[i])+iv.mpf([-r.b,r.b]) for i in range(2)]
    jet=objective_jet(X[0],X[1],iv.mpf(z['epsilon']))
    H=jet.H
    det=H[0][0]*H[1][1]-H[0][1]*H[1][0]
    tr=H[0][0]+H[1][1]
    lower=(det.a/tr.b).a
    upper=tr.b
    hess.append(dict(label=z['label'],branch=z['branch'],box_radius_u=z['radius_u'],
                     lambda_min_lower=str(lower),lambda_max_upper=str(upper),
                     determinant_lower=str(det.a),trace_upper=str(tr.b),
                     verified=bool(det.a>0 and H[0][0].a>0 and lower>0)))

cases=[]
for rec in ref['records']:
    label=rec['label']
    a=next(z for z in local if z['label']==label and z['branch']=='A')
    b=next(z for z in local if z['label']==label and z['branch']=='B')
    oa=next(z for z in outer if z['label']==label)
    ca=next(z for z in audit if z['label']==label)
    ia=[z for z in inner if z['label']==label]
    cases.append(dict(label=label,epsilon=rec['epsilon'],A_box=a['x0'],B_box=b['x0'],
                      A_J_enclosure=a['J_root_enclosure'],B_J_enclosure=b['J_root_enclosure'],
                      A_Krawczyk=a['krawczyk_inclusion'],B_Krawczyk=b['krawczyk_inclusion'],
                      A_H_positive=a['positive_definite'],B_H_positive=b['positive_definite'],
                      outer=oa,inner=ia,cell_audit=ca,
                      global_result=bool(oa['full_torus_interval_exclusion'] and ca['complete']
                                         and all(z['complete'] for z in ia) and a['krawczyk_inclusion']
                                         and b['krawczyk_inclusion'] and a['positive_definite']
                                         and b['positive_definite'])))

out=dict(title='N=21 finite-record 3-tone to 2-tone globality log',
         model='x_t=2*cos(2*t/21)-epsilon*exp(i*10*t/21), t=-10,...,10; u_j=21*nu_j',
         arithmetic=dict(reference_dps=110,local_and_crossing_interval_dps=90,
                         outer_interval_dps=65,inner_interval_dps=45,
                         mpmath=mpmath.__version__,numpy=numpy.__version__,scipy=scipy.__version__,
                         directed_rounding='mpmath.iv via libmpi round_floor/round_ceiling',
                         independent_library_verification='not_performed'),
         domain=dict(frequency_torus='ordered two frequencies modulo 2*pi',
                     parameterization='m=(nu1+nu2)/2; h=(nu2-nu1)/2; m in [-pi,pi], h in [0,pi/2]',
                     outer_cell_inflation_rad='1e-12',
                     inner_box_radius_u='2.00000001',
                     coalescent_strip_h='[0,0.03]',
                     coalescent_u_separation='[0,1.26]'),
         cases=cases,hessian_eigenvalue_bounds=hess,
         crossing=cross,crossing_global_bracket=crossglobal,
         cell_log='interval_boxes.csv',inner_cell_log='stage2_inner_cells.csv',
         verified=bool(all(z['global_result'] for z in cases)
                       and cross['bracket_certified'] and crossglobal['full_bracket_global']
                       and all(z['verified'] for z in hess)))
(ROOT/'validated_globality_certificate.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print('verified',out['verified'],'cases',len(cases),'hessian',len(hess))
