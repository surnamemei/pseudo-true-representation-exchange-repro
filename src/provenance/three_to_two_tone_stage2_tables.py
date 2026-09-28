"""Stage-2 scaling and regime tables with explicit validation status."""
import csv,json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent
ref=json.loads((ROOT/'stage2_reference.json').read_text())
lam=json.loads((ROOT/'stage2_lambda_largeN.json').read_text())

def save(path,rows):
    with (ROOT/path).open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

scaling=[]
for r in csv.DictReader((ROOT/'scaling_fit.csv').open(newline='')):
    n=int(r['N']);d=float(r['d_NDelta']);v=float(r['epsilon_cross']);pred=float(r['epsilon_tangent'])
    certified=(n==21 and d==2.0)
    value=ref['epsilon_cross'] if certified else r['epsilon_cross']
    scaling.append(dict(N=n,d_NDelta=d,Delta_per_sample=d/n,b_Nomega3=10.,
                        weak_phase_rad='pi',epsilon_cross=value,epsilon_float64=r['epsilon_cross'],
                        lambda_tangent=float(r['lambda_tangent']),epsilon_leading=pred,
                        relative_error=abs(float(value)-pred)/float(value),
                        evidence=('interval_certified_global_baseline' if certified else
                                  'float64_local_branch_crossing_not_global_certificate')))
save('finiteN_scaling.csv',scaling)

ns=np.array([v['N'] for v in lam[-4:]],float)
ys=np.array([v['lambda'] for v in lam[-4:]],float)
fit=np.polyfit(1/ns**2,ys,1)
coeff=[]
for r in lam:
    coeff.append(dict(N=r['N'],lambda_N=r['lambda'],lambda_infinity_estimate=fit[1],
                      N2_correction_fit=fit[0],geometry='strong u=+-d; weak u=10; phase=pi',
                      derivation='finite-sum real tangent-space projected residual',
                      coefficient_status='numerical root of explicit finite-sum equation'))
save('asymptotic_coefficients.csv',coeff)

reg=[]
for r in csv.DictReader((ROOT/'crossing_map.csv').open(newline='')):
    if r['strong_amp_ratio']!='1.0':continue
    isbaseline=(r['b_Nomega3']=='10.0' and r['phase_rad']=='3.141592653589793')
    status='certified_noncoalescent_exchange_at_baseline' if isbaseline else (
        'candidate_noncoalescent_exchange_unvalidated' if r['status']=='candidate' else
        'no_A_B_crossing_located_unclassified')
    reg.append(dict(N=int(r['N']),d_NDelta=float(r['d_NDelta']),
                    b_Nomega3=float(r['b_Nomega3']),r_location=float(r['r_location']),
                    phase_rad=float(r['phase_rad']),epsilon_cross=(ref['epsilon_cross'] if isbaseline else r['epsilon_cross']),
                    region_status=status,coalescent_region_boundary='not_determined',
                    globality='interval_certified_baseline_only' if isbaseline else 'not_certified'))
save('regime_boundaries.csv',reg)
print('scaling',len(scaling),'lambda',len(coeff),'regime',len(reg),'lambda_infty',fit[1])
