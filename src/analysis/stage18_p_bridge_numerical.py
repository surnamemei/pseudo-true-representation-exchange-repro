"""High-precision p-continuation; no sampled row is a global certificate."""
import csv
import json
from pathlib import Path

from mpmath import mp
from three_to_two_tone_stage3_full_crossings import crossing, D
from three_to_two_tone_stage4_continuation import statistics

ROOT=Path(__file__).resolve().parent
mp.dps=90

probe={}
for name in ('globality_margins.csv','stage4_bridge_third_branch_probe.csv'):
    with (ROOT/name).open(newline='') as f:
        for row in csv.DictReader(f):
            if row.get('z_NDelta') or row.get('z'):
                probe[str(mp.mpf(row.get('z_NDelta') or row.get('z')))]=row

z_grid=[mp.mpf(x) for x in ('0.001','0.003','0.01','0.02','0.05')]
z_grid += [mp.mpf(k)/20 for k in range(2,41)]
z_grid=sorted(set(z_grid))
rows=[];history=[]
for z in z_grid:
    p=z*z
    seed=None
    if len(history)>=2:
        p0,x0=history[-2];p1,x1=history[-1]
        seed=tuple(x1[k]+(x1[k]-x0[k])*(p-p1)/(p1-p0) for k in range(5))
    fallback=False
    try:q=crossing(z,seed=seed)
    except (ValueError,ZeroDivisionError):
        q=crossing(z);fallback=True
    x=tuple(mp.mpf(q[k]) for k in ('u_A1','u_A2','u_B1','u_B2','epsilon_cross'))
    history.append((p,x))
    e=x[4]
    A=statistics(x[0],x[1],e,z)
    B=statistics(x[2],x[3],e,z)
    gA=21**2-D(x[1]-x[0],21)**2
    gB=21**2-D(x[3]-x[2],21)**2
    key=str(z)
    third=probe.get(key,{})
    raw_third=third.get('third_margin') or third.get('third_candidate_margin') or ''
    def s(v,d=35):return mp.nstr(v,d)
    rows.append(dict(
        p=s(p,20),z=s(z,20),epsilon_cross=s(e,65),lambda_cross=s(e/p,50),
        u_A1=s(x[0]),u_A2=s(x[1]),u_B1=s(x[2]),u_B2=s(x[3]),
        chartS_v_A=s(x[0]),chartS_kappa_A=s(x[1]/p),
        chartS_kappa_B=s(x[2]/p),chartS_v_B=s(x[3]),
        A_amp1=s(A['amp1']),A_amp2=s(A['amp2']),
        B_amp1=s(B['amp1']),B_amp2=s(B['amp2']),
        A_J=s(A['J']),B_J=s(B['J']),
        A_hessian_min=s(A['hmin']),B_hessian_min=s(B['hmin']),
        A_hessian_max=s(A['hmax']),B_hessian_max=s(B['hmax']),
        A_gram_determinant=s(gA),B_gram_determinant=s(gB),
        A_within_fit_separation=s(A['sep_u']),B_within_fit_separation=s(B['sep_u']),
        AB_frequency_pair_distance=s(mp.sqrt(sum((x[k]-x[k+2])**2 for k in range(2)))),
        slope_dJ_depsilon=s(A['dJ_de']-B['dJ_de']),
        third_candidate_margin=str(raw_third),
        third_status=('archived_numerical_probe' if raw_third else 'not_searched'),
        secant_predictor_used=seed is not None,tangent_fallback=fallback,
        status='high_precision_numerical_only'))
    print('p',s(p,9),'lambda',s(e/p,14),'A_Hmin',s(A['hmin'],7),
          'B_Hmin',s(B['hmin'],7),flush=True)

out=ROOT/'p_bridge_numerical.csv'
with out.open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print(json.dumps(dict(rows=len(rows),p_min=rows[0]['p'],p_max=rows[-1]['p'],
                      fallbacks=sum(r['tangent_fallback'] for r in rows),
                      third_probe_rows=sum(r['third_status']=='archived_numerical_probe' for r in rows)),indent=2))
