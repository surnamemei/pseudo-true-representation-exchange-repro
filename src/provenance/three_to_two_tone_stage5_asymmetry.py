"""One-factor asymmetry sweep and compact two-parameter existence map.

Every non-baseline globality statement is numerical (grid plus local refinement).
"""
import csv
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from three_to_two_tone_stage5_core import times, signal, refine, grid_minima

ROOT = Path(__file__).resolve().parent
N, Z = 21, 2.0
S = times(N)
A0 = np.array([-4.772447353918614, .605170320276862])
B0 = np.array([-.057653561926680, 12.46341935680166])
EC0 = .2481906301722774


def evaluate(e, ratio=1, phase=np.pi, b=10, shift=0, seeds=None):
    y = signal(N, Z, b, e, phase, ratio, shift)
    aseed, bseed = seeds if seeds is not None else (A0, B0)
    a = refine(aseed, y, S)
    q = refine(bseed, y, S)
    return y, a, q


def one_case(family, value, ratio=1, phase=np.pi, b=10, shift=0,
             global_search=True):
    args = dict(ratio=ratio, phase=phase, b=b, shift=shift)
    probe = []
    for e in np.linspace(.07, .56, 15):
        _, a, q = evaluate(e, **args)
        # Retain only genuine distinct positive-curvature branches.
        valid = (np.linalg.norm(a['u']-q['u'])>.5 and
                 min(a['hessian_eigs'])>1e-6 and min(q['hessian_eigs'])>1e-6 and
                 a['converged'] and q['converged'])
        probe.append((e, a['J']-q['J'] if valid else np.nan))
    brackets = [(probe[i][0],probe[i+1][0]) for i in range(len(probe)-1)
                if np.isfinite(probe[i][1]) and np.isfinite(probe[i+1][1])
                and probe[i][1]*probe[i+1][1]<0]
    base = dict(family=family, value=value, ratio=ratio, phase=phase, b=b,
                strong_shift=shift, epsilon_c='', crossing_type='no_valid_bracket',
                globality='unresolved', third_margin='', u_A1='', u_A2='',
                u_B1='', u_B2='', fitted_separation_min='', branch_distance='',
                hessian_min_A='', hessian_min_B='', gram_cond_max='',
                deltaJ_slope='', objective_at_crossing='')
    if not brackets:
        return base
    # Select root closest to the baseline after perturbation.
    roots=[]
    for lo,hi in brackets:
        try:
            roots.append(brentq(lambda e: evaluate(e, **args)[1]['J']-
                                evaluate(e, **args)[2]['J'],lo,hi,xtol=1e-11))
        except ValueError:
            pass
    if not roots:
        return base
    e=min(roots,key=lambda q:abs(q-EC0))
    y,a,q=evaluate(e, **args)
    h=1e-4
    slope=((evaluate(e+h,**args)[1]['J']-evaluate(e+h,**args)[2]['J'])-
           (evaluate(e-h,**args)[1]['J']-evaluate(e-h,**args)[2]['J']))/(2*h)
    base.update(epsilon_c=e, crossing_type='noncoalescent_transverse' if
                abs(slope)>1e-3 and min(a['hessian_eigs'])>1e-6 and
                min(q['hessian_eigs'])>1e-6 else 'degenerate_or_uncertain',
                u_A1=a['u'][0],u_A2=a['u'][1],u_B1=q['u'][0],u_B2=q['u'][1],
                fitted_separation_min=min(np.diff(a['u'])[0],np.diff(q['u'])[0]),
                branch_distance=np.linalg.norm(a['u']-q['u']),
                hessian_min_A=min(a['hessian_eigs']),
                hessian_min_B=min(q['hessian_eigs']),
                gram_cond_max=max(a['gram_cond'],q['gram_cond']),
                deltaJ_slope=slope,objective_at_crossing=a['J'])
    if global_search:
        candidates=grid_minima(y,S,grid_size=192,top=24)
        top2=candidates[:2]
        matches=all(any(np.linalg.norm(fit['u']-p['u'])<.05 for p in top2)
                    for fit in (a,q)) if len(top2)==2 else False
        others=[p['J'] for p in candidates if all(np.linalg.norm(p['u']-f['u'])>.05
                 for f in (a,q))]
        base['third_margin']=min(others)-a['J'] if others else ''
        base['globality']='numerical_top_two_grid192' if matches else 'third_branch_or_grid_unresolved'
    return base


rows=[]
for ratio in (.8,.9,1,1.1,1.2):
    rows.append(one_case('strong_amplitude_ratio',ratio,ratio=ratio))
for phase in (0,np.pi/4,np.pi/2,3*np.pi/4,np.pi):
    rows.append(one_case('weak_phase',phase,phase=phase))
for b in (8,9,10,11,12):
    rows.append(one_case('weak_location',b,b=b))
for shift in (-.4,-.2,0,.2,.4):
    rows.append(one_case('strong_location_shift',shift,shift=shift))
with (ROOT/'asymmetry_sweep.csv').open('w',newline='') as file:
    w=csv.DictWriter(file,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
for row in rows:
    print(row['family'],row['value'],row['epsilon_c'],row['crossing_type'],
          row['globality'],row['third_margin'],flush=True)

mapping=[]
for ratio in (.8,.9,1,1.1,1.2):
    for b in (8,9,10,11,12):
        row=one_case('map_ratio_b',f'{ratio},{b}',ratio=ratio,b=b)
        row['classification']=('noncoalescent_branch_exchange' if
             row['crossing_type']=='noncoalescent_transverse' and
             row['globality']=='numerical_top_two_grid192' else
             'no_exchange_or_unresolved')
        mapping.append(row)
        print('map',ratio,b,row['classification'],row['epsilon_c'],flush=True)
with (ROOT/'existence_region.csv').open('w',newline='') as file:
    w=csv.DictWriter(file,fieldnames=list(mapping[0]));w.writeheader();w.writerows(mapping)
