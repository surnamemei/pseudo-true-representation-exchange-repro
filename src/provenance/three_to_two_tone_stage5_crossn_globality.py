"""Independent numerical full-torus checks of selected cross-N equal-cost roots."""
import csv
from pathlib import Path
import numpy as np
from three_to_two_tone_stage5_core import signal,times,grid_minima,refine

ROOT=Path(__file__).resolve().parent
data=list(csv.DictReader((ROOT/'crossN_scaling.csv').open()))
out=[]
for r in data:
    if r['status']!='numerical_stationary_crossing' or float(r['z']) not in (.6,1.,2.):continue
    n=int(r['N']);z=float(r['z']);e=float(r['epsilon_cross']);s=times(n)
    x=signal(n,z,10,e)
    a=refine([float(r['u_A1']),float(r['u_A2'])],x,s)
    b=refine([float(r['u_B1']),float(r['u_B2'])],x,s)
    roots=grid_minima(x,s,grid_size=256,top=32)
    closest_A=min((np.linalg.norm(q['u']-a['u']) for q in roots),default=np.inf)
    closest_B=min((np.linalg.norm(q['u']-b['u']) for q in roots),default=np.inf)
    other=[q['J'] for q in roots if np.linalg.norm(q['u']-a['u'])>.05 and
                                     np.linalg.norm(q['u']-b['u'])>.05]
    top2=all(any(np.linalg.norm(q['u']-v['u'])<.05 for q in roots[:2]) for v in (a,b))
    out.append(dict(N=n,z=z,epsilon_cross=e,grid_size=256,n_distinct_roots=len(roots),
                    labelled_A_found=closest_A<.05,labelled_B_found=closest_B<.05,
                    A_B_are_top_two=top2,third_margin=min(other)-a['J'] if other else '',
                    A_hessian_min=min(a['hessian_eigs']),B_hessian_min=min(b['hessian_eigs']),
                    A_B_distance=np.linalg.norm(a['u']-b['u']),
                    status='numerical_only'))
    print(n,z,top2,out[-1]['third_margin'],flush=True)
with (ROOT/'crossN_globality_checks.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(out[0]));w.writeheader();w.writerows(out)
