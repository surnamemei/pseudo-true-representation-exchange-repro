"""Five fixed-phase N=21 tangent diagnostics (numerical, not certified)."""
from __future__ import annotations
import csv
from pathlib import Path
import numpy as np
from mpmath import mp
from scipy.optimize import minimize_scalar
from stage20_phase_regimes import P,H,K_coal,N

mp.dps=60
R=lambda v,mu:P(v)+mu*H(v)
eq=lambda a,b,mu:(mp.diff(lambda v:R(v,mu),a),mp.diff(lambda v:R(v,mu),b),R(a,mu)-R(b,mu))
a,b,mu=mp.findroot(eq,(mp.mpf('1.6'),mp.mpf('10.3'),mp.mpf('.0984')),tol=mp.mpf('1e-45'))
star=R(a,mu)
grid=np.linspace(-np.pi*N,np.pi*N,8001)
values=np.array([float(R(mp.mpf(float(v)),mu)) if abs(v)>.03 else np.inf for v in grid])
roots=[]
for i in range(1,len(grid)-1):
    if abs(grid[i])>.1 and values[i]<values[i-1] and values[i]<values[i+1]:
        fit=minimize_scalar(lambda v:float(R(mp.mpf(v),mu)),bounds=(grid[i-1],grid[i+1]),method='bounded')
        roots.append((fit.x,fit.fun))
third=min(cost for v,cost in roots if abs(v-float(a))>.1 and abs(v-float(b))>.1)
rows=[]
for ratio in ('.05','.10','.20','.30','.40'):
    d=mp.sin(mp.pi*mp.mpf(ratio))
    fac=4*d*d
    rows.append(dict(delta_over_pi=ratio,sin_delta=mp.nstr(d,18),
        epsilon_over_z_critical=mp.nstr(2*d*mp.sqrt(mu),18),
        normalized_mu_c=mp.nstr(mu,24),v_A=mp.nstr(a,18),v_B=mp.nstr(b,18),
        A_B_cost=mp.nstr(fac*star,18),
        coalescent_gap=mp.nstr(fac*(mu*K_coal()-star),18),
        sampled_third_regular_gap=mp.nstr(fac*(mp.mpf(third)-star),18),
        Rvv_A=mp.nstr(fac*mp.diff(lambda v:R(v,mu),a,2),18),
        Rvv_B=mp.nstr(fac*mp.diff(lambda v:R(v,mu),b,2),18),
        globality_status='full-period numerical scan; not interval certified'))
path=Path(__file__).resolve().parent/'fixed_delta_phase_map.csv'
with path.open('w',newline='',encoding='utf-8') as f:
    writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
print(path)
