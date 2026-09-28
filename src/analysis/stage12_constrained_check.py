import stage12_branch_continuation as s
import numpy as np,json,csv
from scipy.optimize import minimize_scalar
q=json.load(open('stage12_continuation_checks.json'));ec=float(q['epsilon_cross']);grid=np.linspace(-21*np.pi,21*np.pi,4097);vals=[]
for x in grid:
 try:vals.append(s.calc([x,10],ec)[0])
 except np.linalg.LinAlgError:vals.append(float('inf'))
mins=[]
for k in range(1,len(grid)-1):
 if vals[k]<vals[k-1] and vals[k]<vals[k+1] and not grid[k-1]<10<grid[k+1]:
  a=minimize_scalar(lambda v:s.calc([v,10],ec)[0],bounds=(grid[k-1],grid[k+1]),method='bounded',options={'xatol':1e-12});mins.append([float(a.fun),float(a.x)])
q['satellite_fixed_b_full_period_numerical_minima']=sorted(mins);rs=list(csv.DictReader(open('branch_continuation.csv')))
q['minimum_sampled_curvature']={b:min(float(r['Hessian_min']) for r in rs if r['branch']==b) for b in ('A','B')}
json.dump(q,open('stage12_continuation_checks.json','w'),indent=2)
print(q['minimum_sampled_curvature']);print(sorted(mins)[:3])
