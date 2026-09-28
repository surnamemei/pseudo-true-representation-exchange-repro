"""Adaptive epsilon/spatial replay of the archived N=21 inner cover.

Requires a separate Krawczyk proof for each epsilon slab. Every accepted
terminal is either excluded by cost/gradient/Krawczyk or contained in that
slab's certified A/B root box.
"""
import csv
import json
import sys
from decimal import Decimal
from pathlib import Path

import numpy as np
from mpmath import iv,mp

import stage17_parametric_jet as p
root_suffix=sys.argv[1] if len(sys.argv)>1 else ''
visit_budget=int(sys.argv[2]) if len(sys.argv)>2 else 1600000
sys.argv=['stage17_epsilon_quadratic_trial.py','0.0025','0']
import stage17_epsilon_quadratic_trial as q

ROOT=Path(__file__).resolve().parent
iv.dps=80;mp.dps=85
q.it.N=21
roots=json.loads((ROOT/f'stage17_epsilon_root_tiling{root_suffix}.json').read_text())
assert roots['all_pass']
bins=roots['records']
num_bins=roots['slabs']
archived=[r for r in csv.DictReader((ROOT/'stage2_inner_cells.csv').open())
          if r['label']=='crossing']
cache={}
leaves=[];unresolved=[];visited=0;max_u_depth=0

def box_data(row):
    key=tuple(row[k] for k in ('u1_lo','u1_hi','u2_lo','u2_hi'))
    if key in cache:return cache[key]
    a,b,c,d=(Decimal(v) for v in key)
    X=[iv.mpf([str(a),str(b)]),iv.mpf([str(c),str(d)])]
    mid=[iv.mpf(str((a+b)/2)),iv.mpf(str((c+d)/2))]
    Y=[X[k]-mid[k] for k in range(2)]
    j0=p.coefficients(*mid)
    jX=p.coefficients(*X)
    lower=[];gcoeff=[[],[]]
    for k in range(3):
        v=j0[k].v
        for i in range(2):
            v+=j0[k].g[i]*Y[i]
            for j in range(2):
                v+=iv.mpf('0.5')*jX[k].H[i][j]*Y[i]*Y[j]
        lower.append(v)
        for i in range(2):
            gcoeff[i].append(j0[k].g[i]+sum(jX[k].H[i][j]*Y[j] for j in range(2)))
    out=dict(X=X,mid=mid,Y=Y,j0=j0,jX=jX,lower=lower,gcoeff=gcoeff)
    cache[key]=out
    return out

def inside_root(row,branch,k):
    X=bins[k]['branches'][branch]['X']
    def endpoint(s):return Decimal(s.strip('[]').split(',')[0].strip())
    slack=Decimal('1e-20')
    return all(Decimal(row[f'u{i}_lo'])>=endpoint(X[i-1][0])+slack and
               Decimal(row[f'u{i}_hi'])<=endpoint(X[i-1][1])-slack
               for i in (1,2))

def ebox(i,j):
    return iv.mpf([bins[i]['epsilon_interval'][0],bins[j-1]['epsilon_interval'][1]])

def classify(row,branch,i,j):
    if j==i+1 and inside_root(row,branch,i):return 'inside_Krawczyk_box',''
    data=box_data(row)
    E=ebox(i,j)
    elo=bins[i]['epsilon_interval'][0]
    ehi=bins[j-1]['epsilon_interval'][1]
    for w,U in q.incumbents.items():
        diff=q.qmin(tuple(data['lower'][k]-U[k] for k in range(3)),elo,ehi)
        if diff.a>0:return 'objective_excluded',str(diff.a)
    for k in range(2):
        g=q.val(data['gcoeff'][k],E)
        if g.a>0 or g.b<0:
            return 'gradient_excluded',str(g.a if g.a>0 else -g.b)
    J0=q.val(data['j0'],E);JX=q.val(data['jX'],E)
    try:
        hf=np.array([[float(J0.H[k][j].mid) for j in range(2)] for k in range(2)])
        inv=np.linalg.inv(hf)
        C=[[iv.mpf(repr(float(inv[k,j]))) for j in range(2)] for k in range(2)]
        M=[[iv.mpf(int(k==j))-sum(C[k][r]*JX.H[r][j] for r in range(2)) for j in range(2)] for k in range(2)]
        K=[data['mid'][k]-sum(C[k][r]*J0.g[r] for r in range(2))
           +sum(M[k][r]*data['Y'][r] for r in range(2)) for k in range(2)]
        if any(K[k].b<data['X'][k].a or K[k].a>data['X'][k].b for k in range(2)):
            return 'krawczyk_excluded',''
    except (np.linalg.LinAlgError,ZeroDivisionError,ValueError):
        pass
    return None,''

def split_u(row):
    a,b,c,d=(Decimal(row[k]) for k in ('u1_lo','u1_hi','u2_lo','u2_hi'))
    if b-a>=d-c:
        mid=(a+b)/2
        return dict(row,u1_hi=str(mid)),dict(row,u1_lo=str(mid))
    mid=(c+d)/2
    return dict(row,u2_hi=str(mid)),dict(row,u2_lo=str(mid))

def replay(row,branch,parent,i,j,depth_u=0):
    global visited,max_u_depth
    visited+=1
    pred,margin=classify(row,branch,i,j)
    if pred:
        leaves.append(dict(parent_index=parent,branch=branch,epsilon_slabs=[i,j],
                           u_depth=depth_u,predicate=pred,margin_lower=margin,
                           cell={k:row[k] for k in ('u1_lo','u1_hi','u2_lo','u2_hi')}))
        return
    if visited>=visit_budget:
        unresolved.append(dict(parent_index=parent,branch=branch,epsilon_slabs=[i,j],
                               u_depth=depth_u,reason='budget'))
        return
    if j-i>1:
        m=(i+j)//2
        replay(row,branch,parent,i,m,depth_u)
        replay(row,branch,parent,m,j,depth_u)
        return
    if depth_u>=20:
        unresolved.append(dict(parent_index=parent,branch=branch,epsilon_slabs=[i,j],
                               u_depth=depth_u,reason='u_depth'))
        return
    max_u_depth=max(max_u_depth,depth_u+1)
    a,b=split_u(row)
    replay(a,branch,parent,i,j,depth_u+1)
    replay(b,branch,parent,i,j,depth_u+1)

for k,row in enumerate(archived):
    replay(row,row['branch'],k,0,num_bins)
    if (k+1)%50==0:print('parents',k+1,'visited',visited,'leaves',len(leaves),'unresolved',len(unresolved),flush=True)
out=dict(epsilon_interval=roots['interval'],archived_terminal_rows=len(archived),
         visited=visited,leaves=len(leaves),unresolved=unresolved,
         max_u_depth=max_u_depth,unique_frequency_boxes=len(cache),
         complete=not unresolved,
         predicate_counts={kind:sum(z['predicate']==kind for z in leaves) for kind in
                           ('inside_Krawczyk_box','objective_excluded','gradient_excluded','krawczyk_excluded')},
         records=leaves)
(ROOT/f'stage17_epsilon_inner_tiled{root_suffix}.json').write_text(json.dumps(out,indent=2))
print(json.dumps({k:v for k,v in out.items() if k!='records'},indent=2))
