"""Replay the complete N=21 crossing inner cover on a short epsilon interval."""
import csv
import json
import sys
from decimal import Decimal

import numpy as np
from mpmath import iv

radius=sys.argv[1] if len(sys.argv)>1 else '0.00001'
limit_arg=int(sys.argv[2]) if len(sys.argv)>2 else None
max_depth_arg=int(sys.argv[3]) if len(sys.argv)>3 else 16
sys.argv=['stage17_epsilon_quadratic_trial.py',radius,'0']
import stage17_epsilon_quadratic_trial as q
import stage17_parametric_jet as p

iv.dps=80
q.it.N=21
E=iv.mpf([str(q.lo),str(q.hi)])
center=q.ref
root_radii={'A':Decimal('0.0003'),'B':Decimal('0.01')}
rows=[r for r in csv.DictReader((q.ROOT/'stage2_inner_cells.csv').open())
      if r['label']=='crossing']
leaves=[];unresolved=[];visited=0

def inside_root(row,branch):
    if q.epsrad>q.mp.mpf('0.00001'):
        return False
    x=[Decimal(str(v)) for v in center[branch+'_at_crossing']]
    rad=root_radii[branch]
    return (Decimal(row['u1_lo'])>=x[0]-rad and Decimal(row['u1_hi'])<=x[0]+rad
            and Decimal(row['u2_lo'])>=x[1]-rad and Decimal(row['u2_hi'])<=x[1]+rad)

def classify(row,branch):
    if inside_root(row,branch):return 'inside_Krawczyk_box',''
    a,b,c,d=(Decimal(row[k]) for k in ('u1_lo','u1_hi','u2_lo','u2_hi'))
    X=[iv.mpf([str(a),str(b)]),iv.mpf([str(c),str(d)])]
    mid=[iv.mpf(str((a+b)/2)),iv.mpf(str((c+d)/2))]
    Y=[X[i]-mid[i] for i in range(2)]
    j0=p.coefficients(*mid)
    jX=p.coefficients(*X)
    lower=[];gcoeff=[[],[]]
    for j in range(3):
        jj0,jjX=j0[j],jX[j]
        v=jj0.v
        for i in range(2):
            v+=jj0.g[i]*Y[i]
            for k in range(2):
                v+=iv.mpf('0.5')*jjX.H[i][k]*Y[i]*Y[k]
        lower.append(v)
        for i in range(2):
            gcoeff[i].append(jj0.g[i]+sum(jjX.H[i][k]*Y[k] for k in range(2)))
    for which,U in q.incumbents.items():
        diff=q.qmin(tuple(lower[j]-U[j] for j in range(3)),q.lo,q.hi)
        if diff.a>0:return 'objective_excluded',str(diff.a)
    for i in range(2):
        g=q.val(gcoeff[i],E)
        if g.a>0 or g.b<0:
            return 'gradient_excluded',str(g.a if g.a>0 else -g.b)
    J0=q.val(j0,E);JX=q.val(jX,E)
    try:
        hf=np.array([[float(J0.H[i][j].mid) for j in range(2)] for i in range(2)])
        inv=np.linalg.inv(hf)
        C=[[iv.mpf(repr(float(inv[i,j]))) for j in range(2)] for i in range(2)]
        M=[[iv.mpf(int(i==j))-sum(C[i][k]*JX.H[k][j] for k in range(2)) for j in range(2)] for i in range(2)]
        K=[mid[i]-sum(C[i][k]*J0.g[k] for k in range(2))+sum(M[i][k]*Y[k] for k in range(2)) for i in range(2)]
        if any(K[i].b<X[i].a or K[i].a>X[i].b for i in range(2)):
            return 'krawczyk_excluded',''
    except (np.linalg.LinAlgError,ZeroDivisionError,ValueError):
        pass
    return None,''

def split(row):
    a,b,c,d=(Decimal(row[k]) for k in ('u1_lo','u1_hi','u2_lo','u2_hi'))
    if b-a>=d-c:
        mid=(a+b)/2
        return dict(row,u1_hi=str(mid)),dict(row,u1_lo=str(mid))
    mid=(c+d)/2
    return dict(row,u2_hi=str(mid)),dict(row,u2_lo=str(mid))

def replay(row,branch,parent,depth=0):
    global visited
    visited+=1
    pred,margin=classify(row,branch)
    if pred:
        leaves.append(dict(parent_index=parent,branch=branch,depth=depth,
                           predicate=pred,margin_lower=margin,
                           cell={k:row[k] for k in ('u1_lo','u1_hi','u2_lo','u2_hi')}))
        return
    if depth>=max_depth_arg or visited>100000:
        unresolved.append(dict(parent_index=parent,branch=branch,depth=depth,
                               cell={k:row[k] for k in ('u1_lo','u1_hi','u2_lo','u2_hi')}))
        return
    l,r=split(row)
    replay(l,branch,parent,depth+1)
    replay(r,branch,parent,depth+1)

limit=limit_arg if limit_arg is not None else len(rows)
for i,row in enumerate(rows[:limit]):
    replay(row,row['branch'],i)
    if (i+1)%100==0:print('parents',i+1,'visited',visited,'unresolved',len(unresolved),flush=True)
out=dict(epsilon_interval=[str(q.lo),str(q.hi)],archived_terminal_rows=len(rows),
         replayed=limit,visited=visited,terminal_leaves=len(leaves),unresolved=unresolved,
         complete=limit==len(rows) and not unresolved,
         predicate_counts={k:sum(z['predicate']==k for z in leaves) for k in
                           ('inside_Krawczyk_box','objective_excluded','gradient_excluded','krawczyk_excluded')},
         leaves=leaves)
(q.ROOT/f'stage17_epsilon_inner_{radius.replace(".","p")}.json').write_text(json.dumps(out,indent=2))
print(json.dumps({k:v for k,v in out.items() if k!='leaves'},indent=2))
