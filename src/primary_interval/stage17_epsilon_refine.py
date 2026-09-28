"""Targeted epsilon-quadratic refinement of archived outer terminal cells.

This checks only outer exclusion. Inner neighborhoods and parametric roots are
separate proof obligations and are never inferred from this output.
"""
import json
import sys
from decimal import Decimal

epsrad=sys.argv[1] if len(sys.argv)>1 else '0.0025'
limit=int(sys.argv[2]) if len(sys.argv)>2 else 40
sys.argv=['stage17_epsilon_quadratic_trial.py',epsrad,'0']
import stage17_epsilon_quadratic_trial as q

leaves=[];unresolved=[];visited=0
def predicate(row):
    coef,theta=q.cell_coeff(row)
    if coef is None:return False,'denominator',None
    Jmax=q.qmax(coef,q.lo,q.hi)
    Jmin=q.qmin(coef,q.lo,q.hi)
    Emax=q.qmax(q.E,q.lo,q.hi)
    Emin=q.qmin(q.E,q.lo,q.hi)
    if Jmin.a<=(theta.b**2*Emax).b:
        return False,'nonpositive_residual_bound',None
    penalty=2*theta.b*q.iv.sqrt(Jmax*Emax)-theta.a**2*Emin
    for branch,U in q.incumbents.items():
        diff=q.qmin(tuple(coef[k]-U[k] for k in range(3)),q.lo,q.hi)
        if diff.a>penalty.b:return True,branch,str((diff-penalty).a)
    return False,'quadratic_gap',None

def refine(row,depth,parent):
    global visited
    visited+=1
    ok,why,margin=predicate(row)
    if ok:
        leaves.append(dict(parent_index=parent,cell={k:row[k] for k in ('m_lo','m_hi','h_lo','h_hi')},depth=depth,branch_incumbent=why,margin_lower=margin))
        return
    if depth>=20 or visited>100000:
        unresolved.append(dict(parent_index=parent,cell={k:row[k] for k in ('m_lo','m_hi','h_lo','h_hi')},depth=depth,reason=why))
        return
    a,b,c,d=(Decimal(row[k]) for k in ('m_lo','m_hi','h_lo','h_hi'))
    if b-a>=d-c:
        mid=(a+b)/2
        r1=dict(row,m_hi=str(mid));r2=dict(row,m_lo=str(mid))
    else:
        mid=(c+d)/2
        r1=dict(row,h_hi=str(mid));r2=dict(row,h_lo=str(mid))
    refine(r1,depth+1,parent)
    refine(r2,depth+1,parent)

for i,row in enumerate(q.rows[:limit]):
    refine(row,0,i)
    if (i+1)%10==0:print('parents',i+1,'visited',visited,'unresolved',len(unresolved),flush=True)
out=dict(epsilon_interval=[str(q.lo),str(q.hi)],parents=limit,outer_total=len(q.rows),
         visited=visited,leaves=len(leaves),unresolved=unresolved,
         max_depth=max([z['depth'] for z in leaves] or [0]),
         complete_outer=limit==len(q.rows) and not unresolved,
         proof_scope='outer only',refined_leaves=leaves)
(q.ROOT/f'stage17_epsilon_outer_refinement_{epsrad.replace(".","p")}_{limit}.json').write_text(json.dumps(out,indent=2))
print(json.dumps({k:v for k,v in out.items() if k!='refined_leaves'},indent=2))
