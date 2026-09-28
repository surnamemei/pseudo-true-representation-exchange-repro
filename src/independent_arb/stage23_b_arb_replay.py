"""Independent Arb reevaluation of the proposed parametric local obligations.

This backend imports no mpmath evaluator. Failed local boxes cannot support a
global certificate; no global theorem is inferred from a failed inclusion.
"""
import sys,json,time
from pathlib import Path
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'stage12_deps'))
from flint import arb,ctx
from stage16_independent_arb import Dual,sinc_deriv
ctx.dps=85
def hull(x):return arb(x[0]).union(arb(x[1]))
def kd(q,k=0):
    ds=[sinc_deriv(q.x,k+j) for j in range(3)]
    return Dual(ds[0],tuple(ds[1]*g for g in q.g),
         tuple(tuple(ds[1]*q.h[i][j]+ds[2]*q.g[i]*q.g[j] for j in range(2)) for i in range(2)))
def tangent(v,l,b):
    v=Dual.variable(v,0);l=Dual.variable(l,1);b=Dual.const(b)
    m2=arb(1)/12;m4=arb(1)/80
    d,dp=kd(v),kd(v,1);q0=Dual.const(-m2)-l*kd(b);q1=2*l*kd(b,1)
    g=1-d*d-dp*dp/m2;h=kd(v,2)-l*kd(b-v)-d*q0+dp*q1/(2*m2)
    gain=h*h/g;cost=Dual.const(m4)-2*l*kd(b,2)+l*l-q0*q0-q1*q1/(4*m2)-gain
    return cost,g,h/g,gain
def system(x,b):
    a,ga,ba,pa=tangent(x[0],x[2],b);c,gc,bc,pc=tangent(x[1],x[2],b)
    return [a.g[0],c.g[0],pc.x-pa.x],[[a.h[0][0],arb(0),a.h[0][1]],[arb(0),c.h[0][0],c.h[0][1]],[a.g[0],-c.g[0],a.g[1]-c.g[1]]]
def run():
    start=time.perf_counter();primary=json.loads((ROOT/'b_interval_certificate.json').read_text());rows=[]
    for r in primary['attempts']:
        local=r['local'];x=[hull(q) for q in local['X']];x0=[hull(q) for q in local['center']];b=hull(r['b_interval'])
        ci=[[hull(q) for q in row] for row in local['preconditioner']]
        f,_=system(x0,b);_,j=system(x,b)
        mat=[[arb(int(i==k))-sum(ci[i][h]*j[h][k] for h in range(3)) for k in range(3)] for i in range(3)]
        kk=[x0[i]-sum(ci[i][k]*f[k] for k in range(3))+sum(mat[i][k]*(x[k]-x0[k]) for k in range(3)) for i in range(3)]
        good=all(v.is_finite() for v in kk)
        if good:good=all(kk[i].lower()>x[i].lower() and kk[i].upper()<x[i].upper() for i in range(3))
        rows.append(dict(half_width=r['half_width'],b_interval=r['b_interval'],inclusion=good,K=[str(v) for v in kk]))
        print(r['half_width'],'Arb inclusion',good,flush=True)
    out=dict(backend='python-flint 0.9.0/Arb',decimal_precision=ctx.dps,seconds=time.perf_counter()-start,
             attempts=rows,all_proposed_boxes_fail=not any(r['inclusion'] for r in rows),
             global_replay_performed=False,new_theorem=False)
    (ROOT/'b_interval_arb_replay.json').write_text(json.dumps(out,indent=2))
if __name__=='__main__':run()
