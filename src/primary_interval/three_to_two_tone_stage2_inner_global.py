"""Interval branch-and-bound inside the two radius-2 frequency boxes.

This tries to exclude all other stationary points in those boxes. If the
search exhausts its budget, it reports unresolved cells instead of a proof.
"""

import csv
import json
from decimal import Decimal,getcontext
from pathlib import Path

import numpy as np
from mpmath import iv

import three_to_two_tone_branch_study as st
from three_to_two_tone_stage2_interval import objective_jet

ROOT=Path(__file__).resolve().parent
iv.dps=45;getcontext().prec=110
ref=json.loads((ROOT/'stage2_reference.json').read_text(encoding='utf-8'))
local=json.loads((ROOT/'stage2_local_interval_checks.json').read_text(encoding='utf-8'))


def dec_point_interval(lo,hi):
    return iv.mpf([str(lo),str(hi)])


def in_tiny(cell,ref_center,radius):
    (a,b,c,d)=cell
    u1,u2=map(Decimal,ref_center)
    return a>=u1-radius and b<=u1+radius and c>=u2-radius and d<=u2+radius


def process(label,which,e_str,float_center,best_upper,max_cells=100000,
            known_center_override=None,local_radius_override=None):
    e=iv.mpf(e_str)+iv.mpf(['-1e-95','1e-95'])
    f1,f2=[Decimal(repr(float(z))) for z in float_center]
    rad=Decimal('2.00000001')
    root=(f1-rad,f1+rad,f2-rad,f2+rad,0)
    known_center=(ref['A_at_crossing'] if which=='A' else ref['B_at_crossing'])
    if label in ('below','above'):
        known_center=next(r[which] for r in ref['records'] if r['label']==label)
    if known_center_override is not None:
        known_center=known_center_override
    # Stay strictly inside the larger interval-certified Krawczyk boxes.
    local_radius=Decimal('0.000299999999') if which=='A' else Decimal('0.009999999')
    if local_radius_override is not None:
        local_radius=Decimal(str(local_radius_override))
    stack=[root];log=[];nvisited=0;unresolved=[]
    while stack:
        cell=stack.pop();a,b,c,d,depth=cell;nvisited+=1
        if in_tiny((a,b,c,d),known_center,local_radius):
            log.append(dict(label=label,branch=which,depth=depth,u1_lo=str(a),u1_hi=str(b),
                            u2_lo=str(c),u2_hi=str(d),status='inside_Krawczyk_box',
                            J_lower='',gradient_exclusion=''))
            continue
        X=[dec_point_interval(a,b),dec_point_interval(c,d)]
        mid=[iv.mpf(str((a+b)/2)),iv.mpf(str((c+d)/2))]
        J0=objective_jet(mid[0],mid[1],e)
        JX=objective_jet(X[0],X[1],e)
        Y=[X[i]-mid[i] for i in range(2)]
        grad=[J0.g[i]+sum(JX.H[i][k]*Y[k] for k in range(2)) for i in range(2)]
        Jlower=J0.v
        for i in range(2):Jlower+=J0.g[i]*Y[i]
        for i in range(2):
            for k in range(2):Jlower+=iv.mpf('0.5')*JX.H[i][k]*Y[i]*Y[k]
        if Jlower.a>best_upper:
            status='objective_excluded';detail=''
        elif grad[0].a>0 or grad[0].b<0:
            status='gradient_excluded';detail='g1'
        elif grad[1].a>0 or grad[1].b<0:
            status='gradient_excluded';detail='g2'
        else:
            status=''
            try:
                hf=np.array([[float(J0.H[i][j].mid) for j in range(2)] for i in range(2)])
                cf=np.linalg.inv(hf)
                C=[[iv.mpf(repr(float(cf[i,j]))) for j in range(2)] for i in range(2)]
                M=[[iv.mpf(int(i==j))-sum(C[i][k]*JX.H[k][j] for k in range(2))
                    for j in range(2)] for i in range(2)]
                K=[mid[i]-sum(C[i][k]*J0.g[k] for k in range(2))
                   +sum(M[i][k]*Y[k] for k in range(2)) for i in range(2)]
                if any(K[i].b<X[i].a or K[i].a>X[i].b for i in range(2)):
                    status='krawczyk_excluded';detail='K disjoint'
            except (np.linalg.LinAlgError,ZeroDivisionError,ValueError):
                pass
        if not status:
            if nvisited>=max_cells:
                unresolved.append(cell)
                break
            if b-a>=d-c:
                mid=(a+b)/2
                stack.append((a,mid,c,d,depth+1));stack.append((mid,b,c,d,depth+1))
            else:
                mid=(c+d)/2
                stack.append((a,b,c,mid,depth+1));stack.append((a,b,mid,d,depth+1))
            continue
        log.append(dict(label=label,branch=which,depth=depth,u1_lo=str(a),u1_hi=str(b),
                        u2_lo=str(c),u2_hi=str(d),status=status,
                        J_lower=str(Jlower.a),gradient_exclusion=detail))
        if nvisited%2000==0:
            print(label,which,'visited',nvisited,'stack',len(stack),flush=True)
    complete=not unresolved and not stack
    return complete,nvisited,log,len(unresolved)+len(stack)


if __name__=='__main__':
    results=[];allrows=[]
    from mpmath import mp
    from three_to_two_tone_stage2_krawczyk import check
    mp.dps=90
    for rec in ref['records']:
        label=rec['label'];e=rec['epsilon']
        eps_float=float(e)
        (ua,ja),(ub,jb),_=st.two_branches(21,2,10,eps_float,3.141592653589793)
        _,JAcert=check(label,e,'A',return_taylor=True)
        _,JBcert=check(label,e,'B',return_taylor=True)
        best_upper=min(JAcert.b,JBcert.b)
        iv.dps=45
        for which,center in [('A',ua),('B',ub)]:
            done,nv,rows,remaining=process(label,which,e,center,best_upper,max_cells=10000)
            results.append(dict(label=label,branch=which,complete=done,visited=nv,
                                terminal=len(rows),remaining=remaining,
                                incumbent_upper=str(best_upper),
                                objective_excluded=sum(r['status']=='objective_excluded' for r in rows),
                                gradient_excluded=sum(r['status']=='gradient_excluded' for r in rows),
                                krawczyk_excluded=sum(r['status']=='krawczyk_excluded' for r in rows),
                                inside_tiny=sum(r['status']=='inside_Krawczyk_box' for r in rows)))
            allrows.extend(rows)
            print(results[-1],flush=True)
    (ROOT/'stage2_inner_global_summary.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
    with (ROOT/'stage2_inner_cells.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(allrows[0]));w.writeheader();w.writerows(allrows)
