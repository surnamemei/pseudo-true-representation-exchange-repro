"""Directed-interval verification of full-torus projector-distance cells.

The float64 adaptive search proposes a complete partition. Every excluded
cell is rechecked by mpmath.iv; accepted cells must lie in the two boxes
passed to the independent inner interval search.
"""

import csv
import json
from decimal import Decimal,getcontext
from pathlib import Path

import numpy as np
from mpmath import iv

import three_to_two_tone_branch_study as st
from three_to_two_tone_stage1_global import adaptive
from three_to_two_tone_stage2_interval import data_energy

ROOT=Path(__file__).resolve().parent
iv.dps=65;getcontext().prec=100
ref=json.loads((ROOT/'stage2_reference.json').read_text(encoding='utf-8'))


def sinc(z):return iv.sin(z)/z
def sinc_prime(z):return (z*iv.cos(z)-iv.sin(z))/(z*z)


def center_projected_J(m,h,e,E,strong_d=2):
    h1=iv.mpf(2)-e;h2=iv.mpf(0)
    c1sq=iv.mpf(1);c2sq=iv.mpf(0)
    w1psq=iv.mpf(0);w2psq=iv.mpf(0)
    for t in range(1,11):
        tt=iv.mpf(t);z=h*tt;co=iv.cos(z);si=iv.sin(z)
        sc=sinc(z);sp=sinc_prime(z)
        p=iv.cos((iv.mpf(strong_d)/21-m)*tt)+iv.cos((-iv.mpf(strong_d)/21-m)*tt)
        p-=e*iv.cos((iv.mpf(10)/21-m)*tt)
        q=iv.sin((iv.mpf(strong_d)/21-m)*tt)+iv.sin((-iv.mpf(strong_d)/21-m)*tt)
        q-=e*iv.sin((iv.mpf(10)/21-m)*tt)
        h1+=2*co*p
        h2+=2*tt*sc*q
        c1sq+=2*co*co
        c2sq+=2*tt*tt*sc*sc
        w1psq+=2*tt*tt*si*si
        w2psq+=2*tt**4*sp*sp
    J=E-h1*h1/c1sq-h2*h2/c2sq
    return J,iv.sqrt(c1sq),iv.sqrt(c2sq),iv.sqrt(w1psq),iv.sqrt(w2psq)


def verify_cell(z,e,E,incumbent_upper,inflate=Decimal('1e-12'),strong_d=2):
    dep,mlo,mhi,hlo,hhi,j_float,lb_float,status=z
    a,b,c,d=map(lambda v:Decimal(repr(float(v))),(mlo,mhi,hlo,hhi))
    mc=(a+b)/2;hc=(c+d)/2
    rm=(b-a)/2+inflate;rh=(d-c)/2+inflate
    m=iv.mpf(str(mc));h=iv.mpf(str(hc))
    rmi=iv.mpf(str(rm));rhi=iv.mpf(str(rh))
    Jc,c1,c2,w1p,w2p=center_projected_J(m,h,e,E,strong_d=strong_d)
    t2=iv.sqrt(iv.mpf(sum(t**4 for t in range(-10,11))))
    t3=iv.sqrt(iv.mpf(sum(t**6 for t in range(-10,11))))
    D1=w1p*rhi+iv.mpf('0.5')*t2*rhi*rhi
    D2=w2p*rhi+t3*rhi*rhi/6
    if not (c1.a>D1.b and c2.a>D2.b):
        return False,None,'denominator_nonpositive'
    theta=iv.sqrt(2)*10*rmi+iv.sqrt((D1/(c1-D1))**2+(D2/(c2-D2))**2)
    if Jc.a<=0:return False,None,'center_J_nonpositive'
    residual=iv.sqrt(Jc).a-iv.sqrt(E).b*theta.b
    bound=iv.mpf(0) if residual<=0 else residual*residual
    passed=bound.a>incumbent_upper
    return bool(passed),bound,'ok' if passed else 'bound_not_above_incumbent'


def run_one(rec):
    label=rec['label'];e0=rec['epsilon']
    e=iv.mpf(e0)+iv.mpf(['-1e-95','1e-95'])
    E=data_energy(e)
    # The reference objective is a rigorous candidate upper after an
    # independent Krawczyk enclosure. Use a conservative explicit decimal
    # threshold for the float proposal, then verify against that upper.
    from three_to_two_tone_stage2_krawczyk import check
    _,ja=check(label,e0,'A',return_taylor=True)
    _,jb=check(label,e0,'B',return_taylor=True)
    incumbent_upper=min(ja.b,jb.b)
    iv.dps=65
    (ua,fa),(ub,fb),x=st.two_branches(21,2,10,float(e0),np.pi)
    threshold=float(incumbent_upper)+1e-6
    complete,summary,leaves,un=adaptive(x,(ua,ub),threshold,21,
                                        radius_u=2.,maxdepth=38,write_rows=True)
    rows=[];fail=[]
    for z in leaves:
        dep,mlo,mhi,hlo,hhi,jf,lbf,status=z
        if status=='excluded':
            passed,bb,why=verify_cell(z,e,E,incumbent_upper)
            if not passed:fail.append((z,why))
            bound_text='' if bb is None else str(bb.a)
        else:
            passed=True;bound_text='';why='inside_A_or_B_box'
        rows.append(dict(label=label,depth=dep,m_lo=repr(float(mlo)),m_hi=repr(float(mhi)),
                         h_lo=repr(float(hlo)),h_hi=repr(float(hhi)),status=status,
                         proposed_J_lower=repr(float(lbf)),interval_J_lower=bound_text,
                         verified=passed,reason=why))
    return dict(label=label,proposal_complete=complete,terminal_cells=len(leaves),
                excluded_cells=sum(z[-1]=='excluded' for z in leaves),
                accepted_cells=sum(z[-1]!='excluded' for z in leaves),
                failed_interval_cells=len(fail),max_depth=summary[-1][0],
                incumbent_upper=str(incumbent_upper),
                accepted_box_radius_u='2.00000001',
                full_torus_interval_exclusion=complete and len(fail)==0),rows,fail


if __name__=='__main__':
    summaries=[];allrows=[];allfail=[]
    for rec in ref['records']:
        summary,rows,fail=run_one(rec)
        summaries.append(summary);allrows.extend(rows)
        allfail.extend((rec['label'],z,why) for z,why in fail)
        print(summary,flush=True)
    (ROOT/'stage2_outer_global_summary.json').write_text(json.dumps(summaries,indent=2),encoding='utf-8')
    with (ROOT/'interval_boxes.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(allrows[0]));w.writeheader();w.writerows(allrows)
    (ROOT/'stage2_outer_failures.json').write_text(json.dumps(allfail,indent=2,default=str),encoding='utf-8')
