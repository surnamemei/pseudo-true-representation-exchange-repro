"""Independent N=31 interval certificate; floating search only proposes cells.

The interval predicates and local searches below decide every assertion.  No
Stage-2 (N=21) certificate data are used as a premise.
"""
from __future__ import annotations

import csv
import json
from decimal import Decimal, getcontext
from pathlib import Path

import numpy as np
from mpmath import iv, mp

import three_to_two_tone_branch_study as study
import three_to_two_tone_stage2_interval as interval
from three_to_two_tone_stage1_global import adaptive
from three_to_two_tone_stage2_inner_global import process
from three_to_two_tone_stage3_full_crossings import crossing, D

ROOT = Path(__file__).resolve().parent
N = 31
iv.dps = 90
mp.dps = 100
getcontext().prec = 110
interval.N = N


def bounds(x):
    return [str(x.a), str(x.b)]


def sup(x):
    return abs(x).b


def local_check(center, eps, radius='1e-10'):
    """Krawczyk existence/uniqueness and positive Hessian on one box."""
    e = iv.mpf(eps)
    x0 = [iv.mpf(v) for v in center]
    r = iv.mpf(radius)
    X = [v + iv.mpf([-r.b, r.b]) for v in x0]
    j0 = interval.objective_jet(*x0, e)
    jx = interval.objective_jet(*X, e)
    H = jx.H
    hf = np.array([[float(j0.H[i][j].mid) for j in range(2)] for i in range(2)])
    cf = np.linalg.inv(hf)
    C = [[iv.mpf(repr(float(cf[i, j]))) for j in range(2)] for i in range(2)]
    M = [[iv.mpf(int(i == j)) - sum(C[i][k]*H[k][j] for k in range(2))
          for j in range(2)] for i in range(2)]
    Y = [X[i]-x0[i] for i in range(2)]
    K = [x0[i]-sum(C[i][k]*j0.g[k] for k in range(2))
         +sum(M[i][k]*Y[k] for k in range(2)) for i in range(2)]
    included = all(K[i].a > X[i].a and K[i].b < X[i].b for i in range(2))
    contract = max(sum(sup(M[i][j]) for j in range(2)) for i in range(2))
    det = H[0][0]*H[1][1]-H[0][1]*H[1][0]
    positive = H[0][0].a > 0 and det.a > 0
    J = j0.v
    for i in range(2):
        J += j0.g[i]*Y[i]
    for i in range(2):
        for k in range(2):
            J += iv.mpf('0.5')*H[i][k]*Y[i]*Y[k]
    out = dict(center=center, radius_u=radius, X=[bounds(v) for v in X],
               K=[bounds(v) for v in K], inclusion=bool(included),
               contraction_upper=str(contract), contraction_lt_one=bool(contract < 1),
               H11_lower=str(H[0][0].a), Hessian_det_lower=str(det.a),
               positive_Hessian=bool(positive), J=bounds(J))
    return out, J, X


def dJ_de(u1, u2, e):
    def di(v):
        return iv.mpf(1)+sum(2*iv.cos(iv.mpf(k)*v/N) for k in range(1, (N-1)//2+1))
    g=di(u2-u1); den=iv.mpf(N*N)-g*g
    p1=di(u1+2)+di(u1-2); p2=di(u2+2)+di(u2-2)
    z1=di(u1-10); z2=di(u2-10)
    h1=p1-e*z1; h2=p2-e*z2
    Eder=sum(2*e-4*iv.cos(iv.mpf(2*t)/N)*iv.cos(iv.mpf(10*t)/N)
             for t in range(-(N-1)//2, (N-1)//2+1))
    capder=-2*N*(h1*z1+h2*z2)+2*g*(z1*h2+h1*z2)
    return Eder-capder/den


def gram_amplitude_bounds(X, eps):
    """Interval 2-tone Gram determinant and reduced LS amplitudes."""
    def di(v):
        return iv.mpf(1)+sum(2*iv.cos(iv.mpf(k)*v/N)
                             for k in range(1,(N-1)//2+1))
    u1,u2=X;e=iv.mpf(eps)
    g=di(u2-u1);den=iv.mpf(N*N)-g*g
    h1=di(u1+2)+di(u1-2)-e*di(u1-10)
    h2=di(u2+2)+di(u2-2)-e*di(u2-10)
    a1=(N*h1-g*h2)/den;a2=(N*h2-g*h1)/den
    return dict(Gram_determinant=bounds(den),Gram_positive=bool(den.a>0),
                amplitudes=[bounds(a1),bounds(a2)])


def center_projected(m, h, e, E):
    """Continuous confluent basis, including the exactly coalescent limit."""
    h1=iv.mpf(2)-e; h2=iv.mpf(0)
    c1sq=iv.mpf(1); c2sq=iv.mpf(0)
    w1psq=iv.mpf(0); w2psq=iv.mpf(0)
    for t in range(1, (N-1)//2+1):
        tt=iv.mpf(t); z=h*tt
        co=iv.cos(z); si=iv.sin(z)
        sc=iv.sin(z)/z
        sp=(z*iv.cos(z)-iv.sin(z))/(z*z)
        p=iv.cos((iv.mpf(2)/N-m)*tt)+iv.cos((-iv.mpf(2)/N-m)*tt)
        p-=e*iv.cos((iv.mpf(10)/N-m)*tt)
        q=iv.sin((iv.mpf(2)/N-m)*tt)+iv.sin((-iv.mpf(2)/N-m)*tt)
        q-=e*iv.sin((iv.mpf(10)/N-m)*tt)
        h1+=2*co*p; h2+=2*tt*sc*q
        c1sq+=2*co*co; c2sq+=2*tt*tt*sc*sc
        w1psq+=2*tt*tt*si*si; w2psq+=2*tt**4*sp*sp
    J=E-h1*h1/c1sq-h2*h2/c2sq
    return J,iv.sqrt(c1sq),iv.sqrt(c2sq),iv.sqrt(w1psq),iv.sqrt(w2psq)


def directed_cell(cell,e,E,incumbent,inflate=Decimal('1e-12')):
    dep,mlo,mhi,hlo,hhi,*_=cell
    a,b,c,d=map(lambda v:Decimal(repr(float(v))), (mlo,mhi,hlo,hhi))
    mc=(a+b)/2;hc=(c+d)/2
    rm=(b-a)/2+inflate;rh=(d-c)/2+inflate
    m=iv.mpf(str(mc));h=iv.mpf(str(hc))
    rmi=iv.mpf(str(rm));rhi=iv.mpf(str(rh))
    Jc,c1,c2,w1p,w2p=center_projected(m,h,e,E)
    tmax=(N-1)//2
    t2=iv.sqrt(iv.mpf(sum(t**4 for t in range(-tmax,tmax+1))))
    t3=iv.sqrt(iv.mpf(sum(t**6 for t in range(-tmax,tmax+1))))
    D1=w1p*rhi+iv.mpf('0.5')*t2*rhi*rhi
    D2=w2p*rhi+t3*rhi*rhi/6
    if not (c1.a>D1.b and c2.a>D2.b):
        return False,None,'basis_denominator_nonpositive'
    theta=iv.sqrt(2)*tmax*rmi+iv.sqrt((D1/(c1-D1))**2+(D2/(c2-D2))**2)
    if Jc.a <= 0:
        return False,None,'center_J_nonpositive'
    residual=iv.sqrt(Jc).a-iv.sqrt(E).b*theta.b
    bound=iv.mpf(0) if residual<=0 else residual*residual
    return bool(bound.a>incumbent),bound,'ok' if bound.a>incumbent else 'bound_not_above_incumbent'


def accepted_contained(cell, centers, inflate=Decimal('1e-12')):
    _,mlo,mhi,hlo,hhi,*_=cell
    a,b,c,d=map(lambda v:Decimal(repr(float(v))), (mlo,mhi,hlo,hhi))
    a-=inflate;b+=inflate;c-=inflate;d+=inflate
    rad=Decimal('2.00000001')
    # Both ordered endpoint-frequency coordinates must lie inside one inner box.
    u1lo=N*(a-d);u1hi=N*(b-c);u2lo=N*(a+c);u2hi=N*(b+d)
    for center in centers:
        q1,q2=map(lambda x:Decimal(repr(float(x))), center)
        if (u1lo>=q1-rad and u1hi<=q1+rad and
                u2lo>=q2-rad and u2hi<=q2+rad):
            return True
    return False


def main():
    seed=(-4.732308155500305,0.6135741125336811,
          -0.05731377959637947,12.440325828105324,0.24914308985250178)
    root=crossing('2',N=N,b=10,seed=seed)
    centers=[[root['u_A1'],root['u_A2']],[root['u_B1'],root['u_B2']]]
    ec=mp.mpf(root['epsilon_cross'])
    er=[mp.nstr(ec-mp.mpf('1e-14'),90), mp.nstr(ec+mp.mpf('1e-14'),90)]
    eps=iv.mpf(er)
    checks=[];endpoint_rows=[]
    for label,e in [('lower',er[0]),('upper',er[1]),('bracket',er)]:
        ja,Ja,Xa=local_check(centers[0],e)
        jb,Jb,Xb=local_check(centers[1],e)
        diff=Ja-Jb
        sign='negative' if diff.b<0 else 'positive' if diff.a>0 else 'undetermined'
        checks.append(dict(label=label,A=ja,B=jb,DeltaJ=bounds(diff),sign=sign))
        endpoint_rows.append(dict(endpoint=label,epsilon=str(e),DeltaJ_lower=str(diff.a),
                                  DeltaJ_upper=str(diff.b),sign=sign,
                                  A_Krawczyk=ja['inclusion'],B_Krawczyk=jb['inclusion'],
                                  A_H_positive=ja['positive_Hessian'],
                                  B_H_positive=jb['positive_Hessian']))
    with (ROOT/'N31_crossing_certificate.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(endpoint_rows[0]));w.writeheader();w.writerows(endpoint_rows)
    slope=dJ_de(*Xa,eps)-dJ_de(*Xb,eps)
    local_valid=all(c[k]['inclusion'] and c[k]['positive_Hessian']
                    for c in checks for k in ('A','B'))
    bracket_valid=(checks[0]['sign']=='negative' and checks[1]['sign']=='positive'
                   and slope.a>0 and local_valid)
    print('N31 local',local_valid,'bracket',bracket_valid,'slope',bounds(slope),flush=True)
    if not bracket_valid:
        raise RuntimeError('N31 local crossing bracket did not validate')
    _,ja_bracket,_=local_check(centers[0],er)
    _,jb_bracket,_=local_check(centers[1],er)
    incumbent=min(ja_bracket.b,jb_bracket.b)
    large=[]
    for center,radius in zip(centers,('0.0003','0.01')):
        q,_,X=local_check(center,er,radius)
        q.update(gram_amplitude_bounds(X,er))
        large.append(q)
    large_valid=all(q['inclusion'] and q['contraction_lt_one'] and
                    q['positive_Hessian'] and q['Gram_positive'] for q in large)
    print('large-box certification',large_valid,flush=True)
    if not large_valid:
        raise RuntimeError('N31 large root boxes did not validate')
    # Full-domain float search proposes a partition, then all leaves are
    # independently checked with directed interval arithmetic.
    iv.dps=70
    (ua,_),(ub,_),x=study.two_branches(N,2,10,float(ec),np.pi)
    threshold=float(incumbent)+1e-6
    complete,summary,leaves,un=adaptive(x,(ua,ub),threshold,N,
                                        radius_u=2.,maxdepth=40,write_rows=True)
    print('proposal',complete,'leaves',len(leaves),'unresolved',len(un),flush=True)
    E=interval.data_energy(eps)
    rows=[];failed=[];minmargin=None
    for idx,cell in enumerate(leaves):
        depth,mlo,mhi,hlo,hhi,jf,lbf,status=cell
        if status=='excluded':
            ok,bound,reason=directed_cell(cell,eps,E,incumbent)
            if ok:
                margin=bound.a-incumbent
                minmargin=margin if minmargin is None else min(minmargin,margin)
            else: failed.append((idx,reason))
            ilb='' if bound is None else str(bound.a)
        elif status in ('near_A','near_B'):
            ok=accepted_contained(cell,(ua,ub))
            reason='inside_inner_box' if ok else 'accepted_not_contained'
            ilb=''
            if not ok:failed.append((idx,reason))
        else:
            ok=False;reason='proposal_unresolved';ilb='';failed.append((idx,reason))
        rows.append(dict(index=idx,depth=depth,m_lo=repr(float(mlo)),m_hi=repr(float(mhi)),
                         h_lo=repr(float(hlo)),h_hi=repr(float(hhi)),status=status,
                         proposed_J_lower=repr(float(lbf)),interval_J_lower=ilb,
                         verified=ok,reason=reason))
        if idx and idx%4000==0:print('verified',idx,'/',len(leaves),flush=True)
    with (ROOT/'N31_interval_cells.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    print('outer failures',len(failed),'minmargin',minmargin,flush=True)
    iv.dps=45
    inner=[]
    for which,centerfloat,centerhigh in [('A',ua,centers[0]),('B',ub,centers[1])]:
        valid,visited,cells,remaining=process('crossing',which,er,centerfloat,incumbent,
                                              max_cells=50000,
                                              known_center_override=centerhigh,
                                              local_radius_override='0.000299999999' if which=='A' else '0.009999999')
        inner.append(dict(branch=which,complete=valid,visited=visited,
                          terminal=len(cells),remaining=remaining,
                          status_counts={s:sum(r['status']==s for r in cells)
                                         for s in set(r['status'] for r in cells)}))
        print('inner',inner[-1],flush=True)
    certificate=dict(N=N,z=2,geometry=dict(strong_frequencies_u=[-2,2],weak_frequency_u=10,
                                            physical_frequencies_omega=['-2/31','2/31','10/31'],
                                            weak_phase='pi',strong_amplitudes=[1,1]),
                     high_precision_numerical_root=root,epsilon_bracket=er,
                     crossing_checks=checks,large_root_boxes=large,
                     slope_interval=bounds(slope),
                     crossing_bracket_certified=bracket_valid,
                     outer=dict(proposal_complete=complete,leaf_count=len(leaves),
                                excluded=sum(r['status']=='excluded' for r in rows),
                                accepted=sum(r['status'] in ('near_A','near_B') for r in rows),
                                failed=len(failed),first_failures=failed[:30],
                                min_objective_margin=str(minmargin),
                                complete=complete and not failed),inner=inner,
                     full_global_certificate=bool(bracket_valid and large_valid and complete and not failed
                                                  and all(v['complete'] for v in inner)))
    (ROOT/'N31_global_certificate.json').write_text(json.dumps(certificate,indent=2),encoding='utf-8')
    print('CERTIFICATE',certificate['full_global_certificate'],flush=True)


if __name__=='__main__':main()
