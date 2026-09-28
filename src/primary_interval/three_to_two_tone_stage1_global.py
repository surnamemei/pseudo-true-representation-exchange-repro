"""Adaptive projection-distance lower bound on the full unordered frequency torus.

The bound is analytic but evaluated in ordinary float64, so this is a
numerically certified-style exclusion, not a formal outward-rounded proof.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

import three_to_two_tone_branch_study as st
import three_to_two_tone_explore as ex

ROOT=Path(__file__).resolve().parent


def sinc_prime(z):
    z=np.asarray(z)
    out=np.empty_like(z,dtype=float)
    small=np.abs(z)<1e-3
    q=z[small]
    out[small]=-q/3+q**3/30-q**5/840
    q=z[~small]
    out[~small]=(q*np.cos(q)-np.sin(q))/q**2
    return out


def evaluate(cells,x,n=21,chunk=12000):
    """Return J(center) and rigorous-form lower bounds for each m,h cell."""
    s=np.arange(n)-(n-1)/2
    t2norm=np.linalg.norm(s*s);t3norm=np.linalg.norm(s**3)
    xnorm=np.linalg.norm(x);xnorm2=xnorm*xnorm
    mlo,mhi,hlo,hhi=cells.T
    mc=(mlo+mhi)/2;hc=(hlo+hhi)/2
    rm=(mhi-mlo)/2;rh=(hhi-hlo)/2
    outJ=np.empty(len(cells));outL=np.empty(len(cells))
    for k in range(0,len(cells),chunk):
        sl=slice(k,min(k+chunk,len(cells)))
        hm=hc[sl,None]*s[None,:]
        w1=np.cos(hm)
        w2=s[None,:]*np.sinc(hm/np.pi)
        c1=np.linalg.norm(w1,axis=1)
        c2=np.linalg.norm(w2,axis=1)
        phase=np.exp(1j*mc[sl,None]*s[None,:])
        q1=phase*w1/c1[:,None]
        q2=phase*w2/c2[:,None]
        h1=q1.conj()@x
        h2=q2.conj()@x
        jj=np.maximum(0,xnorm2-np.abs(h1)**2-np.abs(h2)**2)
        outJ[sl]=jj
        w1p=-s[None,:]*np.sin(hm)
        w2p=s[None,:]**2*sinc_prime(hm)
        dw1=np.linalg.norm(w1p,axis=1)*rh[sl]+.5*t2norm*rh[sl]**2
        dw2=np.linalg.norm(w2p,axis=1)*rh[sl]+t3norm/6*rh[sl]**2
        with np.errstate(divide='ignore',invalid='ignore'):
            hh=np.sqrt((dw1/(c1-dw1))**2+(dw2/(c2-dw2))**2)
        hh=np.where((dw1<c1)&(dw2<c2),hh,np.inf)
        # Q is orthonormal (even/odd centered vectors), hence
        # ||P-Pc|| <= ||Q-Qc||_F. The m term uses |t|<= (n-1)/2.
        theta=np.sqrt(2)*(n-1)/2*rm[sl]+hh
        outL[sl]=np.maximum(0,np.sqrt(jj)-xnorm*theta-1e-12)**2
    return outJ,outL


def inside_neighborhood(cells,centers,n,radius_u):
    m=(cells[:,0]+cells[:,1])/2*n
    h=(cells[:,2]+cells[:,3])/2*n
    rad=((cells[:,1]-cells[:,0])+(cells[:,3]-cells[:,2]))/2*n
    u1=m-h;u2=m+h
    hit=np.full(len(cells),-1,dtype=int)
    for j,(a,b) in enumerate(centers):
        good=(np.abs(u1-a)+rad<=radius_u)&(np.abs(u2-b)+rad<=radius_u)
        hit[good]=j
    return hit


def adaptive(x,branches,threshold,n=21,radius_u=3.0,maxdepth=24,
             max_active=1200000,write_rows=True):
    cells=np.array([[-np.pi,np.pi,0,np.pi/2]],dtype=float)
    depths=np.array([0],dtype=int)
    summary=[];leaves=[]
    completed=True
    for dep in range(maxdepth+1):
        jj,lb=evaluate(cells,x,n)
        hit=inside_neighborhood(cells,branches,n,radius_u)
        excluded=lb>threshold
        accepted=(hit>=0)&~excluded
        unresolved=~excluded&~accepted
        if write_rows:
            for ix in np.flatnonzero(excluded|accepted|((dep==maxdepth)&unresolved)):
                leaves.append((dep,*cells[ix],jj[ix],lb[ix],
                               'excluded' if excluded[ix] else
                               ('near_A' if hit[ix]==0 else 'near_B') if accepted[ix] else 'unresolved'))
        summary.append((dep,len(cells),int(excluded.sum()),int(accepted.sum()),
                        int(unresolved.sum()),float(np.min(lb)),float(np.max(lb))))
        print('depth',dep,'cells',len(cells),'excluded',excluded.sum(),
              'near',accepted.sum(),'unresolved',unresolved.sum(),flush=True)
        if not unresolved.any():break
        if dep==maxdepth or int(unresolved.sum())*2>max_active:
            completed=False;break
        old=cells[unresolved]
        rm=(old[:,1]-old[:,0])/2;rh=(old[:,3]-old[:,2])/2
        split_m=np.sqrt(2)*(n-1)/2*rm>=13*rh
        child=np.repeat(old,2,axis=0)
        midm=(old[:,0]+old[:,1])/2
        midh=(old[:,2]+old[:,3])/2
        ids=np.flatnonzero(split_m)
        child[2*ids,1]=midm[ids]
        child[2*ids+1,0]=midm[ids]
        ids=np.flatnonzero(~split_m)
        child[2*ids,3]=midh[ids]
        child[2*ids+1,2]=midh[ids]
        cells=child
    return completed,summary,leaves,cells[unresolved] if not completed else np.empty((0,4))


def objective_crosscheck(n=21,count=3000):
    st.set_n(n)
    x=ex.record(2,10,st.crossing(n,2,10,np.pi)[0],np.pi)
    rng=np.random.default_rng(19)
    m=rng.uniform(-np.pi,np.pi,count)
    h=rng.uniform(0,np.pi/2,count)
    h[:20]=np.geomspace(1e-10,1e-3,20)
    cells=np.column_stack((m,m,h,h))
    jj,_=evaluate(cells,x,n)
    direct=np.array([ex.cost([mm-hh,mm+hh],x) for mm,hh in zip(m,h)])
    return float(np.max(np.abs(jj-direct))),float(np.quantile(np.abs(jj-direct),.99))


if __name__=='__main__':
    n,d,b,phase=21,2,10,np.pi
    ec=st.crossing(n,d,b,phase)[0]
    print('projection versus QR',objective_crosscheck(),flush=True)
    rows=[]
    for label,eps in [('below',ec-.035),('crossing',ec),('above',ec+.035)]:
        (ua,ja),(ub,jb),x=st.two_branches(n,d,b,eps,phase)
        threshold=min(ja,jb)+.4
        complete,summary,leaves,unresolved=adaptive(
            x,(ua,ub),threshold,n,radius_u=5.0,maxdepth=24,
            max_active=1200000)
        print(label,'complete',complete,'leaves',len(leaves),
              'unresolved',len(unresolved),flush=True)
        for z in leaves:
            depth,mlo,mhi,hlo,hhi,jc,lb,status=z
            rows.append(dict(epsilon_label=label,epsilon=eps,depth=depth,
                             m_lo=mlo,m_hi=mhi,h_lo=hlo,h_hi=hhi,
                             J_center=jc,J_lower_bound=lb,
                             threshold=threshold,status=status))
        if not complete:
            for q in unresolved:
                rows.append(dict(epsilon_label=label,epsilon=eps,depth=-1,
                                 m_lo=q[0],m_hi=q[1],h_lo=q[2],h_hi=q[3],
                                 J_center=np.nan,J_lower_bound=np.nan,
                                 threshold=threshold,status='unresolved_budget'))
    with (ROOT/'globality_cells.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
