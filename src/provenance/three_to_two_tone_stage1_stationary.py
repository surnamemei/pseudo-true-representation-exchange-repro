"""Stationary-point search for the centered two-tone profiled objective."""

import csv
from pathlib import Path

import numpy as np
from scipy.optimize import root

import three_to_two_tone_branch_study as st
import three_to_two_tone_explore as ex
from three_to_two_tone_stage1_global import sinc_prime

ROOT=Path(__file__).resolve().parent


def objective_gradient(mh,x,n=21):
    m,h=mh
    t=np.arange(n)-(n-1)/2
    zz=h*t
    w1=np.cos(zz)
    w2=t*np.sinc(zz/np.pi)
    d1=-t*np.sin(zz)
    d2=t*t*sinc_prime(zz)
    c1=np.linalg.norm(w1);c2=np.linalg.norm(w2)
    q1=np.exp(1j*m*t)*w1/c1
    q2=np.exp(1j*m*t)*w2/c2
    dh1=np.exp(1j*m*t)*(d1-w1*np.dot(w1,d1)/c1**2)/c1
    dh2=np.exp(1j*m*t)*(d2-w2*np.dot(w2,d2)/c2**2)/c2
    cc=np.array([np.vdot(q1,x),np.vdot(q2,x)])
    cm=np.array([np.vdot(1j*t*q1,x),np.vdot(1j*t*q2,x)])
    ch=np.array([np.vdot(dh1,x),np.vdot(dh2,x)])
    jj=float(np.vdot(x,x).real-np.sum(abs(cc)**2))
    grad=np.array([-2*np.vdot(cc,cm).real,-2*np.vdot(cc,ch).real])
    return jj,grad


def canonical(m,h):
    h=abs(h)%(np.pi)
    if h>np.pi/2:
        h=np.pi-h;m+=np.pi
    m=(m+np.pi)%(2*np.pi)-np.pi
    return np.array([m,h])


def hessian(mh,x,n=21,step=1e-5):
    a=np.array([step,0]);b=np.array([0,step])
    m1=(objective_gradient(mh+a,x,n)[1]-objective_gradient(mh-a,x,n)[1])/(2*step)
    m2=(objective_gradient(mh+b,x,n)[1]-objective_gradient(mh-b,x,n)[1])/(2*step)
    return (np.column_stack((m1,m2))+np.column_stack((m1,m2)).T)/2


def search(x,n=21,nm=64,nh=32):
    seeds=[(m,h) for m in np.linspace(-np.pi,np.pi,nm,endpoint=False)
                 for h in np.linspace(.025,np.pi/2-.025,nh)]
    found=[]
    for m,h in seeds:
        res=root(lambda q:objective_gradient(q,x,n)[1],[m,h],method='hybr',
                 options={'xtol':1e-10,'maxfev':140})
        if np.linalg.norm(objective_gradient(res.x,x,n)[1])>1e-6:
            continue
        q=canonical(*res.x)
        if q[1]<1e-5 or q[1]>np.pi/2-1e-5:
            continue
        if any(np.linalg.norm(q-v['q'])<2e-4 for v in found):
            continue
        jj,grad=objective_gradient(q,x,n)
        eig=np.linalg.eigvalsh(hessian(q,x,n))
        cl='minimum' if eig[0]>1e-6 else ('maximum' if eig[1]<-1e-6 else 'saddle')
        uu=np.sort(np.array([q[0]-q[1],q[0]+q[1]])*n)
        found.append({'q':q,'J':jj,'grad':float(np.linalg.norm(grad)),
                      'eig':eig,'class':cl,'u':uu})
    found.sort(key=lambda r:r['J'])
    return found,len(seeds)


if __name__=='__main__':
    n,d,b=21,2,10
    ec=st.crossing(n,d,b,np.pi)[0]
    rows=[]
    for label,eps in [('below',ec-.035),('crossing',ec),('above',ec+.035)]:
        st.set_n(n);x=ex.record(d,b,eps,np.pi)
        out,seeds=search(x,n)
        print(label,'seeds',seeds,'roots',len(out),
              'minima',sum(r['class']=='minimum' for r in out),
              'best',[(r['class'],r['J'],r['u']) for r in out[:5]],flush=True)
        for k,r in enumerate(out):
            rows.append(dict(epsilon_label=label,epsilon=eps,rank_by_J=k+1,
                             m=r['q'][0],h=r['q'][1],u1=r['u'][0],u2=r['u'][1],
                             J=r['J'],grad_norm=r['grad'],
                             hess_eig_min=r['eig'][0],hess_eig_max=r['eig'][1],
                             classification=r['class'],deterministic_seeds=seeds))
    with (ROOT/'stationary_points.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
