"""Exploratory tangent landscape; certification is separate."""
import csv
import numpy as np
from scipy.optimize import minimize_scalar
from mpmath import mp
from pathlib import Path

N=21
s=np.arange(-10,11,dtype=float)/N
S2=sum(s*s);S4=sum(s**4);S6=sum(s**6)
def D(x,k=0):
    t=s.reshape((1,-1))
    x=np.asarray(x).reshape((-1,1))
    q=x*t
    return np.sum(((-1)**(k//2)*t**k*np.cos(q)) if k%2==0 else ((-1)**((k+1)//2)*t**k*np.sin(q)),axis=1)
def cost(v,l):
    v=np.atleast_1d(v)
    d=D(v);dp=D(v,1);d2=D(v,2)
    q0=-S2-l*D([10])[0];q1=2*l*D([10],1)[0]
    C=S4-2*l*D([10],2)[0]+N*l*l-q0*q0/N-q1*q1/(4*S2)
    A=N-d*d/N-dp*dp/S2
    B=d2-l*D(10-v)-d*q0/N+dp*q1/(2*S2)
    out=C-B*B/A
    out[np.abs(v)<.03]=np.nan
    return out
mp.dps=70
ss=[mp.mpf(t)/N for t in range(-10,11)]
def dm(x,k=0):
    return mp.fsum((1j*t)**k*mp.e**(1j*x*t) for t in ss).real
SM2=mp.fsum(t*t for t in ss);SM4=mp.fsum(t**4 for t in ss);SM6=mp.fsum(t**6 for t in ss)
def cm(v,l):
    v=mp.mpf(v);l=mp.mpf(l)
    d=dm(v);dp=dm(v,1)
    q0=-SM2-l*dm(10);q1=2*l*dm(10,1)
    C=SM4-2*l*dm(10,2)+N*l*l-q0*q0/N-q1*q1/(4*SM2)
    A=N-d*d/N-dp*dp/SM2
    B=dm(v,2)-l*dm(10-v)-d*q0/N+dp*q1/(2*SM2)
    return C-B*B/A
M=dm(10,3)+SM4/SM2*dm(10,1)
V3=SM6-SM4*SM4/SM2
K=3*M/V3
V2=SM4-SM2*SM2/N
Q=dm(10,2)+SM2*dm(10)/N
def coal(l):
    l=mp.mpf(l)
    q0=-SM2-l*dm(10);q1=2*l*dm(10,1)
    qnorm=SM4-2*l*dm(10,2)+N*l*l
    q2=-SM4+l*dm(10,2)-SM2*q0/N
    return qnorm-q0*q0/N-q1*q1/(4*SM2)-q2*q2/V2
def strong(l):return mp.mpf(l)**2*(N-dm(10)**2/N-dm(10,1)**2/SM2)
def main():
    print('K',mp.nstr(K,25),'M',mp.nstr(M,25),'V3',mp.nstr(V3,25),'lambda_gamma_zero',mp.nstr(-V2/Q,25),flush=True)
    rows=[];grid=np.linspace(-66,66,2641)
    for l in [0,.00001,.0001,.001,.002,.005,.01,.015,.02,.025,.03,.035,.04,.045,.05,.055,.06,.065,.0659804111269914,.07,.075,.08,.09,.1,.12,.15,.2]:
        cc=cost(grid,l)
        cand=[]
        for i in range(1,len(grid)-1):
            if np.isfinite(cc[i]) and cc[i]<cc[i-1] and cc[i]<cc[i+1]:
                z=minimize_scalar(lambda v:cost([v],l)[0],bounds=(grid[i-1],grid[i+1]),method='bounded',options={'xatol':1e-13})
                cand.append((float(z.fun),float(z.x)))
        if l>0:
            try:
                if l<.005:
                    v=mp.findroot(lambda x:mp.diff(lambda y:cm(y,l),x),(K*l*.8,K*l*1.2),tol=mp.mpf('1e-55'))
                    cand.append((float(cm(v,l)),float(v)))
            except Exception as e:print('small root failure',l,e,flush=True)
        cand=sorted(cand)[:4]
        r=dict(lam=l,C_coal=float(coal(l)),C_strong=float(strong(l)),global_cost=min(float(coal(l)),cand[0][0] if cand else float('inf')),
               winner='coal' if not cand or float(coal(l))<cand[0][0] else 'regular',v1=cand[0][1] if cand else '',R1=cand[0][0] if cand else '',v2=cand[1][1] if len(cand)>1 else '',R2=cand[1][0] if len(cand)>1 else '',v3=cand[2][1] if len(cand)>2 else '',R3=cand[2][0] if len(cand)>2 else '')
        rows.append(r);print(r,flush=True)
    with open('stage13_explore.csv','w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
if __name__=='__main__':main()
