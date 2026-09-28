"""Symmetric midpoint expansion and numerical finite-N validation."""
from __future__ import annotations
import csv
import json
from pathlib import Path
from mpmath import mp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from stage14_continuum import D as Dinf, R as Rinf, MU2, MU4

ROOT=Path(__file__).resolve().parent
mp.dps=70
ref=json.loads((ROOT/'stage14_continuum_candidate.json').read_text(encoding='utf-8'))
va0=mp.mpf(ref['v_A']);vb0=mp.mpf(ref['v_B']);l0=mp.mpf(ref['lambda'])

def Dt(q,t,k=0):
    f=lambda x:Dinf(x)*(1+x*x*t/24+7*x**4*t*t/5760)
    return mp.diff(f,q,k) if k else f(q)

def Rt(v,l,t):
    m2=mp.mpf(1)/12-t/12
    m4=mp.mpf(1)/80-t/24+7*t*t/240
    d=Dt(v,t);dp=Dt(v,t,1)
    q0=-m2-l*Dt(10,t);q1=2*l*Dt(10,t,1)
    A=1-d*d-dp*dp/m2
    B=Dt(v,t,2)-l*Dt(10-v,t)-d*q0+dp*q1/(2*m2)
    C=m4-2*l*Dt(10,t,2)+l*l-q0*q0-q1*q1/(4*m2)
    return C-B*B/A

def DN(q,N,k=0):
    q=mp.mpf(q);N=mp.mpf(N)
    def d(x):
        if abs(x)<mp.mpf('1e-30'):return mp.mpf(1)
        return mp.sin(x/2)/(N*mp.sin(x/(2*N)))
    return mp.diff(d,q,k) if k else d(q)

def RN(v,l,N):
    n=mp.mpf(N);m2=(1-1/n**2)/12
    m4=(n*n-1)*(3*n*n-7)/(240*n**4)
    d=DN(v,N);dp=DN(v,N,1)
    q0=-m2-l*DN(10,N);q1=2*l*DN(10,N,1)
    A=1-d*d-dp*dp/m2
    B=DN(v,N,2)-l*DN(10-v,N)-d*q0+dp*q1/(2*m2)
    C=m4-2*l*DN(10,N,2)+l*l-q0*q0-q1*q1/(4*m2)
    return C-B*B/A

def crossN(N):
    fa=lambda va,vb,l:mp.diff(lambda v:RN(v,l,N),va)
    fb=lambda va,vb,l:mp.diff(lambda v:RN(v,l,N),vb)
    fc=lambda va,vb,l:RN(va,l,N)-RN(vb,l,N)
    return mp.findroot((fa,fb,fc),(va0,vb0,l0),tol=mp.mpf('1e-50'),maxsteps=40)

def main():
    LA=lambda v:Rt(v,l0,mp.mpf(0))
    F=lambda l,t:Rt(va0,l,t)-Rt(vb0,l,t)
    Fl=mp.diff(lambda l:F(l,0),l0)
    Ft=mp.diff(lambda t:F(l0,t),0)
    a=-Ft/Fl
    def gj(v):
        return (mp.diff(lambda l:mp.diff(lambda x:Rt(x,l,0),v),l0)*a
                +mp.diff(lambda t:mp.diff(lambda x:Rt(x,l0,t),v),0))
    ga=gj(va0);gb=gj(vb0)
    ha=mp.diff(lambda v:Rt(v,l0,0),va0,2)
    hb=mp.diff(lambda v:Rt(v,l0,0),vb0,2)
    ua=-ga/ha;ub=-gb/hb
    Fll=mp.diff(lambda l:F(l,0),l0,2)
    Flt=mp.diff(lambda l:mp.diff(lambda t:F(l,t),0),l0)
    Ftt=mp.diff(lambda t:F(l0,t),0,2)
    bb=-(Fll*a*a/2+Flt*a+Ftt/2-ga*ga/(2*ha)+gb*gb/(2*hb))/Fl
    theory={
      'lambda_infinity':mp.nstr(l0,60),
      'a_1_over_N2':mp.nstr(a,55),
      'b_1_over_N4':mp.nstr(bb,55),
      'v_A_1_over_N2':mp.nstr(ua,40),
      'v_B_1_over_N2':mp.nstr(ub,40),
      'method':'midpoint Euler-Maclaurin/Dirichlet expansion and implicit crossing derivatives',
      'status':'analytical expansion conditional on continuum nondegenerate global crossing; coefficients numerically evaluated'
    }
    (ROOT/'stage14_asymptotic_coefficients.json').write_text(json.dumps(theory,indent=2),encoding='utf-8')
    finite=json.loads((ROOT/'stage12_generalN'/'N21_reference.json').read_text(encoding='utf-8')) if False else None
    known={}
    with (ROOT/'generalN_certified_instances.csv').open(encoding='utf-8') as f:
        for row in csv.DictReader(f): known[int(row['N'])]=row['lambda_N']
    rows=[]
    for N in [11,15,21,31,41,101,201,501,1001]:
        if N in known:
            l=mp.mpf(known[N]);origin='existing certified'
        else:
            va,vb,l=crossN(N);origin='new numerical validation'
        t=mp.mpf(1)/N**2
        pred2=l0+a*t;pred4=pred2+bb*t*t
        rows.append({'N':N,'lambda_N':mp.nstr(l,42),'origin':origin,
                     'lambda_N_minus_lambda_inf':mp.nstr(l-l0,30),
                     'N2_scaled_difference':mp.nstr((l-l0)*N*N,30),
                     'N2_prediction':mp.nstr(pred2,30),
                     'N4_prediction':mp.nstr(pred4,30),
                     'N4_prediction_error':mp.nstr(l-pred4,25)})
        print('N',N,'lambda',mp.nstr(l,25),flush=True)
    with (ROOT/'lambdaN_asymptotics.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    fig,ax=plt.subplots(figsize=(7.1,4.1))
    xx=np.array([1/r['N']**2 for r in rows])
    yy=np.array([float(r['lambda_N']) for r in rows])
    grid=np.linspace(0,max(xx),300)
    ax.plot(grid,float(l0)+float(a)*grid+float(bb)*grid*grid,color='#222222',label=r'$\lambda_\infty+a/N^2+b/N^4$')
    ax.scatter(xx[:5],yy[:5],color='#1464a0',s=36,label='interval-certified finite-N crossings')
    ax.scatter(xx[5:],yy[5:],color='#d06322',s=36,label='numerical validation')
    ax.scatter([0],[float(l0)],marker='*',s=110,color='black',label='continuum limit')
    ax.set(xlabel=r'$1/N^2$',ylabel=r'critical $\lambda_N$')
    ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.tight_layout();fig.savefig(ROOT/'continuum_vs_finiteN_figure.pdf');plt.close(fig)
    print(json.dumps(theory,indent=2),flush=True)

if __name__=='__main__':main()
