"""Quadratic-in-epsilon outer-cell experiment, with directed interval coefficients.

This is a diagnostic until the complete outer and inner covers pass.
"""
import csv
import json
import re
import sys
from decimal import Decimal, getcontext
from pathlib import Path

from mpmath import iv, mp

import three_to_two_tone_stage2_interval as it

ROOT = Path(__file__).resolve().parent
iv.dps = 70
mp.dps = 80
getcontext().prec = 90
it.N = 21

def interval(x):
    return [str(mp.mpf(x.a)), str(mp.mpf(x.b))]

def num(s):
    return mp.mpf(re.findall(r"[-+]?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?", s)[0])

def coefficients(eval_at):
    m = eval_at(iv.mpf(-1))
    z = eval_at(iv.mpf(0))
    p = eval_at(iv.mpf(1))
    return ((p+m-2*z)/2, (p-m)/2, z)

def val(c, e):
    return (c[0]*e+c[1])*e+c[2]

def qmin(c, lo, hi):
    lo, hi = iv.mpf(str(lo)), iv.mpf(str(hi))
    mid=(lo+hi)/2
    deriv=2*c[0]*mid+c[1]
    whole_deriv=2*c[0]*(lo+(hi-lo)*iv.mpf([0,1]))+c[1]
    if whole_deriv.a>0:
        return val(c,lo).a
    if whole_deriv.b<0:
        return val(c,hi).a
    rad=(hi-lo).b/2
    return (val(c,mid)-abs(deriv).b*rad-abs(c[0]).b*rad*rad).a

def qmax(c, lo, hi):
    return -qmin(tuple(-a for a in c),lo,hi)

E=coefficients(it.data_energy)
ref=json.loads((ROOT/'stage2_reference.json').read_text())
ec0,ec1=json.loads((ROOT/'stage2_crossing_global.json').read_text())['epsilon_bracket']
ec=(mp.mpf(ec0)+mp.mpf(ec1))/2
epsrad=mp.mpf(sys.argv[1]) if len(sys.argv)>1 else mp.mpf('0.0025')
lo,hi=ec-epsrad,ec+epsrad

incumbents={}
for branch in ('A','B'):
    u=ref[branch+'_at_crossing']
    incumbent=coefficients(lambda e:it.objective_gradient_interval(iv.mpf(str(u[0])),iv.mpf(str(u[1])),e)[0])
    incumbents[branch]=incumbent

def sinc(z): return iv.sin(z)/z
def sinc_prime(z): return (z*iv.cos(z)-iv.sin(z))/(z*z)

def cell_coeff(row):
    a,b,c,d=(Decimal(row[k]) for k in ('m_lo','m_hi','h_lo','h_hi'))
    m=iv.mpf(str((a+b)/2)); h=iv.mpf(str((c+d)/2))
    rm=iv.mpf(str((b-a)/2+Decimal('1e-12')))
    rh=iv.mpf(str((d-c)/2+Decimal('1e-12')))
    h10=iv.mpf(2);h11=iv.mpf(-1)
    h20=iv.mpf(0);h21=iv.mpf(0)
    c1sq=iv.mpf(1);c2sq=iv.mpf(0)
    w1psq=iv.mpf(0);w2psq=iv.mpf(0)
    for k in range(1,11):
        t=iv.mpf(k);z=h*t
        co=iv.cos(z);si=iv.sin(z)
        sc=sinc(z);sp=sinc_prime(z)
        p=iv.cos((iv.mpf(2)/21-m)*t)+iv.cos((-iv.mpf(2)/21-m)*t)
        q=iv.sin((iv.mpf(2)/21-m)*t)+iv.sin((-iv.mpf(2)/21-m)*t)
        pw=iv.cos((iv.mpf(10)/21-m)*t)
        qw=iv.sin((iv.mpf(10)/21-m)*t)
        h10+=2*co*p;h11-=2*co*pw
        h20+=2*t*sc*q;h21-=2*t*sc*qw
        c1sq+=2*co*co;c2sq+=2*t*t*sc*sc
        w1psq+=2*t*t*si*si;w2psq+=2*t**4*sp*sp
    c=(E[2]-h10*h10/c1sq-h20*h20/c2sq)
    b=(E[1]-2*h10*h11/c1sq-2*h20*h21/c2sq)
    a=(E[0]-h11*h11/c1sq-h21*h21/c2sq)
    t2=iv.sqrt(iv.mpf(sum(t**4 for t in range(-10,11))))
    t3=iv.sqrt(iv.mpf(sum(t**6 for t in range(-10,11))))
    D1=iv.sqrt(w1psq)*rh+iv.mpf('0.5')*t2*rh*rh
    D2=iv.sqrt(w2psq)*rh+t3*rh*rh/6
    if not (iv.sqrt(c1sq).a>D1.b and iv.sqrt(c2sq).a>D2.b):
        return None,None
    theta=iv.sqrt(2)*10*rm+iv.sqrt((D1/(iv.sqrt(c1sq)-D1))**2+(D2/(iv.sqrt(c2sq)-D2))**2)
    return (a,b,c),theta

rows=[r for r in csv.DictReader((ROOT/'interval_boxes.csv').open(encoding='utf-8'))
      if r['label']=='crossing' and r['status']=='excluded']
rows.sort(key=lambda r:num(r['interval_J_lower']))
limit=int(sys.argv[2]) if len(sys.argv)>2 else len(rows)
out=[]
for i,row in enumerate(rows[:limit]):
    coef,theta=cell_coeff(row)
    if coef is None:
        out.append(dict(rank=i,passed=False,reason='denominator',cell={k:row[k] for k in ('m_lo','m_hi','h_lo','h_hi')}))
        continue
    Jmax=qmax(coef,lo,hi)
    Emax=qmax(E,lo,hi)
    Emin=qmin(E,lo,hi)
    penalty=2*theta.b*iv.sqrt(Jmax*Emax)-theta.a**2*Emin
    checks={}
    for branch,U in incumbents.items():
        q=tuple(coef[k]-U[k] for k in range(3))
        diff=qmin(q,lo,hi)
        checks[branch]=dict(diff_lower=str(diff),penalty_upper=str(penalty.b),passed=bool(diff.a>penalty.b))
    passed=any(v['passed'] for v in checks.values())
    if not passed:
        out.append(dict(rank=i,passed=False,reason='quadratic_gap',checks=checks,theta_upper=str(mp.mpf(theta.b)),cell={k:row[k] for k in ('m_lo','m_hi','h_lo','h_hi')}))
    if (i+1)%2000==0:
        print(i+1,'failures',len(out),flush=True)
result=dict(epsilon_interval=[str(lo),str(hi)],tested=limit,failures=len(out),
            first_failures=out[:100],method='fixed-center exact epsilon quadratics plus projector variation bound',
            complete_outer=limit==len(rows) and len(out)==0)
(ROOT/'stage17_epsilon_quadratic_trial.json').write_text(json.dumps(result,indent=2))
print(json.dumps(dict(epsilon_interval=result['epsilon_interval'],tested=limit,failures=len(out),first_failures=out[:3]),indent=2))
