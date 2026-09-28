"""110-digit baseline and branch references for interval validation."""

import json
from pathlib import Path
import mpmath as mp

mp.mp.dps=110
N=21;d=mp.mpf(2);b=mp.mpf(10)
s=[mp.mpf(t)/N for t in range(-10,11)]

def D(u):return mp.fsum(mp.cos(u*v) for v in s)

def objective(a,c,e):
    g=D(c-a)
    h1=D(a+d)+D(a-d)-e*D(a-b)
    h2=D(c+d)+D(c-d)-e*D(c-b)
    norm=mp.fsum(4*mp.cos(d*v)**2+e*e-4*e*mp.cos(d*v)*mp.cos(b*v) for v in s)
    return norm-(N*(h1*h1+h2*h2)-2*g*h1*h2)/(N*N-g*g)

def grad(a,c,e):
    return mp.diff(lambda v:objective(v,c,e),a),mp.diff(lambda v:objective(a,v,e),c)

def branch(e,guess):
    return mp.findroot(lambda a,c:grad(a,c,e),guess,tol=mp.mpf('1e-101'),
                       maxsteps=40,solver='mdnewton')

def five(a1,a2,c1,c2,e):
    return (*grad(a1,a2,e),*grad(c1,c2,e),
            objective(a1,a2,e)-objective(c1,c2,e))

sol=mp.findroot(five,(-4.77244735,.60517032,-.05765356,12.46341936,.24819063),
                tol=mp.mpf('1e-100'),maxsteps=40,solver='mdnewton')
a1,a2,c1,c2,ec=sol
eps=[('below',ec-mp.mpf('0.035')),
     ('crossing',ec),('above',ec+mp.mpf('0.035'))]
records=[]
for label,e in eps:
    aa=branch(e,(a1,a2));bb=branch(e,(c1,c2))
    ja=objective(*aa,e);jb=objective(*bb,e)
    records.append(dict(label=label,epsilon=mp.nstr(e,105),
                        A=[mp.nstr(v,105) for v in aa],
                        B=[mp.nstr(v,105) for v in bb],
                        J_A=mp.nstr(ja,105),J_B=mp.nstr(jb,105),
                        DeltaJ=mp.nstr(ja-jb,105)))
slope=mp.diff(lambda e:objective(a1,a2,e)-objective(c1,c2,e),ec)
output=dict(precision_dps=mp.mp.dps,epsilon_cross=mp.nstr(ec,105),
            A_at_crossing=[mp.nstr(a1,105),mp.nstr(a2,105)],
            B_at_crossing=[mp.nstr(c1,105),mp.nstr(c2,105)],
            DeltaJ_slope=mp.nstr(slope,105),records=records)
Path(__file__).with_name('stage2_reference.json').write_text(json.dumps(output,indent=2),encoding='utf-8')
print('ec',mp.nstr(ec,45),'slope',mp.nstr(slope,30),flush=True)
for rec in records:
    print(rec['label'],'JA',rec['J_A'][:34],'JB',rec['J_B'][:34],flush=True)
