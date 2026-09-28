"""Independent symbolic Taylor/implicit-system audit of lambda_N coefficients."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import mpmath as mp
import sympy as sp

ROOT=Path(__file__).resolve().parent
mp.mp.dps=75
v,lam=sp.symbols('v lam', real=True)
ref=json.loads((ROOT/'stage14_continuum_candidate.json').read_text())
va,vb,l0=(mp.mpf(ref[k]) for k in ('v_A','v_B','lambda'))


def ds(q):
    s=2*sp.sin(q/2)/q
    return [s,q*q*s/24,7*q**4*s/5760]


def coeffs(q):
    if q == v:
        return ds(q),[sp.diff(d,q) for d in ds(q)],[sp.diff(d,q,2) for d in ds(q)]
    u=sp.symbols('u',real=True)
    return ([d.subs(u,q) for d in ds(u)],
            [sp.diff(d,u).subs(u,q) for d in ds(u)],
            [sp.diff(d,u,2).subs(u,q) for d in ds(u)])


d,dp,dpp=coeffs(v)
db,dbp,dbpp=coeffs(sp.Integer(10))
dshift=ds(10-v)
m2=[sp.Rational(1,12),-sp.Rational(1,12),sp.Integer(0)]
m4=[sp.Rational(1,80),-sp.Rational(1,24),sp.Rational(7,240)]
inv_m2=[sp.Integer(12)]*3
q0=[-m2[i]-lam*db[i] for i in range(3)]
q1=[2*lam*dbp[i] for i in range(3)]


def conv2(a,b,i):
    return sum(a[j]*b[i-j] for j in range(i+1))


def conv3(a,b,c,i):
    return sum(a[j]*b[k]*c[i-j-k] for j in range(i+1) for k in range(i-j+1))


A=[(sp.Integer(1) if i==0 else 0)-conv2(d,d,i)-conv3(inv_m2,dp,dp,i)
   for i in range(3)]
B=[dpp[i]-lam*dshift[i]-conv2(d,q0,i)+conv3(inv_m2,dp,q1,i)/2
   for i in range(3)]
C=[m4[i]-2*lam*dbpp[i]+(lam*lam if i==0 else 0)
   -conv2(q0,q0,i)-conv3(inv_m2,q1,q1,i)/4 for i in range(3)]
R0=C[0]-B[0]**2/A[0]
R1=C[1]-2*B[0]*B[1]/A[0]+B[0]**2*A[1]/A[0]**2
R2=(C[2]-(B[1]**2+2*B[0]*B[2])/A[0]
    +2*B[0]*B[1]*A[1]/A[0]**2
    +B[0]**2*(A[2]/A[0]**2-A[1]**2/A[0]**3))


def eval_expr(expr, vv):
    f=sp.lambdify((v,lam),expr,modules='mpmath',cse=True)
    return f(vv,l0)


def difference(expr):
    return eval_expr(expr,va)-eval_expr(expr,vb)


print('symbolic expressions built',flush=True)
slope=difference(sp.diff(R0,lam))
ft=difference(R1)
a=-ft/slope
hA=eval_expr(sp.diff(R0,v,2),va)
hB=eval_expr(sp.diff(R0,v,2),vb)
gA=eval_expr(sp.diff(R0,v,lam),va)*a+eval_expr(sp.diff(R1,v),va)
gB=eval_expr(sp.diff(R0,v,lam),vb)*a+eval_expr(sp.diff(R1,v),vb)
driftA=-gA/hA
driftB=-gB/hB
fll=difference(sp.diff(R0,lam,2))
flt=difference(sp.diff(R1,lam))
ftt2=difference(R2)
b=-(fll*a*a/2+flt*a+ftt2-gA*gA/(2*hA)+gB*gB/(2*hB))/slope
out={'method':'SymPy exact symbolic coefficient algebra and joint-root perturbation; mpmath evaluation',
     'centered_grid':'odd N, s_n=n/N, t=N^-2',
     'kernel_coefficients':['d(q)','q^2*d(q)/24','7*q^4*d(q)/5760'],
     'lambda_infinity':mp.nstr(l0,65),'slope':mp.nstr(slope,60),
     'a2':mp.nstr(a,60),'a4':mp.nstr(b,60),
     'vA_Nminus2':mp.nstr(driftA,50),'vB_Nminus2':mp.nstr(driftB,50),
     'curvature_A':mp.nstr(hA,50),'curvature_B':mp.nstr(hB,50),
     'components':{'F_t':mp.nstr(ft,50),'F_lambda_lambda':mp.nstr(fll,50),
                   'F_lambda_t':mp.nstr(flt,50),'F_tt_over_2':mp.nstr(ftt2,50),
                   'G_A':mp.nstr(gA,50),'G_B':mp.nstr(gB,50)}}
(ROOT/'stage16_symbolic_asymptotics.json').write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2),flush=True)

records={}
with (ROOT/'lambdaN_largeN_data.csv').open() as f:
    for row in csv.DictReader(f): records[int(row['N'])]=row
with (ROOT/'independent_N0_replay.csv').open() as f:
    pass
rows=[]
for n in (11,15,21,31,41,101,201,501,1001):
    row=records[n]
    exact=mp.mpf(row['lambda_N'])
    p0=l0;p2=l0+a/n**2;p4=p2+b/n**4
    r2=exact-p2;r4=exact-p4
    rows.append({'N':n,'origin':row['origin'],'lambda_N':mp.nstr(exact,45),
                 'lambda_infinity_only':mp.nstr(p0,45),
                 'Nminus2_prediction':mp.nstr(p2,45),
                 'Nminus4_prediction':mp.nstr(p4,45),
                 'residual_after_Nminus2':mp.nstr(r2,40),
                 'residual_after_Nminus4':mp.nstr(r4,40),
                 'N4_scaled_residual':mp.nstr(r2*n**4,35),
                 'N6_scaled_residual':mp.nstr(r4*n**6,35)})
with (ROOT/'stage16_residual_scaling.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
