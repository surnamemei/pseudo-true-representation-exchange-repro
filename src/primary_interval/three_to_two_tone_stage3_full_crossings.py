"""Multiprecision finite-record branch crossings for theory comparison.

These roots are numerical unless separately covered by the Stage-2 interval
global certificate. No globality is inferred from this continuation.
"""
import csv,json
from pathlib import Path
from mpmath import mp

ROOT=Path(__file__).resolve().parent
ref=json.loads((ROOT/'stage3_tangent_reference.json').read_text())
baseline=json.loads((ROOT/'stage2_reference.json').read_text())
mp.dps=75

def D(v,N):
    return 1+sum(2*mp.cos(mp.mpf(k)*v/N) for k in range(1,(N-1)//2+1))
def D1(v,N):
    return sum(-2*mp.mpf(k)/N*mp.sin(mp.mpf(k)*v/N) for k in range(1,(N-1)//2+1))

def value_grad(u1,u2,e,d,N=21,b=10):
    x=[2*mp.cos(d*mp.mpf(t)/N)-e*mp.e**(1j*b*mp.mpf(t)/N)
       for t in range(-(N-1)//2,(N+1)//2)]
    E=mp.re(sum(mp.conj(v)*v for v in x))
    h1=D(u1-d,N)+D(u1+d,N)-e*D(u1-b,N)
    h2=D(u2-d,N)+D(u2+d,N)-e*D(u2-b,N)
    hp1=D1(u1-d,N)+D1(u1+d,N)-e*D1(u1-b,N)
    hp2=D1(u2-d,N)+D1(u2+d,N)-e*D1(u2-b,N)
    g=D(u2-u1,N);gp=D1(u2-u1,N)
    den=N*N-g*g
    cap=N*(h1*h1+h2*h2)-2*g*h1*h2
    cap1=2*N*h1*hp1+2*gp*h1*h2-2*g*hp1*h2
    cap2=2*N*h2*hp2-2*gp*h1*h2-2*g*h1*hp2
    den1=2*g*gp;den2=-2*g*gp
    J=E-cap/den
    return J,-(cap1*den-cap*den1)/den**2,-(cap2*den-cap*den2)/den**2

def crossing(d,N=21,b=10,seed=None):
    d=mp.mpf(str(d));p=d*d
    if seed is None:
        a=(mp.mpf(ref['v_A']),mp.mpf(ref['A_coeff'][2])*p)
        z=(mp.mpf(ref['B_coeff'][2])*p,mp.mpf(ref['v_B']))
        e=mp.mpf(ref['lambda_N'])*p
        seed=(*a,*z,e)
    def eq(a1,a2,b1,b2,eps):
        A=value_grad(a1,a2,eps,d,N,b);B=value_grad(b1,b2,eps,d,N,b)
        return A[1],A[2],B[1],B[2],A[0]-B[0]
    sol=mp.findroot(eq,seed,tol=mp.mpf('1e-62'),maxsteps=60,solver='mdnewton')
    a1,a2,b1,b2,e=sol
    A=value_grad(a1,a2,e,d,N,b);B=value_grad(b1,b2,e,d,N,b)
    lam=mp.mpf(ref['lambda_N'])
    lam4=mp.mpf(ref['lambda_quartic'])
    pred4=lam*p+lam4*p*p
    return dict(N=N,d_NDelta=mp.nstr(d,30),epsilon_cross=mp.nstr(e,70),
                epsilon_leading=mp.nstr(lam*p,70),
                epsilon_through_quartic=mp.nstr(pred4,70),
                remainder=mp.nstr(e-lam*p,70),
                remainder_over_d4=mp.nstr((e-lam*p)/p**2,50),
                relative_error=mp.nstr(abs(e-lam*p)/e,30),
                remainder_after_quartic=mp.nstr(e-pred4,60),
                remainder_after_quartic_over_d6=mp.nstr((e-pred4)/p**3,40),
                relative_error_quartic=mp.nstr(abs(e-pred4)/e,30),
                u_A1=mp.nstr(a1,60),u_A2=mp.nstr(a2,60),
                u_B1=mp.nstr(b1,60),u_B2=mp.nstr(b2,60),
                J_A=mp.nstr(A[0],60),J_B=mp.nstr(B[0],60),
                Rscaled_A=mp.nstr(A[0]/p**2,40),Rscaled_B=mp.nstr(B[0]/p**2,40),
                globality=('interval_certified' if d==2 and N==21 else 'not_certified'),
                numerical_precision_dps=mp.dps)

if __name__=='__main__':
    rows=[]
    for d in ('0.15','0.2','0.3','0.4','0.5','0.8','1.0','1.5','2.0'):
        z=crossing(d)
        rows.append(z)
        print(d,z['epsilon_cross'][:28],z['remainder_over_d4'][:18],flush=True)
    with (ROOT/'lambdaN_theory_vs_certified.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
