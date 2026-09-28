"""High-precision finite-N tangent problem for the small-spacing theorem."""
import csv,json
from pathlib import Path
from mpmath import mp

ROOT=Path(__file__).resolve().parent
mp.dps=85

def tangent(v,lam,N=21,b=10,phi=None):
    if phi is None:phi=mp.pi
    s=[mp.mpf(t)/N for t in range(-(N-1)//2,(N+1)//2)]
    z=[-u*u+lam*mp.e**(1j*(b*u+phi)) for u in s]
    cols=[[mp.mpc(1) for u in s],[1j for u in s],[2j*u for u in s],
          [mp.e**(1j*v*u) for u in s],[1j*mp.e**(1j*v*u) for u in s]]
    M=mp.matrix([[mp.re(sum(mp.conj(cols[i][k])*cols[j][k] for k in range(N)))
                  for j in range(5)] for i in range(5)])
    q=mp.matrix([mp.re(sum(mp.conj(cols[i][k])*z[k] for k in range(N))) for i in range(5)])
    c=mp.lu_solve(M,q)
    model=[sum(c[i]*cols[i][k] for i in range(5)) for k in range(N)]
    residual=[z[k]-model[k] for k in range(N)]
    cost=mp.re(sum(mp.conj(w)*w for w in residual))
    beta=mp.mpc(c[3],c[4])
    dv=-2*mp.re(sum(mp.conj(residual[k])*beta*1j*s[k]*cols[3][k] for k in range(N)))
    dl=2*mp.re(sum(mp.conj(residual[k])*mp.e**(1j*(b*s[k]+phi)) for k in range(N)))
    return cost,dv,dl,[c[i] for i in range(5)],M

def solve(N=21,b=10,phi=None,guess=(-4.15,12.47,.066)):
    if phi is None:phi=mp.pi
    def eq(a,z,l):
        A=tangent(a,l,N,b,phi);B=tangent(z,l,N,b,phi)
        return A[1],B[1],A[0]-B[0]
    va,vb,lam=mp.findroot(eq,tuple(mp.mpf(z) for z in guess),
                          tol=mp.mpf('1e-75'),maxsteps=35,solver='mdnewton')
    A=tangent(va,lam,N,b,phi);B=tangent(vb,lam,N,b,phi)
    haa=mp.diff(lambda u:tangent(u,lam,N,b,phi)[1],va)
    hbb=mp.diff(lambda u:tangent(u,lam,N,b,phi)[1],vb)
    trans=A[2]-B[2]
    def cubic_cost(v,coeff):
        alpha=mp.mpc(coeff[0],coeff[1]);kappa=coeff[2]
        beta=mp.mpc(coeff[3],coeff[4])
        s=[mp.mpf(t)/N for t in range(-(N-1)//2,(N+1)//2)]
        r0=[-u*u+lam*mp.e**(1j*(b*u+phi))-alpha-2j*kappa*u
            -beta*mp.e**(1j*v*u) for u in s]
        r1=[u**4/12-1j*kappa*u*alpha+kappa*kappa*u*u for u in s]
        return 2*mp.re(sum(mp.conj(r0[k])*r1[k] for k in range(N)))
    TA=cubic_cost(va,A[3]);TB=cubic_cost(vb,B[3])
    lambda1=-(TA-TB)/trans
    return dict(N=N,b=b,phi=mp.nstr(phi,80),lambda_N=mp.nstr(lam,80),
                v_A=mp.nstr(va,80),v_B=mp.nstr(vb,80),
                R_A=mp.nstr(A[0],80),R_B=mp.nstr(B[0],80),
                Rvv_A=mp.nstr(haa,80),Rvv_B=mp.nstr(hbb,80),
                delta_R_lambda=mp.nstr(trans,80),
                T_A=mp.nstr(TA,80),T_B=mp.nstr(TB,80),
                lambda_quartic=mp.nstr(lambda1,80),
                A_coeff=[mp.nstr(z,80) for z in A[3]],
                B_coeff=[mp.nstr(z,80) for z in B[3]],
                A_Gram_det=mp.nstr(mp.det(A[4]),80),
                B_Gram_det=mp.nstr(mp.det(B[4]),80))

if __name__=='__main__':
    result=solve()
    (ROOT/'stage3_tangent_reference.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))
