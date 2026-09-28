"""Interval Krawczyk proof for the N=21 symmetric tangent crossing."""
import json
from pathlib import Path
import numpy as np
from mpmath import iv
from three_to_two_tone_stage2_interval import Jet

ROOT=Path(__file__).resolve().parent
iv.dps=90
ref=json.loads((ROOT/'stage3_tangent_reference.json').read_text())
N=21;B=10
S2=iv.mpf(N*N-1)/(12*N)
S4=iv.mpf((N*N-1)*(3*N*N-7))/(240*N**3)

def D(z):
    out=Jet.c(1)
    for k in range(1,11):out+=2*(z*iv.mpf(k)/N).cos()
    return out
def D1(z):
    out=Jet.c(0)
    for k in range(1,11):out+=(z*iv.mpf(k)/N).sin()*(-2*iv.mpf(k)/N)
    return out
def D2(z):
    out=Jet.c(0)
    for k in range(1,11):out+=(z*iv.mpf(k)/N).cos()*(-2*(iv.mpf(k)/N)**2)
    return out

def tangent_jet(v,lam):
    v=Jet.var(v,0);lam=Jet.var(lam,1)
    d=D(v);dp=D1(v)
    q0=-(Jet.c(S2)+lam*D(Jet.c(B)))
    q1=2*lam*D1(Jet.c(B))
    q3=D2(v)-lam*D(Jet.c(B)-v)
    gram=Jet.c(N)-d*d/N-dp*dp/S2
    mixed=q3-d*q0/N+dp*q1/(2*S2)
    energy=Jet.c(S4)-2*lam*D2(Jet.c(B))+lam*lam*N
    R=energy-q0*q0/N-q1*q1/(4*S2)-mixed*mixed/gram
    beta=mixed/gram
    alpha=(q0-d*beta)/N
    kappa=(q1+2*dp*beta)/(4*S2)
    return R,gram,beta,alpha,kappa

def FJ(a,b,l):
    A=tangent_jet(a,l)[0];B=tangent_jet(b,l)[0]
    F=[A.g[0],B.g[0],A.v-B.v]
    J=[[A.H[0][0],iv.mpf(0),A.H[0][1]],
       [iv.mpf(0),B.H[0][0],B.H[0][1]],
       [A.g[0],-B.g[0],A.g[1]-B.g[1]]]
    return F,J,A,B

def sup_abs(z):return abs(z).b
def interval(z):return f'[{z.a},{z.b}]'

def cubic_envelope(v,lam):
    _,_,beta,alpha,kappa=tangent_jet(v,lam)
    beta=beta.v;alpha=alpha.v;kappa=kappa.v
    acc=iv.mpf(0)
    for t in range(-10,11):
        s=iv.mpf(t)/N
        rre=-s*s-lam*iv.cos(B*s)-alpha-beta*iv.cos(v*s)
        rim=-lam*iv.sin(B*s)-2*kappa*s-beta*iv.sin(v*s)
        wre=s**4/12+kappa*kappa*s*s
        wim=-kappa*alpha*s
        acc+=2*(rre*wre+rim*wim)
    return acc

def prove(radius_v='1e-9',radius_lambda='1e-10'):
    x0=[iv.mpf(ref['v_A']),iv.mpf(ref['v_B']),iv.mpf(ref['lambda_N'])]
    rad=[iv.mpf(radius_v),iv.mpf(radius_v),iv.mpf(radius_lambda)]
    X=[x0[i]+iv.mpf([-rad[i].b,rad[i].b]) for i in range(3)]
    F0,J0,_,_=FJ(*x0)
    FX,JX,A,B=FJ(*X)
    point=np.array([[float(J0[i][j].mid) for j in range(3)] for i in range(3)])
    cinv=np.linalg.inv(point)
    C=[[iv.mpf(repr(float(cinv[i,j]))) for j in range(3)] for i in range(3)]
    M=[[iv.mpf(int(i==j))-sum(C[i][k]*JX[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
    Y=[X[i]-x0[i] for i in range(3)]
    K=[x0[i]-sum(C[i][k]*F0[k] for k in range(3))+sum(M[i][j]*Y[j] for j in range(3)) for i in range(3)]
    included=all(K[i].a>X[i].a and K[i].b<X[i].b for i in range(3))
    contraction=max(sum(sup_abs(M[i][j]) for j in range(3)) for i in range(3))
    ta=tangent_jet(X[0],X[2]);tb=tangent_jet(X[1],X[2])
    TK_A=cubic_envelope(K[0],K[2]);TK_B=cubic_envelope(K[1],K[2])
    _,_,AK,BK=FJ(*K)
    lambda4=-(TK_A-TK_B)/(AK.g[1]-BK.g[1])
    result=dict(precision_dps=iv.dps,N=N,weak_b=10,
                X=[interval(z) for z in X],K=[interval(z) for z in K],
                krawczyk_inclusion=bool(included),contraction_upper=interval(contraction),
                contraction_below_one=bool(contraction<1),
                Rvv_A=interval(A.H[0][0]),Rvv_B=interval(B.H[0][0]),
                delta_R_lambda=interval(A.g[1]-B.g[1]),
                A_gram=interval(ta[1].v),B_gram=interval(tb[1].v),
                beta_A=interval(ta[2].v),beta_B=interval(tb[2].v),
                kappa_A=interval(ta[4].v),kappa_B=interval(tb[4].v),
                T_A=interval(TK_A),T_B=interval(TK_B),
                lambda_quartic=interval(lambda4),
                verified=bool(included and contraction<1 and A.H[0][0].a>0
                              and B.H[0][0].a>0 and (A.g[1]-B.g[1]).a>0
                              and ta[1].v.a>0 and tb[1].v.a>0
                              and ta[2].v.a>0 and tb[2].v.b<0))
    return result

if __name__=='__main__':
    z=prove()
    (ROOT/'stage3_tangent_interval.json').write_text(json.dumps(z,indent=2))
    print(json.dumps(z,indent=2))
