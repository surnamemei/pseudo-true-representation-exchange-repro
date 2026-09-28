"""Directed-interval Krawczyk validation of the continuum A/B crossing."""
import json
from pathlib import Path
import numpy as np
from mpmath import iv, mp
from three_to_two_tone_stage2_interval import Jet

ROOT=Path(__file__).resolve().parent
iv.dps=85;mp.dps=95
ref=json.loads((ROOT/'stage14_continuum_candidate.json').read_text(encoding='utf-8'))
M2=iv.mpf(1)/12;M4=iv.mpf(1)/80; B=Jet.c(10)

def D0(z):return 2*(z/2).sin()/z
def D1(z):return ((z/2).cos()-D0(z))/z
def D2(z):return (-(z/2).sin()/2-2*D1(z))/z

def tangent(v,lam):
    v=Jet.var(v,0);lam=Jet.var(lam,1)
    d=D0(v);dp=D1(v)
    q0=Jet.c(-M2)-lam*D0(B)
    q1=2*lam*D1(B)
    gram=Jet.c(1)-d*d-dp*dp/M2
    mixed=D2(v)-lam*D0(B-v)-d*q0+dp*q1/(2*M2)
    base=Jet.c(M4)-2*lam*D2(B)+lam*lam-q0*q0-q1*q1/(4*M2)
    R=base-mixed*mixed/gram
    return R,gram,mixed/gram

def FJ(a,b,l):
    A=tangent(a,l)[0];BB=tangent(b,l)[0]
    F=[A.g[0],BB.g[0],A.v-BB.v]
    J=[[A.H[0][0],iv.mpf(0),A.H[0][1]],
       [iv.mpf(0),BB.H[0][0],BB.H[0][1]],
       [A.g[0],-BB.g[0],A.g[1]-BB.g[1]]]
    return F,J,A,BB

def bounds(z):return [str(mp.mpf(z.a)),str(mp.mpf(z.b))]

def main():
    x0=[iv.mpf(ref[k]) for k in ('v_A','v_B','lambda')]
    rad=[iv.mpf('1e-15'),iv.mpf('1e-15'),iv.mpf('1e-16')]
    X=[x0[i]+iv.mpf([-rad[i].b,rad[i].b]) for i in range(3)]
    F0,J0,_,_=FJ(*x0)
    FX,JX,A,BB=FJ(*X)
    mid=np.array([[float(J0[i][j].mid) for j in range(3)] for i in range(3)])
    cinv=np.linalg.inv(mid)
    C=[[iv.mpf(repr(float(cinv[i,j]))) for j in range(3)] for i in range(3)]
    T=[[iv.mpf(int(i==j))-sum(C[i][k]*JX[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
    Y=[X[i]-x0[i] for i in range(3)]
    K=[x0[i]-sum(C[i][k]*F0[k] for k in range(3))+sum(T[i][j]*Y[j] for j in range(3)) for i in range(3)]
    included=all(K[i].a>X[i].a and K[i].b<X[i].b for i in range(3))
    contraction=max(sum(abs(T[i][j]).b for j in range(3)) for i in range(3))
    ta=tangent(X[0],X[2]);tb=tangent(X[1],X[2])
    out={
      'precision_dps':iv.dps,
      'X':[bounds(z) for z in X],'K':[bounds(z) for z in K],
      'krawczyk_inclusion':bool(included),
      'contraction_upper':str(contraction),
      'Rvv_A':bounds(A.H[0][0]),'Rvv_B':bounds(BB.H[0][0]),
      'crossing_slope':bounds(A.g[1]-BB.g[1]),
      'Gram_A':bounds(ta[1].v),'Gram_B':bounds(tb[1].v),
      'beta_A':bounds(ta[2].v),'beta_B':bounds(tb[2].v),
      'verified_local':bool(included and contraction<1 and A.H[0][0].a>0 and BB.H[0][0].a>0
                            and (A.g[1]-BB.g[1]).a>0 and ta[1].v.a>0 and tb[1].v.a>0
                            and ta[2].v.a>0 and tb[2].v.b<0)
    }
    (ROOT/'stage14_continuum_local_interval.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
