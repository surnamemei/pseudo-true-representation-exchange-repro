"""Continuum tangent functional and high-precision candidate calculations.

This file produces candidates only. Interval verification is separate.
"""
from __future__ import annotations
import csv
import json
from pathlib import Path
from mpmath import mp
import numpy as np
from scipy.optimize import brentq

ROOT=Path(__file__).resolve().parent
mp.dps=70
MU2=mp.mpf(1)/12
MU4=mp.mpf(1)/80
MU6=mp.mpf(1)/448
V2=MU4-MU2**2
V3=MU6-MU4**2/MU2
B=mp.mpf(10)

def D(x,k=0):
    x=mp.mpf(x)
    if abs(x)<mp.mpf('1'):
        return mp.fsum((-1)**m*x**(2*m-k)/(mp.mpf(2*m+1)*mp.factorial(2*m-k)*2**(2*m))
                       for m in range((k+1)//2,55))
    vals=[2*mp.sin(x/2)/x]
    for j in range(1,k+1):
        vals.append((2*mp.sin(x/2+j*mp.pi/2)/2**j-j*vals[-1])/x)
    return vals[k]

def parts(v,l):
    d=D(v);dp=D(v,1)
    q0=-MU2-l*D(B);q1=2*l*D(B,1)
    A=1-d*d-dp*dp/MU2
    BB=D(v,2)-l*D(B-v)-d*q0+dp*q1/(2*MU2)
    C=MU4-2*l*D(B,2)+l*l-q0*q0-q1*q1/(4*MU2)
    return A,BB,C

def R(v,l):
    A,BB,C=parts(v,l)
    return C-BB*BB/A

def coal(l):
    q0=-MU2-l*D(B)
    q2=-MU4+l*D(B,2)-MU2*q0
    return parts(mp.mpf('2'),l)[2]-q2*q2/V2

def strong(l):
    return l*l*(1-D(B)**2-D(B,1)**2/MU2)

def solve():
    f1=lambda va,vb,l:mp.diff(lambda v:R(v,l),va)
    f2=lambda va,vb,l:mp.diff(lambda v:R(v,l),vb)
    f3=lambda va,vb,l:R(va,l)-R(vb,l)
    va,vb,l=mp.findroot((f1,f2,f3),(-4.07,12.38,.0665),tol=mp.mpf('1e-58'),maxsteps=50)
    A1,B1,_=parts(va,l);A2,B2,_=parts(vb,l)
    out={
      'v_A':mp.nstr(va,60),'v_B':mp.nstr(vb,60),'lambda':mp.nstr(l,60),
      'R_star':mp.nstr(R(va,l),60),
      'beta_A':mp.nstr(B1/A1,60),'beta_B':mp.nstr(B2/A2,60),
      'A_A':mp.nstr(A1,60),'A_B':mp.nstr(A2,60),
      'curvature_A':mp.nstr(mp.diff(lambda v:R(v,l),va,2),60),
      'curvature_B':mp.nstr(mp.diff(lambda v:R(v,l),vb,2),60),
      'slope':mp.nstr(mp.diff(lambda L:R(va,L)-R(vb,L),l),60),
      'coalescent':mp.nstr(coal(l),60),
      'C_base':mp.nstr(parts(va,l)[2],60),
      'strong':mp.nstr(strong(l),60)
    }
    (ROOT/'stage14_continuum_candidate.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
    return out

def landscape(l):
    # Candidate-generating scan, never a global certificate.
    vv=np.linspace(-100,100,4001)
    def grad(x):
        try:return float(mp.diff(lambda y:R(y,l),mp.mpf(x)))
        except ZeroDivisionError:return float('nan')
    gg=np.array([grad(x) if abs(x)>.1 else np.nan for x in vv])
    roots=[]
    for i in range(len(vv)-1):
        if np.isfinite(gg[i]) and np.isfinite(gg[i+1]) and gg[i]*gg[i+1]<0:
            x=brentq(grad,vv[i],vv[i+1],xtol=1e-12)
            h=mp.diff(lambda y:R(y,l),mp.mpf(x),2)
            roots.append({'v':float(x),'kind':'min' if h>0 else 'max',
                          'cost':float(R(mp.mpf(x),l)),'curvature':float(h)})
    roots.sort(key=lambda r:r['v'])
    with (ROOT/'stage14_continuum_candidate_roots.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(roots[0]));w.writeheader();w.writerows(roots)
    return roots

def main():
    out=solve()
    print(json.dumps(out,indent=2),flush=True)
    roots=landscape(mp.mpf(out['lambda']))
    mins=sorted((r for r in roots if r['kind']=='min'),key=lambda r:r['cost'])
    print('stationary',len(roots),'minima',len(mins),'top',mins[:5],flush=True)

if __name__=='__main__':main()
