"""Attempt a uniform-in-t transfer of the frozen continuum partition.

For all t in [0,1/N0^2], this encloses the exact midpoint kernel via
H(w)=sqrt(w)/sin(sqrt(w)), w=q^2*t/4. Trials do not alter Stage-14 files.
"""
from __future__ import annotations
import csv, json, math, sys
from fractions import Fraction
from pathlib import Path
import numpy as np
from mpmath import iv, mp
from three_to_two_tone_stage2_interval import Jet
import stage14_continuum_global_interval as inf
from stage17_kernel_bound import upper_error, decimal_power_ceiling

ROOT=Path(__file__).resolve().parent
iv.dps=65;mp.dps=75
N0=int(sys.argv[1]) if len(sys.argv)>1 else 1000001
assert N0>=201 and N0%2==1
T=iv.mpf([0,1])/(iv.mpf(N0)**2)
ref=json.loads((ROOT/'stage14_continuum_candidate.json').read_text())
leaves=list(csv.DictReader((ROOT/'stage14_continuum_global_partition.csv').open(encoding='utf-8')))

# Coefficients h_m of u/sin(u)=sum h_m u^(2m), from an exact recurrence.
h=[Fraction(1)]
for m in range(1,9):
    h.append(-sum(Fraction((-1)**j,math.factorial(2*j+1))*h[m-j] for j in range(1,m+1)))

kernel_majorant=upper_error(N0)
err_exponent=decimal_power_ceiling(kernel_majorant)
ERR=iv.mpf([f'-1e{err_exponent}',f'1e{err_exponent}'])

def hp(q,k):
    ans=iv.mpf(0)
    for m in range((k+1)//2,9):
        p=2*m
        fall=math.factorial(p)//math.factorial(p-k)
        frac=h[m]
        ans+=iv.mpf(frac.numerator)/frac.denominator*T**m/4**m*fall*q**(p-k)
    return ans+ERR

def D(q,k=0):
    return sum(math.comb(k,j)*inf.D(q,k-j)*hp(q,j) for j in range(k+1))

M2=(1-T)/12
M4=iv.mpf(1)/80-T/24+7*T*T/240

def jeterr():
    e=ERR
    return Jet(e,(e,e),((e,e),(e,e)))

def jh(z,k):
    ans=Jet.c(0)
    for m in range((k+1)//2,9):
        p=2*m
        fall=math.factorial(p)//math.factorial(p-k)
        frac=h[m]
        zp=Jet.c(1)
        for _ in range(p-k):zp=zp*z
        ans+=Jet.c(iv.mpf(frac.numerator)/frac.denominator*T**m/4**m*fall)*zp
    return ans+jeterr()

def JD0(z):return 2*(z/2).sin()/z
def JD1(z):return ((z/2).cos()-JD0(z))/z
def JD2(z):return (-(z/2).sin()/2-2*JD1(z))/z
def JD(z,k):
    dd=[JD0(z),JD1(z),JD2(z)]
    return sum(math.comb(k,j)*dd[k-j]*jh(z,j) for j in range(k+1))

def tangent(v,lam):
    v=Jet.var(v,0);lam=Jet.var(lam,1);weak=Jet.c(10)
    d=JD(v,0);dp=JD(v,1)
    q0=Jet.c(-M2)-lam*JD(weak,0)
    q1=2*lam*JD(weak,1)
    a=Jet.c(1)-d*d-dp*dp/M2
    b=JD(v,2)-lam*JD(weak-v,0)-d*q0+dp*q1/(2*M2)
    c=Jet.c(M4)-2*lam*JD(weak,2)+lam*lam-q0*q0-q1*q1/(4*M2)
    return c-b*b/a,a,b/a

def FJ(a,b,l):
    A=tangent(a,l)[0];B=tangent(b,l)[0]
    F=[A.g[0],B.g[0],A.v-B.v]
    J=[[A.H[0][0],iv.mpf(0),A.H[0][1]],
       [iv.mpf(0),B.H[0][0],B.H[0][1]],
       [A.g[0],-B.g[0],A.g[1]-B.g[1]]]
    return F,J,A,B

def local_check():
    x0=[iv.mpf(ref[k]) for k in ('v_A','v_B','lambda')]
    # Contains the O(t) root shifts at this threshold, with generous slack.
    tmax=mp.mpf(1)/N0**2
    rad=[max(mp.mpf('1e-8'),mp.mpf(2000)*tmax),
         max(mp.mpf('1e-8'),mp.mpf(120)*tmax),
         max(mp.mpf('1e-9'),mp.mpf(4)*tmax)]
    X=[x0[i]+iv.mpf([str(-rad[i]),str(rad[i])]) for i in range(3)]
    F0,J0,_,_=FJ(*x0)
    _,JX,A,B=FJ(*X)
    cinv=np.linalg.inv(np.array([[float(J0[i][j].mid) for j in range(3)] for i in range(3)]))
    C=[[iv.mpf(repr(float(cinv[i,j]))) for j in range(3)] for i in range(3)]
    E=[[iv.mpf(int(i==j))-sum(C[i][k]*JX[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
    Y=[X[i]-x0[i] for i in range(3)]
    K=[x0[i]-sum(C[i][k]*F0[k] for k in range(3))+sum(E[i][j]*Y[j] for j in range(3)) for i in range(3)]
    inclusion=all(K[i].a>X[i].a and K[i].b<X[i].b for i in range(3))
    contraction=max(sum(abs(E[i][j]).b for j in range(3)) for i in range(3))
    ta=tangent(X[0],X[2]);tb=tangent(X[1],X[2])
    signs=(A.H[0][0].a>0 and B.H[0][0].a>0 and
           (A.g[1]-B.g[1]).a>0 and ta[1].v.a>0 and tb[1].v.a>0 and
           ta[2].v.a>0 and tb[2].v.b<0)
    return K,{'inclusion':bool(inclusion),'contraction_upper':str(contraction),
              'signs':bool(signs),'curvature_A_lower':str(mp.mpf(A.H[0][0].a)),
              'curvature_B_lower':str(mp.mpf(B.H[0][0].a)),
              'beta_A_lower':str(mp.mpf(ta[2].v.a)),
              'beta_B_abs_lower':str(-mp.mpf(tb[2].v.b)),
              'slope_lower':str(mp.mpf((A.g[1]-B.g[1]).a)),
              'X':[[str(mp.mpf(x.a)),str(mp.mpf(x.b))] for x in X],
              'K':[[str(mp.mpf(x.a)),str(mp.mpf(x.b))] for x in K]}

def make_core(L):
    q0=-M2-L*D(iv.mpf(10))
    q1=2*L*D(iv.mpf(10),1)
    C=M4-2*L*D(iv.mpf(10),2)+L*L-q0*q0-q1*q1/(4*M2)
    K=46
    # Polynomial coefficients of d_t(q), through q^(K+2).
    dc=[iv.mpf(0)]*(K+3)
    for p in range(0,K+3,2):
        k=p//2
        dc[p]=sum(iv.mpf((-1)**(k-m))/(math.factorial(2*(k-m)+1)*2**(2*(k-m)))
                  *iv.mpf(h[m].numerator)/h[m].denominator*T**m/4**m
                  for m in range(min(k,8)+1))
    dp=[(k+1)*dc[k+1] for k in range(K+1)]
    db=[(-1)**k*D(iv.mpf(10),k)/math.factorial(k) for k in range(K+1)]
    ac=[];bc=[]
    for k in range(K+1):
        ac.append((iv.mpf(1) if k==0 else iv.mpf(0))-
                  sum(dc[j]*dc[k-j]+dp[j]*dp[k-j]/M2 for j in range(k+1)))
        bc.append((k+2)*(k+1)*dc[k+2]-L*db[k]-dc[k]*q0+dp[k]*q1/(2*M2))
    ac=ac[4:45];bc=bc[2:43]
    def poly(coeff,x,r):
        ans=iv.mpf(0)
        for k in reversed(range(r,len(coeff))):
            falling=math.prod(range(k-r+1,k+1))
            ans=ans*x+coeff[k]*falling
        return ans+iv.mpf(['-1e-25','1e-25'])
    def calc(lo,hi):
        x=iv.mpf([str(lo),str(hi)])
        if max(abs(lo),abs(hi))<=1:
            a,a1,a2=[poly(ac,x,j) for j in range(3)]
            b,b1,b2=[poly(bc,x,j) for j in range(3)]
        else:
            d,d1,d2,d3,d4=[D(x,j) for j in range(5)]
            a=1-d*d-d1*d1/M2
            a1=-2*d*d1-2*d1*d2/M2
            a2=-2*(d1*d1+d*d2)-2*(d2*d2+d1*d3)/M2
            b=d2-L*D(10-x)-d*q0+d1*q1/(2*M2)
            b1=d3+L*D(10-x,1)-d1*q0+d2*q1/(2*M2)
            b2=d4-L*D(10-x,2)-d2*q0+d3*q1/(2*M2)
        if a.a<=0:return None
        R=C-b*b/a
        g=-2*b*b1/a+b*b*a1/(a*a)
        H=-2*(b1*b1+b*b2)/a+4*b*b1*a1/(a*a)+b*b*a2/(a*a)-2*b*b*a1*a1/(a*a*a)
        return R,g,H,a
    return calc,C,q0,q1

def gap(z):return max(mp.mpf(z.a),-mp.mpf(z.b))

def main():
    X,local=local_check()
    print('local',local['inclusion'],local['signs'],flush=True)
    calc,C,q0,q1=make_core(X[2])
    U=calc(float(ref['v_A']),float(ref['v_A']))[0].b
    failed=[];failed_total=0;counts={};margins={}
    for i,r in enumerate(leaves):
        lo,hi=float(r['lo']),float(r['hi']);p=r['predicate']
        z=calc(lo,hi);ok=False
        if z is not None:
            if p=='cost_excluded':
                ok=z[0].a>U
                m=mp.mpf(z[0].a)-mp.mpf(U)
            elif p=='gradient_excluded':
                ok=inf.excludes(z[1]);m=gap(z[1])
            elif p=='monotone_derivative_no_root':
                if inf.excludes(z[2]):
                    ga,gb=calc(lo,lo)[1],calc(hi,hi)[1]
                    ok=(ga.a>0 and gb.a>0) or (ga.b<0 and gb.b<0)
                m=gap(z[2])
            elif p.startswith('unique_root_'):
                if inf.excludes(z[2]):
                    ga,gb=calc(lo,lo)[1],calc(hi,hi)[1]
                    ok=(ga.b<0 and gb.a>0) or (ga.a>0 and gb.b<0)
                m=gap(z[2])
            if p not in margins or m<mp.mpf(margins[p]):margins[p]=str(m)
        counts[p]=counts.get(p,0)+1
        if not ok:
            failed_total+=1
            if len(failed)<30:failed.append({'index':i,'cell':[lo,hi],'predicate':p})
        if (i+1)%5000==0:print('checked',i+1,'failures',len(failed),flush=True)
    V=iv.mpf(100)
    q2=-M4+X[2]*D(iv.mpf(10),2)-M2*q0
    coal=C-q2*q2/(M4-M2*M2)
    e0=2*iv.pi/V;e1=iv.pi/V;e2=iv.pi/(2*V)
    AA=1-e0*e0-e1*e1/M2
    BB=e2+X[2].b*(2*iv.pi/(V-10))+e0*abs(q0).b+e1*abs(q1).b/(2*M2.a)
    tail=C-BB*BB/AA
    root_cells=[r for r in leaves if r['predicate'].startswith('unique_root_')]
    linked=all(float(r['lo'])<float(X[0 if '_A_' in r['predicate'] else 1].a) and
               float(X[0 if '_A_' in r['predicate'] else 1].b)<float(r['hi']) for r in root_cells)
    verified=(local['inclusion'] and local['signs'] and failed_total==0 and linked and
              coal.a>U and tail.a>U)
    out={'N0_interval_trial':N0,'t_range':['0',str(mp.mpf(T.b))],
         'uniform_kernel_remainder_enclosure':f'1e{err_exponent} for derivatives through order 46 on |q|<=110; near-zero deflation 1e-25',
         'local':local,'predicate_counts':counts,'failed_total':failed_total,
         'first_failures':failed[:10],'minimum_predicate_margins':margins,
         'trial_cost_upper':str(mp.mpf(U)),
         'coalescent_margin_lower':str(mp.mpf(coal.a)-mp.mpf(U)),
         'tail_margin_lower':str(mp.mpf(tail.a)-mp.mpf(U)),
         'root_boxes_linked':bool(linked),'verified':bool(verified)}
    out['kernel_majorant']=str(mp.mpf(kernel_majorant.numerator)/kernel_majorant.denominator)
    out['kernel_error_enclosure_power']=err_exponent
    (ROOT/f'stage17_interval_trial_N{N0}.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in out.items() if k!='local'},indent=2),flush=True)

if __name__=='__main__':main()
