"""Complete continuum tangent root isolation on [-100,100] plus an analytic tail bound."""
from __future__ import annotations
import csv
import json
import math
from pathlib import Path
from mpmath import iv,mp

ROOT=Path(__file__).resolve().parent
iv.dps=65;mp.dps=75
loc=json.loads((ROOT/'stage14_continuum_local_interval.json').read_text(encoding='utf-8'))
L=iv.mpf(loc['K'][2])
M2=iv.mpf(1)/12;M4=iv.mpf(1)/80;M6=iv.mpf(1)/448

def D(x,k=0):
    if max(abs(mp.mpf(x.a)),abs(mp.mpf(x.b)))<=1:
        out=iv.mpf(0)
        for m in range((k+1)//2,38):
            out+=(-1)**m*x**(2*m-k)/(iv.mpf(2*m+1)*math.factorial(2*m-k)*2**(2*m))
        # all derivative tails for k<=5 on |x|<=1 are below 1e-55
        return out+iv.mpf(['-1e-50','1e-50'])
    out=2*iv.sin(x/2)/x
    for j in range(1,k+1):
        out=(2*iv.sin(x/2+j*iv.pi/2)/2**j-j*out)/x
    return out

q0=-M2-L*D(iv.mpf(10))
q1=2*L*D(iv.mpf(10),1)
C=M4-2*L*D(iv.mpf(10),2)+L*L-q0*q0-q1*q1/(4*M2)

K=46
dc=[iv.mpf(0)]*(K+3)
for m in range((K+3)//2):
    k=2*m
    if k<len(dc):dc[k]=iv.mpf((-1)**m)/(math.factorial(k+1)*2**k)
dp=[(k+1)*dc[k+1] for k in range(K+1)]
db=[(-1)**k*D(iv.mpf(10),k)/math.factorial(k) for k in range(K+1)]
ac=[];bc=[]
for k in range(K+1):
    ac.append((iv.mpf(1) if k==0 else iv.mpf(0))-sum(dc[j]*dc[k-j]+dp[j]*dp[k-j]/M2 for j in range(k+1)))
    bc.append((k+2)*(k+1)*dc[k+2]-L*db[k]-dc[k]*q0+dp[k]*q1/(2*M2))
ac=ac[4:45];bc=bc[2:43]

def poly(coeff,x,r):
    out=iv.mpf(0)
    for k in reversed(range(r,len(coeff))):
        out=out*x+coeff[k]*math.factorial(k)/math.factorial(k-r)
    return out+iv.mpf(['-1e-30','1e-30'])

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

def excludes(z):return z.a>0 or z.b<0
def bounds(z):return [str(mp.mpf(z.a)),str(mp.mpf(z.b))]

def main():
    stack=[(k/4,(k+1)/4) for k in range(-400,400)]
    roots=[];leaves=[];un=[];vis=0
    while stack:
        a,b=stack.pop();vis+=1
        z=calc(a,b)
        if z is not None and excludes(z[1]):
            leaves.append((a,b,'gradient_excluded'));continue
        if z is not None and excludes(z[2]):
            ga=calc(a,a)[1];gb=calc(b,b)[1]
            if (ga.b<0 and gb.a>0) or (ga.a>0 and gb.b<0):
                lo,hi=a,b
                for _ in range(31):
                    mid=(lo+hi)/2
                    gm=calc(mid,mid)[1]
                    if not excludes(gm):break
                    if (gm.a>0)==(ga.a>0):lo=mid
                    else:hi=mid
                rv=calc(lo,hi)
                roots.append({'lo':lo,'hi':hi,'v':(lo+hi)/2,'kind':'min' if z[2].a>0 else 'max',
                              'R_lo':str(mp.mpf(rv[0].a)),'R_hi':str(mp.mpf(rv[0].b)),
                              'H_lo':str(mp.mpf(rv[2].a)),'H_hi':str(mp.mpf(rv[2].b))})
                leaves.append((a,b,'unique_root'));continue
        if b-a<1e-10:
            un.append([a,b]);continue
        m=(a+b)/2;stack.extend([(a,m),(m,b)])
        if vis%5000==0:print('visited',vis,'pending',len(stack),flush=True)
    roots.sort(key=lambda z:z['v'])
    mins=[r for r in roots if r['kind']=='min']
    cross_upper=max(mp.mpf(x) for r in mins for x in [r['R_hi']]
                    if (abs(r['v']+4.06)<.5 or abs(r['v']-12.42)<.5))
    other=[r for r in mins if abs(r['v']+4.06)>=.5 and abs(r['v']-12.42)>=.5]
    other_gap=min(mp.mpf(r['R_lo'])-cross_upper for r in other)
    # Coalescent endpoint:
    V2=M4-M2*M2
    q2=-M4+L*D(iv.mpf(10),2)-M2*q0
    coal=C-q2*q2/V2
    # For |v|>=100, D,D',D'' and D(10-v) are bounded by 2/V,
    # (1+2/V)/V, 1/(2V)+2(1+2/V)/V^2, and 2/(V-10), respectively.
    V=iv.mpf(100)
    dd=iv.mpf(2)/V
    dpb=(1+dd)/V
    d2b=1/(2*V)+2*dpb/V
    dbw=2/(V-10)
    AA=1-dd*dd-dpb*dpb/M2
    BB=d2b+L.b*dbw+dd*abs(q0).b+dpb*abs(q1).b/(2*M2.a)
    tail=C-BB*BB/AA
    out={
      'domain':[-100,100],'precision_dps':iv.dps,'lambda_box':bounds(L),
      'visited':vis,'unresolved':un,'stationary_count':len(roots),'minimum_count':len(mins),
      'crossing_upper_from_scalar_cover':str(cross_upper),
      'other_regular_margin_lower':str(other_gap),
      'coalescent_interval':bounds(coal),
      'coalescent_margin_lower':str(mp.mpf(coal.a)-cross_upper),
      'tail_lower_interval':bounds(tail),
      'tail_margin_lower':str(mp.mpf(tail.a)-cross_upper),
      'tail_bound_method':'summation-by-parts/sinc derivative envelope at |v|>=100',
      'verified_global':bool(not un and len(mins)>=2 and other_gap>0 and coal.a>cross_upper and tail.a>cross_upper)
    }
    (ROOT/'stage14_continuum_global_interval.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
    with (ROOT/'stage14_continuum_global_roots.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(roots[0]));w.writeheader();w.writerows(roots)
    with (ROOT/'stage14_continuum_global_partition.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f);w.writerow(['lo','hi','predicate']);w.writerows(leaves)
    print(json.dumps(out,indent=2),flush=True)

if __name__=='__main__':main()
