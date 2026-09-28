"""Directed interval isolation of the complete N=21 tangent circle.
No frozen certificate is modified. Run: python stage11_tangent_global.py
"""
import csv,json,math
from pathlib import Path
import numpy as np
from mpmath import iv,mp
from three_to_two_tone_stage3_tangent_interval import prove

iv.dps=60;mp.dps=70
ROOT=Path(__file__).resolve().parent
ref=json.loads((ROOT/'stage3_tangent_reference.json').read_text())
local=prove('1e-24','1e-25')
assert local['verified']
L=iv.mpf(ref['lambda_N'])+iv.mpf(['-1e-25','1e-25'])
N=21;s=[iv.mpf(t)/N for t in range(-10,11)]
S2=sum(t*t for t in s);S4=sum(t**4 for t in s)
def dn(x,k=0):
    if k%2==0:return sum((-1)**(k//2)*t**k*iv.cos(x*t) for t in s)
    return sum((-1)**((k+1)//2)*t**k*iv.sin(x*t) for t in s)
q0=-S2-L*dn(10);q1=2*L*dn(10,1)
C=S4-2*L*dn(10,2)+N*L*L-q0*q0/N-q1*q1/(4*S2)
# Deflate the exact v^4 and v^2 zeros before interval evaluation near 0.
K=46
dc=[dn(iv.mpf(0),k)/math.factorial(k) for k in range(K+3)]
db=[(-1)**k*dn(iv.mpf(10),k)/math.factorial(k) for k in range(K+1)]
dp=[(k+1)*dc[k+1] for k in range(K+1)]
ac=[];bc=[]
for k in range(K+1):
    ac.append((N if k==0 else 0)-sum(dc[j]*dc[k-j]/N+dp[j]*dp[k-j]/S2 for j in range(k+1)))
    bc.append((k+2)*(k+1)*dc[k+2]-L*db[k]-dc[k]*q0/N+dp[k]*q1/(2*S2))
ac=ac[4:45];bc=bc[2:43]
def poly(c,x,r):
    out=iv.mpf(0)
    for k in reversed(range(r,len(c))):out=out*x+c[k]*math.factorial(k)/math.factorial(k-r)
    # For |v|<=1 and r<=2, exponential coefficient majorants for A,B
    # give omitted tails < 1e4/(40!) < 1e-43; 1e-35 is conservative.
    return out+iv.mpf(['-1e-35','1e-35'])
def calc(lo,hi):
    x=iv.mpf([str(lo),str(hi)])
    if max(abs(lo),abs(hi))<=1:
        a,a1,a2=[poly(ac,x,k) for k in range(3)]
        b,b1,b2=[poly(bc,x,k) for k in range(3)]
    else:
        d,d1,d2,d3=[dn(x,k) for k in range(4)]
        a=N-d*d/N-d1*d1/S2
        a1=-2*d*d1/N-2*d1*d2/S2
        a2=-2*(d1*d1+d*d2)/N-2*(d2*d2+d1*d3)/S2
        b=d2-L*dn(10-x)-d*q0/N+d1*q1/(2*S2)
        b1=d3+L*dn(10-x,1)-d1*q0/N+d2*q1/(2*S2)
        b2=dn(x,4)-L*dn(10-x,2)-d2*q0/N+d3*q1/(2*S2)
    if a.a<=0:return None
    R=C-b*b/a
    g=-2*b*b1/a+b*b*a1/(a*a)
    H=-2*(b1*b1+b*b2)/a+4*b*b1*a1/(a*a)+b*b*a2/(a*a)-2*b*b*a1*a1/(a*a*a)
    return R,g,H,a
def excludes(x):return x.a>0 or x.b<0
def bounds(x):return [str(mp.mpf(x.a)),str(mp.mpf(x.b))]

def main():
    stack=[(k/4,(k+1)/4) for k in range(-264,264)]
    roots=[];leaves=[];un=[];count=0
    while stack:
        a,b=stack.pop(); count+=1;v=calc(a,b)
        if v is not None and excludes(v[1]):
            leaves.append([a,b,'derivative_excluded']);continue
        if v is not None and excludes(v[2]):
            ga=calc(a,a)[1];gb=calc(b,b)[1]
            if (ga.b<0 and gb.a>0) or (ga.a>0 and gb.b<0):
                lo,hi=a,b
                for _ in range(38):
                    mid=(lo+hi)/2;gm=calc(mid,mid)[1]
                    if not excludes(gm):break
                    if (gm.a>0)==(ga.a>0):lo=mid
                    else:hi=mid
                rv=calc(lo,hi)
                roots.append(dict(lo=lo,hi=hi,v=(lo+hi)/2,kind='min' if v[2].a>0 else 'max',R=bounds(rv[0]),curvature=bounds(rv[2])))
                leaves.append([a,b,'unique_root']);continue
        if b-a<1e-12:un.append([a,b]);continue
        mid=(a+b)/2;stack.extend([(a,mid),(mid,b)])
        if count%5000==0:print('tangent cells',count,'pending',len(stack),flush=True)
    roots.sort(key=lambda r:r['v'])
    mins=[r for r in roots if r['kind']=='min']
    selected=sorted(mins,key=lambda r:float(r['R'][0]))[:2]
    other=[r for r in mins if r not in selected]
    # All central-coalescent O(p) perturbations lie in complex{1,s}+real{s^2}.
    qnorm=S4-2*L*dn(10,2)+N*L*L
    even2=S4-S2*S2/N
    q2= -S4+L*dn(10,2)-S2*q0/N
    confl=qnorm-q0*q0/N-q1*q1/(4*S2)-q2*q2/even2
    fixed=L*L*(N-dn(10)**2/N-dn(10,1)**2/S2)
    costupper=max(mp.mpf(r['R'][1]) for r in selected)
    strict=all(mp.mpf(r['R'][0])>costupper for r in other) and mp.mpf(confl.a)>costupper
    result=dict(N=21,b=10,phase='pi',precision_dps=iv.dps,lambda_box=bounds(L),local_crossing=local,
      domain_cover=[-66,66],period='2*pi*21',visited=count,unresolved=un,stationary_points=roots,
      coalescent_coefficient=bounds(confl),fixed_strong_pair_optimized_amplitudes=bounds(fixed),
      fixed_strong_pair_unit_amplitudes=bounds(N*L*L),weak_location_R=bounds(calc(10,10)[0]),
      other_minimum_margin=str(min(mp.mpf(r['R'][0]) for r in other)-costupper),
      coalescent_margin=str(mp.mpf(confl.a)-costupper),
      tangent_classification_certified=not un and strict,
      global_smallz_theorem_requires='compactification and tangent-cone argument in smallz_global_tangent_classification.md',
      finite_positive_z_radius_certified=False)
    (ROOT/'smallz_global_certificate.json').write_text(json.dumps(result,indent=2))
    with (ROOT/'stage11_tangent_partition.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['lo','hi','predicate']);w.writerows(leaves)
    with (ROOT/'stage11_tangent_roots.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=roots[0]);w.writeheader();w.writerows(roots)
    print(json.dumps({k:v for k,v in result.items() if k not in ('local_crossing','stationary_points')},indent=2),flush=True)
if __name__=='__main__':main()
