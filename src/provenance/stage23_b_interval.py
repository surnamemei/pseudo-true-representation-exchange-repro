"""Bounded attempt at a directed continuum weak-location interval certificate.

Numerical proposals never determine acceptance. mpmath.iv encloses every
predicate. A rejected interval is retained as a failed attempt, not a theorem.
"""
import csv, json, math, time, argparse
from pathlib import Path
import numpy as np
from mpmath import iv, mp
from three_to_two_tone_stage2_interval import Jet
from stage22_b_map import solve

ROOT=Path(__file__).resolve().parent
iv.dps=65; mp.dps=75
M2=iv.mpf(1)/12; M4=iv.mpf(1)/80
def lo(x): return mp.mpf(x.a)
def hi(x): return mp.mpf(x.b)
def bounds(x): return [mp.nstr(lo(x),60),mp.nstr(hi(x),60)]
def away(x): return lo(x)>0 or hi(x)<0
def D(x,k=0):
    if max(abs(lo(x)),abs(hi(x)))<=1:
        return sum((-1)**m*x**(2*m-k)/(iv.mpf(2*m+1)*math.factorial(2*m-k)*2**(2*m))
                   for m in range((k+1)//2,38))+iv.mpf(['-1e-50','1e-50'])
    s,c=iv.sin(x/2),iv.cos(x/2); y=2*s/x
    for j in range(1,k+1): y=(2*(s,c,-s,-c)[j%4]/2**j-j*y)/x
    return y
def kd(x,k=0):
    dd=[D(x.v,k+j) for j in range(3)]
    return Jet(dd[0],tuple(dd[1]*g for g in x.g),
       tuple(tuple(dd[1]*x.H[i][j]+dd[2]*x.g[i]*x.g[j] for j in range(2)) for i in range(2)))
def tangent(v,lam,b):
    v=Jet.var(v,0); lam=Jet.var(lam,1); b=Jet.c(b)
    d,dp=kd(v),kd(v,1)
    q0=Jet.c(-M2)-lam*kd(b);q1=2*lam*kd(b,1)
    g=1-d*d-dp*dp/M2
    h=kd(v,2)-lam*kd(b-v)-d*q0+dp*q1/(2*M2)
    gain=h*h/g
    base=Jet.c(M4)-2*lam*kd(b,2)+lam*lam-q0*q0-q1*q1/(4*M2)
    return base-gain,g,h/g,gain
def system(x,b):
    a,ga,ba,pa=tangent(x[0],x[2],b);c,gc,bc,pc=tangent(x[1],x[2],b)
    f=[a.g[0],c.g[0],pc.v-pa.v]
    j=[[a.H[0][0],iv.mpf(0),a.H[0][1]],
       [iv.mpf(0),c.H[0][0],c.H[0][1]],
       [a.g[0],-c.g[0],a.g[1]-c.g[1]]]
    return f,j,(a,c,ga,gc,ba,bc)
def local(b,center,radius):
    seed=(-4.06,12.424,.066507)
    rr=solve(str(center),seed)['_mp']
    x0=[iv.mpf(mp.nstr(t,65)) for t in rr]
    # Generous frequency boxes; lambda has much less sensitivity to b.
    rad=[iv.mpf(str(radius))*6,iv.mpf(str(radius))*4,iv.mpf(str(radius))*.06]
    X=[x0[i]+iv.mpf([-rad[i].b,rad[i].b]) for i in range(3)]
    f0,j0,_=system(x0,iv.mpf(str(center)))
    cinv=np.linalg.inv(np.array([[float(v.mid) for v in row] for row in j0]))
    ci=[[iv.mpf(repr(float(v))) for v in row] for row in cinv]
    fp,_,_=system(x0,b);_,jx,z=system(X,b)
    t=[[iv.mpf(i==j)-sum(ci[i][k]*jx[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
    kk=[x0[i]-sum(ci[i][k]*fp[k] for k in range(3))+sum(t[i][j]*(X[j]-x0[j]) for j in range(3)) for i in range(3)]
    contraction=max(hi(sum(abs(t[i][j])*rad[j]/rad[i] for j in range(3))) for i in range(3))
    inside=all(lo(kk[i])>lo(X[i]) and hi(kk[i])<hi(X[i]) for i in range(3))
    a,c,ga,gc,ba,bc=z
    margins=dict(curvature=min(lo(a.H[0][0]),lo(c.H[0][0])),
       coefficient=min(lo(ba.v),-hi(bc.v)),slope=lo(a.g[1]-c.g[1]),
       gram=min(lo(ga.v),lo(gc.v)),frequency_separation=min(abs(lo(X[0])),abs(lo(X[1]))))
    out=dict(X=[bounds(x) for x in X],K=[bounds(x) for x in kk],
             center=[bounds(x) for x in x0],preconditioner=[[bounds(x) for x in row] for row in ci],
             contraction=str(contraction),inclusion=inside,margins={k:str(v) for k,v in margins.items()})
    out['pass']=bool(inside and contraction<1 and min(margins.values())>0)
    return out,X,x0
def make_eval(L,b):
    q0=-M2-L*D(b);q1=2*L*D(b,1)
    base=M4-2*L*D(b,2)+L*L-q0*q0-q1*q1/(4*M2)
    dc=[iv.mpf(0)]*49
    for k in range(0,49,2):dc[k]=iv.mpf((-1)**(k//2))/(math.factorial(k+1)*2**k)
    dp=[(k+1)*dc[k+1] for k in range(47)]
    db=[(-1)**k*D(b,k)/math.factorial(k) for k in range(47)]
    ac=[(iv.mpf(k==0))-sum(dc[j]*dc[k-j]+dp[j]*dp[k-j]/M2 for j in range(k+1)) for k in range(47)][4:45]
    bc=[(k+2)*(k+1)*dc[k+2]-L*db[k]-dc[k]*q0+dp[k]*q1/(2*M2) for k in range(47)][2:43]
    def poly(coeff,x,j):
        y=iv.mpf(0)
        for k in reversed(range(j,len(coeff))):y=y*x+coeff[k]*(math.factorial(k)//math.factorial(k-j))
        return y+iv.mpf(['-1e-25','1e-25'])
    def evaluate(a,c):
        x=iv.mpf([str(a),str(c)])
        if max(abs(a),abs(c))<=1:
            g,g1,g2=[poly(ac,x,j) for j in range(3)]
            h,h1,h2=[poly(bc,x,j) for j in range(3)]
        else:
            d,d1,d2,d3,d4=[D(x,j) for j in range(5)]
            g=1-d*d-d1*d1/M2;g1=-2*d*d1-2*d1*d2/M2
            g2=-2*(d1*d1+d*d2)-2*(d2*d2+d1*d3)/M2
            h=d2-L*D(b-x)-d*q0+d1*q1/(2*M2)
            h1=d3+L*D(b-x,1)-d1*q0+d2*q1/(2*M2)
            h2=d4-L*D(b-x,2)-d2*q0+d3*q1/(2*M2)
        if lo(g)<=0:return None
        gain=h*h/g; grad=-2*h*h1/g+h*h*g1/g**2
        curv=-2*(h1*h1+h*h2)/g+4*h*h1*g1/g**2+h*h*g2/g**2-2*h*h*g1*g1/g**3
        return base-gain,grad,curv,gain
    return evaluate,base,q0,q1
def global_check(b,X,x0,limit=16000):
    L=X[2];ev,base,q0,q1=make_eval(L,b)
    u=lo(x0[0]);trial=ev(u,u);trial_gain=trial[3]
    # Broader convex neighborhoods contain all Krawczyk frequency boxes.
    boxes=[(float(lo(x))-0.03,float(hi(x))+0.03) for x in X[:2]]
    convex=[]
    for a,c in boxes:
        z=ev(a,c);convex.append(bool(z and lo(z[2])>0))
    if not all(convex):return dict(pass_global=False,reason='convex_neighborhood_enclosure',boxes=boxes),[]
    cuts=sorted(set([-100.,100.,-1.,1.]+[v for pair in boxes for v in pair]+[k/4 for k in range(-400,401)]))
    stack=[(a,c,0) for a,c in zip(cuts,cuts[1:])]
    leaves=[];un=[];visited=0
    while stack:
        a,c,depth=stack.pop();visited+=1
        ident=next((j for j,(aa,cc) in enumerate(boxes) if a>=aa and c<=cc),None)
        if ident is not None:
            leaves.append(dict(lo=repr(a),hi=repr(c),predicate='root_'+str(ident),margin=''));continue
        z=ev(a,c)
        kind=None;margin=None
        if z:
            gap=trial_gain-z[3]
            if lo(gap)>0:kind='cost';margin=lo(gap)
            elif away(z[1]):kind='gradient';margin=min(abs(lo(z[1])),abs(hi(z[1])))
            elif hi(z[2])<0:kind='concave';margin=-hi(z[2])
        if kind:
            leaves.append(dict(lo=repr(a),hi=repr(c),predicate=kind,margin=str(margin)));continue
        if depth>=22 or visited>=limit:
            un.append([a,c]);
            if visited>=limit:un.extend([[a,c] for a,c,_ in stack]);break
        else:
            m=(a+c)/2;stack.extend([(a,m,depth+1),(m,c,depth+1)])
    q2=-M4+L*D(b,2)-M2*q0
    coal_gap=trial_gain-q2*q2/(M4-M2*M2)
    V=iv.mpf(100);dd=2/V;dp=(1+dd)/V;d2=1/(2*V)+2*dp/V
    aa=1-dd*dd-dp*dp/M2
    hh=d2+abs(L)*2/(V-b)+dd*abs(q0)+dp*abs(q1)/(2*M2)
    tail_gap=trial_gain-hh*hh/aa
    margins=[mp.mpf(r['margin']) for r in leaves if r['predicate']=='cost']
    out=dict(pass_global=not un and lo(coal_gap)>0 and lo(tail_gap)>0,
       reason='complete' if not un else 'enclosure_or_budget',visited=visited,leaves=len(leaves),
       unresolved=un,boxes=boxes,coalescent_margin=str(lo(coal_gap)),tail_margin=str(lo(tail_gap)),
       regular_cost_cell_margin=str(min(margins)) if margins else None,
       convex_neighborhoods=convex)
    return out,leaves
def run():
    t0=time.perf_counter();attempts=[]
    for width in ('0.25','0.10','0.05','0.02','0.01'):
        b=iv.mpf([str(mp.mpf(10)-mp.mpf(width)),str(mp.mpf(10)+mp.mpf(width))])
        loc,X,x0=local(b,10,width)
        rec=dict(half_width=width,b_interval=bounds(b),local=loc)
        print('width',width,'local',loc['pass'],'contraction',loc['contraction'],flush=True)
        if loc['pass']:
            glob,rows=global_check(b,X,x0);rec['global']=glob
            print('global',glob,flush=True)
            if glob['pass_global']:
                with (ROOT/'b_interval_partition.csv').open('w',newline='') as f:
                    w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
                rec['primary_verified']=True;attempts.append(rec);break
        rec['primary_verified']=False;attempts.append(rec)
    out=dict(backend='mpmath.iv/libmpi',precision=iv.dps,elapsed_seconds=time.perf_counter()-t0,
             attempts=attempts,primary_verified=any(r['primary_verified'] for r in attempts),
             independent_replay=False)
    (ROOT/'b_interval_certificate.json').write_text(json.dumps(out,indent=2))
    return out
if __name__=='__main__':run()
