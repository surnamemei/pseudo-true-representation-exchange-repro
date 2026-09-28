"""Parametric interval Krawczyk validation of the 5-equation crossing curve.

This certifies only two local stationary crossing branches. Full-domain
globality is a separate obligation.
"""
import csv,json
from pathlib import Path
import numpy as np
from mpmath import iv,mp
from three_to_two_tone_stage2_interval import Jet
from three_to_two_tone_stage3_full_crossings import crossing

ROOT=Path(__file__).resolve().parent
N=21;B=10
iv.dps=55;mp.dps=75
cache={}

def D(x):
    y=Jet.c(1)
    for k in range(1,11):y+=2*(x*iv.mpf(k)/N).cos()
    return y
def D1(x):
    y=Jet.c(0)
    for k in range(1,11):y+=(x*iv.mpf(k)/N).sin()*(-2*iv.mpf(k)/N)
    return y
def D2(x):
    y=Jet.c(0)
    for k in range(1,11):y+=(x*iv.mpf(k)/N).cos()*(-2*(iv.mpf(k)/N)**2)
    return y

def objective(a,c,e,d):
    a=Jet.var(a,0);c=Jet.var(c,1)
    e=iv.mpf(e);d=iv.mpf(d)
    h1=D(a+d)+D(a-d)-D(a-B)*e
    h2=D(c+d)+D(c-d)-D(c-B)*e
    he1=-D(a-B);he2=-D(c-B)
    hd1=-D1(a-d)+D1(a+d);hd2=-D1(c-d)+D1(c+d)
    hdd1=D2(a-d)+D2(a+d);hdd2=D2(c-d)+D2(c+d)
    g=D(c-a);den=Jet.c(N*N)-g*g
    cap=N*(h1*h1+h2*h2)-2*g*h1*h2
    cape=2*N*(h1*he1+h2*he2)-2*g*(he1*h2+h1*he2)
    capd=2*N*(h1*hd1+h2*hd2)-2*g*(hd1*h2+h1*hd2)
    capdd=2*N*(hd1*hd1+h1*hdd1+hd2*hd2+h2*hdd2)\
        -2*g*(hdd1*h2+2*hd1*hd2+h1*hdd2)
    E=iv.mpf(0);Ee=iv.mpf(0);Ed=iv.mpf(0);Edd=iv.mpf(0)
    for t in range(-10,11):
        s=iv.mpf(t)/N;co=iv.cos(d*s);cb=iv.cos(B*s)
        E+=4*co*co+e*e-4*e*co*cb
        Ee+=2*e-4*co*cb
        Ed+=-8*s*co*iv.sin(d*s)+4*e*s*iv.sin(d*s)*cb
        Edd+=-8*s*s*(co*co-iv.sin(d*s)**2)+4*e*s*s*co*cb
    return (Jet.c(E)-cap/den,Jet.c(Ee)-cape/den,
            Jet.c(Ed)-capd/den,Jet.c(Edd)-capdd/den,den.v)

def system(x,d):
    a,ae,ad,add,ag=objective(x[0],x[1],x[4],d)
    b,be,bd,bdd,bg=objective(x[2],x[3],x[4],d)
    F=[a.g[0],a.g[1],b.g[0],b.g[1],a.v-b.v]
    Z=iv.mpf(0)
    J=[[a.H[0][0],a.H[0][1],Z,Z,ae.g[0]],
       [a.H[1][0],a.H[1][1],Z,Z,ae.g[1]],
       [Z,Z,b.H[0][0],b.H[0][1],be.g[0]],
       [Z,Z,b.H[1][0],b.H[1][1],be.g[1]],
       [a.g[0],a.g[1],-b.g[0],-b.g[1],ae.v-be.v]]
    Fd=[ad.g[0],ad.g[1],bd.g[0],bd.g[1],ad.v-bd.v]
    Fdd=[add.g[0],add.g[1],bdd.g[0],bdd.g[1],add.v-bdd.v]
    return F,J,Fd,Fdd,a,b,ag,bg

def root(d):
    key=str(d)
    if key not in cache:
        q=crossing(key)
        cache[key]=[mp.mpf(q[k]) for k in ('u_A1','u_A2','u_B1','u_B2','epsilon_cross')]
    return cache[key]

def bounds(z):return float(z.a),float(z.b)
def lo(z):return float(z.a)
def hi(z):return float(z.b)

def validate(left,right,inflate=1.4,base_radius='1e-11'):
    left=mp.mpf(str(left));right=mp.mpf(str(right));mid=(left+right)/2
    xl=root(left);xr=root(right);xm=root(mid)
    r=[max(abs(xl[k]-xm[k]),abs(xr[k]-xm[k]))*mp.mpf(str(inflate))
       +mp.mpf(base_radius) for k in range(5)]
    x0=[iv.mpf(mp.nstr(z,65)) for z in xm]
    X=[iv.mpf([mp.nstr(xm[k]-r[k],65),mp.nstr(xm[k]+r[k],65)]) for k in range(5)]
    z=iv.mpf([mp.nstr(left,40),mp.nstr(right,40)])
    Fmid,J0,Fdmid,_,_,_,_,_=system(x0,iv.mpf(mp.nstr(mid,60)))
    _,_,_,Fdd,_,_,_,_=system(x0,z)
    dz=z-iv.mpf(mp.nstr(mid,60))
    F0=[Fmid[i]+Fdmid[i]*dz+iv.mpf('0.5')*Fdd[i]*dz*dz for i in range(5)]
    FX,JX,_,_,A,BB,ga,gb=system(X,z)
    M0=np.array([[float(J0[i][j].mid) for j in range(5)] for i in range(5)])
    try:C0=np.linalg.inv(M0)
    except np.linalg.LinAlgError:return dict(z_lo=str(left),z_hi=str(right),included=False,reason='singular_midpoint')
    C=[[iv.mpf(repr(float(C0[i,j]))) for j in range(5)] for i in range(5)]
    G=[[iv.mpf(int(i==j))-sum(C[i][k]*JX[k][j] for k in range(5))
        for j in range(5)] for i in range(5)]
    Y=[X[k]-x0[k] for k in range(5)]
    # Preserve the common parameter increment before interval evaluation;
    # preconditioning five independent F_j(D) intervals loses cancellations.
    K=[x0[i]-sum(C[i][j]*Fmid[j] for j in range(5))
       -sum(C[i][j]*Fdmid[j] for j in range(5))*dz
       -iv.mpf('0.5')*sum(C[i][j]*Fdd[j] for j in range(5))*dz*dz
       +sum(G[i][j]*Y[j] for j in range(5)) for i in range(5)]
    included=all(K[i].a>X[i].a and K[i].b<X[i].b for i in range(5))
    contraction=max(sum(abs(G[i][j]).b for j in range(5)) for i in range(5))
    ah=A.H;bh=BB.H
    adet=ah[0][0]*ah[1][1]-ah[0][1]*ah[1][0]
    bdet=bh[0][0]*bh[1][1]-bh[0][1]*bh[1][0]
    hpos=ah[0][0].a>0 and adet.a>0 and bh[0][0].a>0 and bdet.a>0
    slope=(A.v-BB.v) # cost not slope; use ae-be at K for slope below
    _,_,_,_,AK,BK,_,_=system(K,z)
    _,AE,_,_,_=objective(K[0],K[1],K[4],z)
    _,BE,_,_,_=objective(K[2],K[3],K[4],z)
    crossing_slope=AE.v-BE.v
    sepA=K[1]-K[0];sepB=K[3]-K[2]
    ok=bool(included and contraction<1 and hpos and crossing_slope.a>0
            and sepA.a>0 and sepB.a>0 and ga.a>0 and gb.a>0)
    return dict(z_lo=str(left),z_hi=str(right),included=bool(included),
                contraction_upper=hi(contraction),A_H11_lower=lo(ah[0][0]),
                A_Hdet_lower=lo(adet),B_H11_lower=lo(bh[0][0]),
                B_Hdet_lower=lo(bdet),slope_lower=lo(crossing_slope),
                A_sep_lower=lo(sepA),B_sep_lower=lo(sepB),
                A_gram_lower=lo(ga),B_gram_lower=lo(gb),
                verified=ok,radii=[mp.nstr(v,12) for v in r],
                K_inside=bool(included),
                K_radius_ratio=[float((max(abs(K[i].a-x0[i].mid),abs(K[i].b-x0[i].mid))/iv.mpf(mp.nstr(r[i],60))).b) for i in range(5)],
                K_width=[float((K[i].b-K[i].a).b) for i in range(5)],
                X_bounds=[[str(X[i].a),str(X[i].b)] for i in range(5)])

if __name__=='__main__':
    rows=[]
    for k in range(10):
        a=mp.mpf('1.9999')+mp.mpf(k)*mp.mpf('0.00001')
        b=a+mp.mpf('0.00001')
        q=validate(a,b,inflate=2)
        rows.append(q)
        print(q['z_lo'],q['z_hi'],q.get('verified'),q.get('contraction_upper'),flush=True)
    (ROOT/'stage4_interval_boxes.json').write_text(json.dumps(rows,indent=2))
