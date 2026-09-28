"""General-phase interval Krawczyk checks and numerical cusp diagnosis."""
import json,csv
from pathlib import Path
import numpy as np
from mpmath import mp,iv
from three_to_two_tone_stage2_interval import Jet
from three_to_two_tone_stage3_tangent_interval import D,D1,D2,S2,S4
from stage11_phase_resolution import dc,dp,b0,b1,aa,cc,hh
mp.dps=65;iv.dps=65
ROOT=Path(__file__).resolve().parent
def dm(v,k=0):return mp.re(sum((1j*mp.mpf(t)/21)**k*mp.exp(1j*v*mp.mpf(t)/21) for t in range(-10,11)))
sm2=-dm(0,2);sm4=dm(0,4)
def rm(v,l,ph):
    c=mp.cos(ph);ss=mp.sin(ph);d=dm(v);d1=dm(v,1)
    q0=-sm2+l*c*dm(10);q1=-2*l*c*dm(10,1)
    if abs(v)<mp.mpf('.5'):
        pol=lambda co:sum(co[k]*v**k for k in range(len(co)))
        a=pol(aa[4:]);b=pol(b0[2:])+l*c*pol(b1[2:]);aI=pol(cc[2:]);bI=l*ss*pol(hh[1:])
    else:
        a=21-d*d/21-d1*d1/sm2;b=dm(v,2)+l*c*dm(10-v)-d*q0/21+d1*q1/(2*sm2)
        aI=21-d*d/21;bI=l*ss*(dm(10-v)-d*dm(10)/21)
    return sm4+2*l*c*dm(10,2)+21*l*l-q0*q0/21-q1*q1/(4*sm2)-(l*ss*dm(10))**2/21-b*b/a-bI*bI/aI
def rj(v,l,ph):
    v=Jet.var(v,0);l=Jet.var(l,1);c=iv.cos(ph);sn=iv.sin(ph)
    d=D(v);dp=D1(v);db=D(Jet.c(10));db1=D1(Jet.c(10))
    q0=-Jet.c(S2)+l*c*db;q1=-2*l*c*db1;q3=D2(v)+l*c*D(10-v)
    a=21-d*d/21-dp*dp/S2;b=q3-d*q0/21+dp*q1/(2*S2)
    ai=21-d*d/21;bi=l*sn*(D(10-v)-d*db/21)
    return Jet.c(S4)+2*l*c*D2(Jet.c(10))+21*l*l-q0*q0/21-q1*q1/(4*S2)-(l*sn*db)*(l*sn*db)/21-b*b/a-bi*bi/ai
def fj(x,ph):
    A=rj(x[0],x[2],ph);B=rj(x[1],x[2],ph)
    return [A.g[0],B.g[0],A.v-B.v],[[A.H[0][0],iv.mpf(0),A.H[0][1]],[iv.mpf(0),B.H[0][0],B.H[0][1]],[A.g[0],-B.g[0],A.g[1]-B.g[1]]],A,B
def ib(x):return [str(mp.mpf(x.a)),str(mp.mpf(x.b))]
def main():
    rows=list(csv.DictReader(open(ROOT/'phase_crossing_map.csv')));checks=[]
    for frac in (1,.75,.5,.4):
        ph=mp.pi*mp.mpf(str(frac));r=min(rows,key=lambda r:abs(float(r['phi_over_pi'])-frac))
        x=mp.findroot(lambda a,b,l:(mp.diff(lambda v:rm(v,l,ph),a),mp.diff(lambda v:rm(v,l,ph),b),rm(a,l,ph)-rm(b,l,ph)),tuple(mp.mpf(r[k]) for k in ('v_A','v_B','lambda_c')),tol=mp.mpf('1e-52'))
        x0=[iv.mpf(str(q)) for q in x];rad=[iv.mpf('1e-12'),iv.mpf('1e-12'),iv.mpf('1e-14')]
        X=[q+iv.mpf([-r.b,r.b]) for q,r in zip(x0,rad)];phi=iv.pi*iv.mpf(str(frac))
        f0,j0,_,_=fj(x0,phi);_,jx,A,B=fj(X,phi)
        cinv=np.linalg.inv([[float(v.mid) for v in row] for row in j0]);C=[[iv.mpf(repr(float(v))) for v in row] for row in cinv]
        M=[[iv.mpf(int(i==j))-sum(C[i][k]*jx[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
        K=[x0[i]-sum(C[i][k]*f0[k] for k in range(3))+sum(M[i][j]*(X[j]-x0[j]) for j in range(3)) for i in range(3)]
        ok=all(K[i].a>X[i].a and K[i].b<X[i].b for i in range(3))
        checks.append(dict(phi_over_pi=frac,root=[str(q) for q in x],box=[ib(q) for q in X],inclusion=ok,curvature_A=ib(A.H[0][0]),curvature_B=ib(B.H[0][0]),slope=ib(A.g[1]-B.g[1]),positive=bool(A.H[0][0].a>0 and B.H[0][0].a>0 and (A.g[1]-B.g[1]).a>0)))
        print('phase interval',frac,checks[-1],flush=True)
    cusp=mp.findroot(lambda v,l,p:tuple(mp.diff(lambda u:rm(u,l,p),v,k) for k in (1,2,3)),(mp.mpf('5.1'),mp.mpf('.058'),mp.mpf('1.06')),tol=mp.mpf('1e-45'),maxsteps=60)
    co=mp.findroot(lambda b,l,p:(mp.diff(lambda v:rm(v,l,p),0),mp.diff(lambda v:rm(v,l,p),b),rm(0,l,p)-rm(b,l,p)),(mp.mpf('8.3'),mp.mpf('.067'),mp.mpf('1.46')),tol=mp.mpf('1e-40'),maxsteps=40)
    out=dict(representative_interval_checks=checks,cusp_numerical=[str(q) for q in cusp],cusp_phi_over_pi=str(cusp[2]/mp.pi),cusp_fourth_derivative=str(mp.diff(lambda v:rm(v,cusp[1],cusp[2]),cusp[0],4)),central_collision_numerical=[str(q) for q in co],central_collision_phi_over_pi=str(co[2]/mp.pi),global_phase_interval_certified=False)
    (ROOT/'stage11_phase_validation.json').write_text(json.dumps(out,indent=2));print(out,flush=True)
if __name__=='__main__':main()
