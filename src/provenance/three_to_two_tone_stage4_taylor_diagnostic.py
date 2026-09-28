"""Vectorized float Taylor diagnostics for *every* outer-unresolved cell.

This is a prioritization tool. Only directed-interval counterparts count as
certificates; finite-difference Hessians here are not rigorous.
"""
import csv,gzip,json
from pathlib import Path
import numpy as np
from three_to_two_tone_stage3_full_crossings import crossing

ROOT=Path(__file__).resolve().parent;OUT=ROOT/'stage4_bottleneck_cells';N=21;B=10

def D(v):
    o=np.ones_like(v,dtype=float)
    for k in range(1,11):o+=2*np.cos(k*v/N)
    return o
def D1(v):
    o=np.zeros_like(v,dtype=float)
    for k in range(1,11):o-=2*k/N*np.sin(k*v/N)
    return o

def valgrad(u1,u2,e,d,E):
    h1=D(u1-d)+D(u1+d)-e*D(u1-B)
    h2=D(u2-d)+D(u2+d)-e*D(u2-B)
    hp1=D1(u1-d)+D1(u1+d)-e*D1(u1-B)
    hp2=D1(u2-d)+D1(u2+d)-e*D1(u2-B)
    g=D(u2-u1);gp=D1(u2-u1)
    den=N*N-g*g
    cap=N*(h1*h1+h2*h2)-2*g*h1*h2
    cap1=2*N*h1*hp1+2*gp*h1*h2-2*g*hp1*h2
    cap2=2*N*h2*hp2-2*gp*h1*h2-2*g*h1*hp2
    den1=2*g*gp;den2=-2*g*gp
    return E-cap/den,-(cap1*den-cap*den1)/den**2,-(cap2*den-cap*den2)/den**2

def one(z):
    p=OUT/f'z_{str(z).replace(".","p")}_unresolved.csv.gz'
    with gzip.open(p,'rt',newline='') as f:rows=list(csv.DictReader(f))
    u1=np.array([float(r['u1_center']) for r in rows]);u2=np.array([float(r['u2_center']) for r in rows])
    w1=np.array([float(r['m_hi'])-float(r['m_lo'])+float(r['h_hi'])-float(r['h_lo']) for r in rows])*N/2
    w2=w1.copy()
    q=crossing(str(z));e=float(q['epsilon_cross']);U=float(q['J_A'])
    s=np.arange(-10,11)/N
    x=2*np.cos(z*s)-e*np.exp(1j*B*s)
    E=float(np.vdot(x,x).real)
    J,g1,g2=valgrad(u1,u2,e,z,E)
    h=1e-4
    _,g1pa,g2pa=valgrad(u1+h,u2,e,z,E)
    _,g1ma,g2ma=valgrad(u1-h,u2,e,z,E)
    _,g1pc,g2pc=valgrad(u1,u2+h,e,z,E)
    _,g1mc,g2mc=valgrad(u1,u2-h,e,z,E)
    H11=(g1pa-g1ma)/(2*h);H22=(g2pc-g2mc)/(2*h)
    H12=((g1pc-g1mc)+(g2pa-g2ma))/(4*h)
    lower=J-abs(g1)*w1-abs(g2)*w2-.5*(abs(H11)*w1*w1+2*abs(H12)*w1*w2+abs(H22)*w2*w2)
    ex=lower>U+1e-8
    # Each gradient sign margin is compared with a Hessian variation proxy.
    gx1=abs(g1)-abs(H11)*w1-abs(H12)*w2
    gx2=abs(g2)-abs(H12)*w1-abs(H22)*w2
    gex=(gx1>1e-9)|(gx2>1e-9)
    out=[]
    for i,r in enumerate(rows):
        out.append([r['cell_id'],repr(float(J[i])),repr(float(g1[i])),repr(float(g2[i])),
                    repr(float(H11[i])),repr(float(H12[i])),repr(float(H22[i])),
                    repr(float(lower[i])),repr(float(lower[i]-U)),bool(ex[i]),bool(gex[i]),
                    repr(float(max(gx1[i],gx2[i])))])
    with gzip.open(OUT/f'z_{str(z).replace(".","p")}_taylor.csv.gz','wt',newline='') as f:
        w=csv.writer(f);w.writerow(['cell_id','J_center_direct','grad_u1','grad_u2',
        'H11_float','H12_float','H22_float','Taylor_lower_float',
        'Taylor_gap_float','Taylor_objective_excluded_float',
        'gradient_excluded_float','max_gradient_sign_margin_float']);w.writerows(out)
    return dict(z=z,cells=len(rows),float_taylor_objective_excluded=int(ex.sum()),
                float_gradient_excluded=int(gex.sum()),
                float_either_excluded=int((ex|gex).sum()),
                float_remaining=int((~(ex|gex)).sum()),
                min_Taylor_gap=float(np.nanmin(lower-U)),
                median_Taylor_gap=float(np.nanmedian(lower-U)),
                max_Taylor_gap=float(np.nanmax(lower-U)),
                median_old_gap=float(np.median([float(r['lower_gap_float']) for r in rows])),
                status='float_diagnostic_not_certificate')

if __name__=='__main__':
    arr=[]
    for z in (.25,.5):
        q=one(z);arr.append(q);print(q,flush=True)
    (ROOT/'stage4_taylor_diagnostic.json').write_text(json.dumps(arr,indent=2))
