"""Directed-interval centered Taylor bound on the 433 hard z=.25 cells.

The quadratic is minimized exactly over each rectangle at 70-digit precision;
interval coefficient radii and interval third derivatives protect that result.
"""
import csv,gzip,json
from pathlib import Path
from mpmath import iv,mp
import three_to_two_tone_stage2_interval as si
from three_to_two_tone_stage4_interval_taylor_cells import setup
from three_to_two_tone_stage4_jet3 import objective_coeffs,max_third_remainder

ROOT=Path(__file__).resolve().parent;OUT=ROOT/'stage4_bottleneck_cells';N=21
iv.dps=45;mp.dps=75

def end(x):return mp.mpf(str(x).strip('[]').split(',')[0])
def coeff(q,k):
    v=q.get(k);lo=end(v.a);hi=end(v.b)
    return (lo+hi)/2,(hi-lo)/2

def quadratic_min(c,w1,w2):
    c00,c10,c01,c20,c11,c02=[c[k] for k in ((0,0),(1,0),(0,1),(2,0),(1,1),(0,2))]
    def f(x,y):return c00+c10*x+c01*y+c20*x*x+c11*x*y+c02*y*y
    xs=(-w1,w1);ys=(-w2,w2)
    cand=[f(x,y) for x in xs for y in ys]
    for x in xs:
        if c02!=0:
            y=-(c01+c11*x)/(2*c02)
            if -w2<=y<=w2:cand.append(f(x,y))
    for y in ys:
        if c20!=0:
            x=-(c10+c11*y)/(2*c20)
            if -w1<=x<=w1:cand.append(f(x,y))
    den=4*c20*c02-c11*c11
    if den!=0:
        x=(c11*c01-2*c02*c10)/den
        y=(c11*c10-2*c20*c01)/den
        if -w1<=x<=w1 and -w2<=y<=w2:cand.append(f(x,y))
    return min(cand)

def quadratic_lower_directed(q,w1,w2):
    """Outward lower bound for the box minimum of the quadratic Taylor part.

    Every box minimum occurs at a corner, an edge critical point, or an
    interior critical point. Including critical points outside the box only
    lowers the result, so no delicate interval feasibility test is needed.
    """
    a,b,c,d,e,f=[q.get(k) for k in ((0,0),(1,0),(0,1),(2,0),(1,1),(0,2))]
    W1=iv.mpf(mp.nstr(w1,75));W2=iv.mpf(mp.nstr(w2,75))
    def val(x,y):return a+b*x+c*y+d*x*x+e*x*y+f*y*y
    cand=[val(x,y) for x in (-W1,W1) for y in (-W2,W2)]
    if not (f.a<=0<=f.b):
        for x in (-W1,W1):
            base=a+b*x+d*x*x;sl=c+e*x
            ystar=-sl/(2*f)
            if not (ystar.b<-W2.b or ystar.a>W2.b):
                cand.append(base-sl*sl/(4*f))
    if not (d.a<=0<=d.b):
        for y in (-W2,W2):
            base=a+c*y+f*y*y;sl=b+e*y
            xstar=-sl/(2*d)
            if not (xstar.b<-W1.b or xstar.a>W1.b):
                cand.append(base-sl*sl/(4*d))
    det=4*d*f-e*e
    if not (det.a<=0<=det.b):
        xstar=(e*c-2*f*b)/det
        ystar=(e*b-2*d*c)/det
        if not (xstar.b<-W1.b or xstar.a>W1.b or
                ystar.b<-W2.b or ystar.a>W2.b):
            cand.append(a-(f*b*b-e*b*c+d*c*c)/det)
    return min(end(v.a) for v in cand)

def one():
    z=.25;si.D_STRONG=iv.mpf('.25');e,U,root_ok=setup(z)
    # setup changes precision internally; restore a practical directed precision.
    iv.dps=45
    p=OUT/'z_0p25_'
    with gzip.open(str(p)+'unresolved.csv.gz','rt',newline='') as f:rows=list(csv.DictReader(f))
    with gzip.open(str(p)+'taylor.csv.gz','rt',newline='') as f:diag=list(csv.DictReader(f))
    ids=[i for i,r in enumerate(diag) if r['Taylor_objective_excluded_float']=='False'
         and r['gradient_excluded_float']=='False']
    res=[]
    for ii,i in enumerate(ids):
        r=rows[i];a,b,c,d=[mp.mpf(r[k]) for k in ('m_lo','m_hi','h_lo','h_hi')]
        u1lo=N*(a-d);u1hi=N*(b-c);u2lo=N*(a+c);u2hi=N*(b+d)
        m1=(u1lo+u1hi)/2;m2=(u2lo+u2hi)/2
        w1=(u1hi-u1lo)/2;w2=(u2hi-u2lo)/2
        X1=iv.mpf([mp.nstr(u1lo,75),mp.nstr(u1hi,75)])
        X2=iv.mpf([mp.nstr(u2lo,75),mp.nstr(u2hi,75)])
        q0=objective_coeffs(iv.mpf(mp.nstr(m1,75)),iv.mpf(mp.nstr(m2,75)),e,iv.mpf('.25'))
        qX=objective_coeffs(X1,X2,e,iv.mpf('.25'))
        ks=((0,0),(1,0),(0,1),(2,0),(1,1),(0,2))
        cr={k:coeff(q0,k) for k in ks}
        c0={k:cr[k][0] for k in ks}
        qmin=quadratic_min(c0,w1,w2)
        qlower=quadratic_lower_directed(q0,w1,w2)
        coeff_err=sum(cr[k][1]*w1**k[0]*w2**k[1] for k in ks)
        rem=max_third_remainder(qX,iv.mpf(mp.nstr(w1,75)),iv.mpf(mp.nstr(w2,75)))
        bound=qlower-end(rem.b)
        gap=bound-end(U.b)
        res.append(dict(cell_id=int(r['cell_id']),quadratic_min=mp.nstr(qmin,35),
                        quadratic_directed_lower=mp.nstr(qlower,35),
                        coefficient_error=mp.nstr(coeff_err,20),
                        third_remainder=mp.nstr(end(rem.b),20),
                        rigorous_lower=mp.nstr(bound,35),
                        rigorous_gap=mp.nstr(gap,30),excluded=bool(gap>0)))
        if ii%50==0:print(ii,'of',len(ids),'excluded',sum(t['excluded'] for t in res),flush=True)
    with gzip.open(str(p)+'centered_interval.csv.gz','wt',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(res[0]));w.writeheader();w.writerows(res)
    out=dict(z=z,cells=len(res),excluded=sum(t['excluded'] for t in res),
             remaining=sum(not t['excluded'] for t in res),
             min_gap=min(float(t['rigorous_gap']) for t in res),
             max_remainder=max(float(t['third_remainder']) for t in res),
             crossing_root_verified=root_ok,
             status='directed_interval_certificate_for_selected_cells_only')
    (OUT/'z_0p25_centered_interval_summary.json').write_text(json.dumps(out,indent=2))
    return out

if __name__=='__main__':print(one(),flush=True)
