"""70-digit, non-interval bounded optimization of the hard z=.25 cells.

The optimizer uses float coordinates for its search but evaluates the exact
finite-record objective and gradient at 70 decimal digits. These results are
diagnostic and do not constitute lower-bound certificates.
"""
import csv,gzip,json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize,root as scipy_root
from mpmath import mp
from three_to_two_tone_stage3_full_crossings import crossing

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'stage4_bottleneck_cells';N=21;B=10
mp.dps=70

def D(v):return 1+sum(2*mp.cos(mp.mpf(k)*v/N) for k in range(1,11))
def D1(v):return sum(-2*mp.mpf(k)/N*mp.sin(mp.mpf(k)*v/N) for k in range(1,11))
def exact(u,e,d,E):
    a,b=map(lambda x:mp.mpf(str(x)),u)
    h1=D(a-d)+D(a+d)-e*D(a-B)
    h2=D(b-d)+D(b+d)-e*D(b-B)
    p1=D1(a-d)+D1(a+d)-e*D1(a-B)
    p2=D1(b-d)+D1(b+d)-e*D1(b-B)
    g=D(b-a);gp=D1(b-a)
    den=N*N-g*g
    cap=N*(h1*h1+h2*h2)-2*g*h1*h2
    c1=2*N*h1*p1+2*gp*h1*h2-2*g*p1*h2
    c2=2*N*h2*p2-2*gp*h1*h2-2*g*h1*p2
    g1=-(c1*den-cap*2*g*gp)/(den*den)
    g2=-(c2*den+cap*2*g*gp)/(den*den)
    return E-cap/den,(g1,g2)

def run(z=.25):
    q=crossing(str(z));e=mp.mpf(q['epsilon_cross']);d=mp.mpf(str(z))
    E=sum((2*mp.cos(d*mp.mpf(t)/N)-e*mp.cos(B*mp.mpf(t)/N))**2
          +(e*mp.sin(B*mp.mpf(t)/N))**2 for t in range(-10,11))
    U=mp.mpf(q['J_A'])
    p=OUT/f'z_{str(z).replace(".","p")}_'
    with gzip.open(str(p)+'unresolved.csv.gz','rt',newline='') as f:rows=list(csv.DictReader(f))
    with gzip.open(str(p)+'taylor.csv.gz','rt',newline='') as f:diag=list(csv.DictReader(f))
    ids=[i for i,r in enumerate(diag) if r['Taylor_objective_excluded_float']=='False'
         and r['gradient_excluded_float']=='False']
    result=[]
    for ix,i in enumerate(ids):
        r=rows[i]
        a,b,c,d0=[mp.mpf(r[k]) for k in ('m_lo','m_hi','h_lo','h_hi')]
        bounds=[(float(N*(a-d0)),float(N*(b-c))),
                (float(N*(a+c)),float(N*(b+d0)))]
        def fg(x):
            v,g=exact(x,e,d,E)
            return float(v),np.array([float(t) for t in g])
        x0=np.array([(lo+hi)/2 for lo,hi in bounds])
        rr=minimize(fg,x0,jac=True,bounds=bounds,method='L-BFGS-B',
                    options={'maxiter':200,'ftol':1e-15,'gtol':1e-12,'maxls':50})
        # A second independent bounded search guards against an optimizer
        # stopping early in the shallow along-valley direction.
        rr2=minimize(fg,x0,jac=True,bounds=bounds,method='SLSQP',
                     options={'maxiter':200,'ftol':1e-14})
        candidates=[rr.x,rr2.x,x0]
        candidates += [np.array([u,v]) for u in bounds[0] for v in bounds[1]]
        # Refine smooth edge extrema at 70 digits if a derivative changes sign.
        for axis in (0,1):
            free=1-axis
            for fixed in bounds[axis]:
                lo,hi=map(lambda x:mp.mpf(repr(x)),bounds[free])
                def derivative(v):
                    x=[None,None];x[axis]=mp.mpf(repr(fixed));x[free]=v
                    return exact(x,e,d,E)[1][free]
                gl,gh=derivative(lo),derivative(hi)
                if gl*gh<0:
                    try:
                        v=mp.findroot(derivative,(lo,hi),solver='anderson',verify=True)
                        if lo<=v<=hi:
                            x=[None,None];x[axis]=mp.mpf(repr(fixed));x[free]=v
                            candidates.append(x)
                    except (ValueError,ZeroDivisionError):pass
        best=min(((exact(x,e,d,E)[0],x) for x in candidates),key=lambda t:t[0])
        # Search the unrestricted stationary equation from every hard cell.
        st=scipy_root(lambda x:fg(x)[1],x0,method='hybr',options={'xtol':1e-11})
        stationary_inside=bool(st.success and all(bounds[j][0]<=st.x[j]<=bounds[j][1] for j in range(2))
                                  and np.linalg.norm(fg(st.x)[1])<1e-8)
        result.append(dict(cell_id=int(r['cell_id']),u1_lo=repr(bounds[0][0]),u1_hi=repr(bounds[0][1]),
            u2_lo=repr(bounds[1][0]),u2_hi=repr(bounds[1][1]),
            minimum_u1=mp.nstr(best[1][0],30),minimum_u2=mp.nstr(best[1][1],30),
            minimum_J=mp.nstr(best[0],45),minimum_gap=mp.nstr(best[0]-U,35),
            stationary_inside=stationary_inside,
            stationary_root_u1=repr(float(st.x[0])),stationary_root_u2=repr(float(st.x[1])),
            stationary_solver_success=bool(st.success),
            lbfgsb_success=bool(rr.success),slsqp_success=bool(rr2.success)))
        if ix%100==0:print(ix,'of',len(ids),flush=True)
    with gzip.open(str(p)+'optimized.csv.gz','wt',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(result[0]));w.writeheader();w.writerows(result)
    gaps=[float(r['minimum_gap']) for r in result]
    out=dict(z=z,cells=len(result),min_gap=min(gaps),median_gap=float(np.median(gaps)),
             max_gap=max(gaps),stationary_inside=sum(r['stationary_inside'] for r in result),
             lbfgsb_failures=sum(not r['lbfgsb_success'] for r in result),
             slsqp_failures=sum(not r['slsqp_success'] for r in result),
             status='non_interval_diagnostic')
    (OUT/'z_0p25_optimized_summary.json').write_text(json.dumps(out,indent=2))
    return out

if __name__=='__main__':print(run(),flush=True)
