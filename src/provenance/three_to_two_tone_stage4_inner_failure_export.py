"""Detailed audit of the 81 A/B-neighborhood parameter-transfer failures."""
import csv,json,math
from pathlib import Path
import numpy as np
from scipy.optimize import minimize,root as scipy_root
from mpmath import iv,mp
import three_to_two_tone_stage2_interval as si
from three_to_two_tone_stage4_taylor_diagnostic import valgrad
from three_to_two_tone_stage4_cell_optimization import exact
from three_to_two_tone_stage4_interval_curve import root
from three_to_two_tone_stage4_targeted_transfer import e0,Utransfer

ROOT=Path(__file__).resolve().parent
iv.dps=50;mp.dps=70

def endpoint(x):return float(str(x).strip('[]').split(',')[0])

def main():
    fails=json.loads((ROOT/'stage4_targeted_transfer_inner.json').read_text())['failures']
    src=[r for r in csv.DictReader((ROOT/'stage2_inner_cells.csv').open()) if r['label']=='crossing']
    lookup={(r['branch'],r['status'],r['depth'],r['u1_lo'],r['u2_lo']):r for r in src}
    z=mp.mpf('1.99995');ec=root(z)[4]
    s=[mp.mpf(t)/21 for t in range(-10,11)]
    E=sum((2*mp.cos(z*t)-ec*mp.cos(10*t))**2+(ec*mp.sin(10*t))**2 for t in s)
    ea=iv.mpf([str(e0-mp.mpf('5e-5')),str(e0+mp.mpf('5e-5'))])
    si.D_STRONG=iv.mpf(['1.9999','2'])
    target=float(Utransfer)
    out=[]
    for i,f in enumerate(fails):
        r=lookup[(f['branch'],f['type'],f['depth'],f['u1_lo'],f['u2_lo'])]
        lo1,hi1,lo2,hi2=map(float,(r['u1_lo'],r['u1_hi'],r['u2_lo'],r['u2_hi']))
        x0=np.array([(lo1+hi1)/2,(lo2+hi2)/2])
        X=[iv.mpf([r['u1_lo'],r['u1_hi']]),iv.mpf([r['u2_lo'],r['u2_hi']])]
        mid=[iv.mpf(repr(float(v))) for v in x0]
        J0=si.objective_jet(mid[0],mid[1],ea)
        JX=si.objective_jet(X[0],X[1],ea)
        Y=[X[k]-mid[k] for k in range(2)]
        L=J0.v+sum(J0.g[k]*Y[k] for k in range(2))
        L+=sum(iv.mpf('.5')*JX.H[k][j]*Y[k]*Y[j] for k in range(2) for j in range(2))
        sep=x0[1]-x0[0]
        gram=441-(1+sum(2*math.cos(k*sep/21) for k in range(1,11)))**2
        gg=math.sqrt(max(0.,441-gram));cond=(21+gg)/(21-gg)
        def fg(u):
            a,b,c=valgrad(np.array(u[0]),np.array(u[1]),float(ec),float(z),float(E))
            return float(a),np.array([float(b),float(c)])
        rr=minimize(fg,x0,jac=True,bounds=[(lo1,hi1),(lo2,hi2)],method='L-BFGS-B',
                    options={'ftol':1e-15,'gtol':1e-11,'maxiter':300})
        candidates=[rr.x,x0]+[np.array([x,y]) for x in (lo1,hi1) for y in (lo2,hi2)]
        best=min(((exact(q,ec,z,E)[0],q) for q in candidates),key=lambda t:t[0])
        sr=scipy_root(lambda u:fg(u)[1],x0,method='hybr',options={'xtol':1e-11})
        inside=bool(sr.success and lo1<=sr.x[0]<=hi1 and lo2<=sr.x[1]<=hi2
                    and np.linalg.norm(fg(sr.x)[1])<1e-7)
        out.append(dict(id=i,branch=f['branch'],old_status=f['type'],depth=r['depth'],
            u1_lo=r['u1_lo'],u1_hi=r['u1_hi'],u2_lo=r['u2_lo'],u2_hi=r['u2_hi'],
            u1_center=x0[0],u2_center=x0[1],width_u1=hi1-lo1,width_u2=hi2-lo2,
            gram_det_center=gram,gram_condition_center=cond,
            coalescent_sep_u=sep,distance_to_sep_0p1_strip_u=max(0,sep-.1),
            archived_objective_lower=r['J_lower'],
            current_interval_Taylor_lower=str(L.a),
            gap_current_lower_to_uniform_candidate=endpoint(L.a)-target,
            center_gradient_u1=str(J0.g[0]),center_gradient_u2=str(J0.g[1]),
            box_min_cost_70digit=mp.nstr(best[0],40),
            box_min_gap_to_pointwise_AB=mp.nstr(best[0]-min(
                exact(root(z)[:2],ec,z,E)[0],exact(root(z)[2:4],ec,z,E)[0]),30),
            stationary_solver_success=bool(sr.success),
            stationary_inside_numeric=inside,
            stationary_root_u1=float(sr.x[0]),stationary_root_u2=float(sr.x[1]),
            diagnosis='parameter_interval_enclosure_or_moving_branch; not_certified'))
    with (ROOT/'stage4_inner_failure_cells.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(out[0]));w.writeheader();w.writerows(out)
    print('cells',len(out),'stationary_inside',sum(t['stationary_inside_numeric'] for t in out),
          'min_numeric_gap',min(float(t['box_min_gap_to_pointwise_AB']) for t in out))

if __name__=='__main__':main()
