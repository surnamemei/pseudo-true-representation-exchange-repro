"""Pointwise full-domain interval checks at selected continuation nodes."""
import json,csv
from pathlib import Path
import numpy as np
from mpmath import iv,mp
import three_to_two_tone_stage2_interval as si
from three_to_two_tone_stage1_global import adaptive
from three_to_two_tone_stage2_outer_global import verify_cell
from three_to_two_tone_stage2_inner_global import process
from three_to_two_tone_stage4_interval_curve import validate,root
import three_to_two_tone_branch_study as st

ROOT=Path(__file__).resolve().parent
iv.dps=65;mp.dps=75

def local_box(center,e,rad):
    x0=[iv.mpf(mp.nstr(center[i],65)) for i in range(2)]
    R=mp.mpf(rad)
    X=[iv.mpf([mp.nstr(center[i]-R,65),mp.nstr(center[i]+R,65)]) for i in range(2)]
    J0=si.objective_jet(x0[0],x0[1],e)
    JX=si.objective_jet(X[0],X[1],e)
    C0=np.linalg.inv(np.array([[float(J0.H[i][j].mid) for j in range(2)] for i in range(2)]))
    C=[[iv.mpf(repr(float(C0[i,j]))) for j in range(2)] for i in range(2)]
    M=[[iv.mpf(int(i==j))-sum(C[i][k]*JX.H[k][j] for k in range(2))
       for j in range(2)] for i in range(2)]
    Y=[X[i]-x0[i] for i in range(2)]
    K=[x0[i]-sum(C[i][k]*J0.g[k] for k in range(2))
       +sum(M[i][j]*Y[j] for j in range(2)) for i in range(2)]
    det=JX.H[0][0]*JX.H[1][1]-JX.H[0][1]*JX.H[1][0]
    return dict(inclusion=all(K[i].a>X[i].a and K[i].b<X[i].b for i in range(2)),
                positive=bool(JX.H[0][0].a>0 and det.a>0),
                det_lower=float(det.a))

def one(z):
    zstr=str(z);si.D_STRONG=iv.mpf(zstr)
    q=validate(zstr,zstr,inflate=2)
    a1,a2,b1,b2,ec=root(mp.mpf(zstr))
    e=iv.mpf(mp.nstr(ec,65))+iv.mpf(['-1e-50','1e-50'])
    A=si.objective_jet(iv.mpf(mp.nstr(a1,65)),iv.mpf(mp.nstr(a2,65)),e)
    B=si.objective_jet(iv.mpf(mp.nstr(b1,65)),iv.mpf(mp.nstr(b2,65)),e)
    U=min(A.v.b,B.v.b)
    from three_to_two_tone_stage2_interval import data_energy
    E=data_energy(e)
    (ua,_),(ub,_),x=st.two_branches(21,float(z),10,float(ec),np.pi,
                                     seed_a=(float(a1),float(a2)),seed_b=(float(b1),float(b2)))
    complete,summary,leaves,_=adaptive(x,(ua,ub),float(U)+1e-6,21,
                                        radius_u=2.,maxdepth=38,write_rows=True)
    bad=[];margins=[]
    for cell in leaves:
        if cell[-1]!='excluded':continue
        ok,bb,why=verify_cell(cell,e,E,U,strong_d=iv.mpf(zstr))
        if not ok:bad.append(why)
        elif bb is not None:margins.append(float(bb.a-U))
    iv.dps=45
    inner=[]
    for lab,center in [('A',(a1,a2)),('B',(b1,b2))]:
        done,visited,rows,remain=process('crossing',lab,
                      [mp.nstr(ec-mp.mpf('1e-50'),70),mp.nstr(ec+mp.mpf('1e-50'),70)],
                      tuple(map(float,center)),U,max_cells=30000,
                      known_center_override=[mp.nstr(v,65) for v in center],
                      local_radius_override='9e-6')
        loc=local_box(center,e,'1e-5')
        inner.append(dict(branch=lab,complete=done,visited=visited,
                          remaining=remain,local=loc))
    ok=bool(q['verified'] and complete and not bad and
            all(v['complete'] and v['local']['inclusion'] and v['local']['positive'] for v in inner))
    return dict(z=zstr,epsilon_cross=mp.nstr(ec,60),crossing_root_verified=q['verified'],
                outer_proposal_complete=complete,outer_cells=len(leaves),
                outer_failures=len(bad),outer_min_margin=min(margins) if margins else None,
                inner=inner,pointwise_global_verified=ok)

if __name__=='__main__':
    results=[]
    for z in (1.0,1.5):
        q=one(z);results.append(q)
        print(z,q['pointwise_global_verified'],q['outer_cells'],q['outer_failures'],
              [(r['branch'],r['complete'],r['visited']) for r in q['inner']],flush=True)
    (ROOT/'stage4_pointwise_global_checks.json').write_text(json.dumps(results,indent=2))
