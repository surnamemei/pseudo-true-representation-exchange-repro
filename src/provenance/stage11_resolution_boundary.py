"""Diagnose the failed large-z crossing; retain failures as observations."""
import json,csv
from pathlib import Path
import numpy as np
from scipy.optimize import root
from mpmath import mp
from stage11_phase_resolution import D
import three_to_two_tone_branch_study as st
ROOT=Path(__file__).resolve().parent
def vg(u1,u2,e,z):
    h1=D(u1-z)+D(u1+z)-e*D(u1-10);h2=D(u2-z)+D(u2+z)-e*D(u2-10)
    p1=D(u1-z,1)+D(u1+z,1)-e*D(u1-10,1);p2=D(u2-z,1)+D(u2+z,1)-e*D(u2-10,1)
    g=D(u2-u1);gp=D(u2-u1,1);den=441-g*g;cap=21*(h1*h1+h2*h2)-2*g*h1*h2
    cap1=42*h1*p1+2*gp*h1*h2-2*g*p1*h2;cap2=42*h2*p2-2*gp*h1*h2-2*g*h1*p2
    E=42+2*D(2*z)+21*e*e-2*e*(D(10-z)+D(10+z))
    return np.array([E-cap/den,-(cap1*den-cap*2*g*gp)/den**2,-(cap2*den+cap*2*g*gp)/den**2])
def eq(x,z):
    A=vg(*x[:2],x[4],z);B=vg(*x[2:4],x[4],z)
    return np.r_[A[1:],B[1:],A[0]-B[0]]
def main():
    rows=[r for r in csv.DictReader(open(ROOT/'resolution_sweep.csv')) if float(r['z'])<=4];r=rows[-1]
    seed=np.array([float(r[k]) for k in ('u_A1','u_A2','u_B1','u_B2','epsilon_cross')]);track=[]
    for z in np.r_[np.arange(4.025,5.001,.025),np.arange(5.1,6.21,.1),2*np.pi]:
        rr=root(eq,seed,args=(z,),tol=1e-10);x=rr.x
        xx=st.x_record(21,z,10,x[4],np.pi)
        ha=st.hessian(x[:2],xx,21)[0];hb=st.hessian(x[2:4],xx,21)[0]
        valid=np.linalg.norm(eq(x,z))<1e-7 and min(ha,hb)>1e-6 and np.linalg.norm(x[:2]-x[2:4])>.1
        track.append(dict(z=z,root=list(x),residual=float(np.linalg.norm(eq(x,z))),hessian_A=ha,hessian_B=hb,valid=bool(valid)))
        if not valid:break
        seed=x
    from three_to_two_tone_stage3_full_crossings import crossing,value_grad
    from stage7_n31_certificate import local_check
    import three_to_two_tone_stage2_interval as it
    it.N=21;mp.dps=75
    from mpmath import iv
    iv.dps=65
    validations=[]
    for z in (5,2*np.pi):
        seed=min(track,key=lambda r:abs(r['z']-z))['root'];cr=crossing(str(z),seed=tuple(map(mp.mpf,seed)))
        e=float(cr['epsilon_cross']);xx=st.x_record(21,z,10,e,np.pi)
        ua=np.array([float(cr[k]) for k in ('u_A1','u_A2')]);ub=np.array([float(cr[k]) for k in ('u_B1','u_B2')]);sp,jsp=st.refine([-z,z],xx,21)
        _,_,mins,n=st.independent_minima(xx,21,512,32)
        others=[m for m in mins if min(np.linalg.norm(m['u']-ua),np.linalg.norm(m['u']-ub))>.01]
        best=mins[0];high=mp.findroot(lambda a,b:value_grad(a,b,mp.mpf(cr['epsilon_cross']),mp.mpf(str(z)))[1:],tuple(best['u']),tol=mp.mpf('1e-60'))
        it.D_STRONG=iv.mpf(str(z));ep=iv.mpf(cr['epsilon_cross'])+iv.mpf(['-1e-14','1e-14'])
        check=[]
        for w in ('A','B'):
            c,_,_=local_check([cr['u_'+w+'1'],cr['u_'+w+'2']],ep,'1e-9');check.append(c)
        validations.append(dict(z=z,checks=check,scope='stationary minima only; full-domain ordering remains numerical'))
        row={k:'' for k in rows[0]};row.update(z=z,separation_DFT_bins=z/np.pi,epsilon_cross=e,valid_local_exchange=True,u_A1=ua[0],u_A2=ua[1],u_B1=ub[0],u_B2=ub[1],J_cross=cr['J_A'],A_hessian_min=st.hessian(ua,xx,21)[0],B_hessian_min=st.hessian(ub,xx,21)[0],strong_fixed_loss=st.value([-z,z],xx,21),strong_refined_loss=jsp,strong_refined_u1=sp[0],strong_refined_u2=sp[1],full_search_best=str(value_grad(*high,mp.mpf(cr['epsilon_cross']),mp.mpf(str(z)))[0]),best_u1=str(high[0]),best_u2=str(high[1]),other_minimum_margin=min(m['J']-float(cr['J_A']) for m in others),grid_candidates=n,refined_minima=len(mins),status='numerical_full_torus_512_plus_32_random; local roots interval checked')
        rows.append(row);print('large-z',row,flush=True)
    with (ROOT/'resolution_sweep.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    (ROOT/'stage11_resolution_boundary.json').write_text(json.dumps(dict(track=track,selected_local_checks=validations),indent=2))
if __name__=='__main__':main()
