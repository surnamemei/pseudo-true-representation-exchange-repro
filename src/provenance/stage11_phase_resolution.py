"""Focused numerical phase and resolution tests; no numerical global proof."""
import csv,json,sys
from pathlib import Path
import numpy as np
from scipy.optimize import root,brentq
from scipy.ndimage import minimum_filter
from mpmath import mp
ROOT=Path(__file__).resolve().parent
S=np.arange(-10,11)/21;S2=np.sum(S*S);S4=np.sum(S**4)
mp.dps=65
ss=[mp.mpf(t)/21 for t in range(-10,11)];sm2=sum(t*t for t in ss)
def dm(k,x=0):return mp.re(sum((1j*t)**k*mp.exp(1j*x*t) for t in ss))
dc=[dm(k)/mp.factorial(k) for k in range(39)]
dp=[(k+1)*dc[k+1] for k in range(38)]
aa=[];cc=[];b0=[];b1=[];hh=[]
for k in range(35):
    aa.append((21 if k==0 else 0)-sum(dc[j]*dc[k-j]/21+dp[j]*dp[k-j]/sm2 for j in range(k+1)))
    cc.append((21 if k==0 else 0)-sum(dc[j]*dc[k-j]/21 for j in range(k+1)))
    b0.append((k+2)*(k+1)*dc[k+2]+sm2*dc[k]/21)
    h=(-1)**k*dm(k,10)/mp.factorial(k)-dc[k]*dm(0,10)/21
    hh.append(h);b1.append(h-dp[k]*dm(1,10)/sm2)
ap=np.array([float(v) for v in aa[4:]]);cp=np.array([float(v) for v in cc[2:]])
b0p=np.array([float(v) for v in b0[2:]]);b1p=np.array([float(v) for v in b1[2:]])
hp=np.array([float(v) for v in hh[1:]])
def D(v,k=0):return np.sum((1j*S)**k*np.exp(1j*v*S)).real if not np.iscomplexobj(v) else np.sum(S**k*np.cos(v*S+k*np.pi/2))
def R(v,l,phi):
    c=np.cos(phi);si=np.sin(phi);d=D(v);dp=D(v,1)
    q0=-S2+l*c*D(10);q1=-2*l*c*D(10,1);q3=D(v,2)+l*c*D(10-v)
    a=21-d*d/21-dp*dp/S2;b=q3-d*q0/21+dp*q1/(2*S2)
    ci=21-d*d/21;bi=l*si*(D(10-v)-d*D(10)/21)
    if abs(v)<2:
        pol=np.polynomial.polynomial.polyval
        a=pol(v,ap);b=pol(v,b0p)+l*c*pol(v,b1p)
        ci=pol(v,cp);bi=l*si*pol(v,hp)
    return S4+2*l*c*D(10,2)+21*l*l-q0*q0/21-q1*q1/(4*S2)-(l*si*D(10))**2/21-b*b/a-bi*bi/ci
def grad(v,l,phi):return R(v+1e-6j,l,phi).imag/1e-6
def curv(v,l,phi):return (grad(v+1e-3,l,phi)-grad(v-1e-3,l,phi))/.002
def eq(x,phi):
    a,b,l=x
    return [grad(a,l,phi),grad(b,l,phi),R(a,l,phi)-R(b,l,phi)]
def write(name,rows):
    with (ROOT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
def phase():
    rows=[];seed=np.array([-4.15344048665,12.4680948443,.065980411127])
    for ph in np.linspace(np.pi,0,65):
        rr=root(eq,seed,args=(ph,),tol=1e-10)
        a,b,l=rr.x;ha=curv(a,l,ph);hb=curv(b,l,ph)
        valid=np.linalg.norm(eq(rr.x,ph))<1e-7 and l>0 and abs(a-b)>.1 and min(ha,hb)>1e-7 and min(abs(a),abs(b))>.02
        if valid:seed=rr.x
        # The deflated formula now permits the entire circle, including v=0.
        grid=np.linspace(-21*np.pi,21*np.pi,5201)
        vals=np.array([grad(v,l,ph) for v in grid]) if valid else np.array([])
        minima=[]
        if valid:
            for i in range(len(grid)-1):
                if vals[i]<0<vals[i+1]:
                    v=brentq(lambda v:grad(v,l,ph),grid[i],grid[i+1]);minima.append((R(v,l,ph),v))
        dl=1e-5
        slope=((R(a,l+dl,ph)-R(b,l+dl,ph))-(R(a,l-dl,ph)-R(b,l-dl,ph)))/(2*dl)
        rows.append(dict(phi=ph,phi_over_pi=ph/np.pi,crossing_found=valid,lambda_c=l if valid else '',v_A=a if valid else '',v_B=b if valid else '',R_A=R(a,l,ph) if valid else '',curvature_A=ha,curvature_B=hb,slope=slope,minima_count=len(minima) if valid else '',best_other_margin=min([c-R(a,l,ph) for c,v in minima if min(abs(v-a),abs(v-b))>.02],default=np.nan) if valid else '',status='numerical_crossing' if valid else 'continuation_failed_not_nonexistence'))
        print('phase',ph/np.pi,valid,a,b,l,flush=True)
    write('phase_crossing_map.csv',rows)
def resolution():
    import three_to_two_tone_branch_study as st
    from three_to_two_tone_stage3_full_crossings import crossing,value_grad
    mp.dps=60
    rows=[];seed=None
    for z in [2,2.5,np.pi,3.5,4,5,2*np.pi]:
        try:
            cr=crossing(str(z),seed=seed)
            seed=tuple(mp.mpf(cr[k]) for k in ('u_A1','u_A2','u_B1','u_B2','epsilon_cross'))
            e=float(seed[4]);ua=np.array([float(v) for v in seed[:2]]);ub=np.array([float(v) for v in seed[2:4]])
            x=st.x_record(21,z,10,e,np.pi)
            A=st.stats(ua,x,21);B=st.stats(ub,x,21)
            _,_,mins,n=st.independent_minima(x,21,grid_size=512,random_count=32)
            usp,jsp=st.refine([-z,z],x,21)
            jf=st.value([-z,z],x,21);J=float(cr['J_A'])
            valid=min(A['hess_min'],B['hess_min'])>1e-7 and np.linalg.norm(ua-ub)>.1 and e>0
            others=[r for r in mins if min(np.linalg.norm(r['u']-ua),np.linalg.norm(r['u']-ub))>.01]
            # Refine the best independent candidate with 60-digit stationary equations.
            best=mins[0];hh=mp.findroot(lambda a,b:value_grad(a,b,mp.mpf(cr['epsilon_cross']),mp.mpf(str(z)))[1:],tuple(best['u']),tol=mp.mpf('1e-48'))
            high=value_grad(*hh,mp.mpf(cr['epsilon_cross']),mp.mpf(str(z)))[0]
            rows.append(dict(z=z,separation_DFT_bins=z/np.pi,epsilon_cross=e,valid_local_exchange=valid,u_A1=ua[0],u_A2=ua[1],u_B1=ub[0],u_B2=ub[1],J_cross=J,A_hessian_min=A['hess_min'],B_hessian_min=B['hess_min'],strong_fixed_loss=jf,strong_refined_loss=jsp,strong_refined_u1=usp[0],strong_refined_u2=usp[1],full_search_best=str(high),best_u1=str(hh[0]),best_u2=str(hh[1]),other_minimum_margin=min([r['J']-J for r in others],default=np.nan),grid_candidates=n,refined_minima=len(mins),status='historical_interval_global' if z==2 else 'numerical_full_torus_512_plus_32_random'))
        except Exception as err:
            print('resolution failure',z,repr(err),flush=True);break
        print('resolution',rows[-1],flush=True)
        write('resolution_sweep.csv',rows)
if __name__=='__main__':
    (phase if len(sys.argv)<2 or sys.argv[1]=='phase' else resolution)()
