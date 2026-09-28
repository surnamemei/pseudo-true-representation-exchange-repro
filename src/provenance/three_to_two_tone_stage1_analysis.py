"""Stage-1 branch continuation, scaling, regime map, and figures."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

import three_to_two_tone_branch_study as st
import three_to_two_tone_explore as ex
from three_to_two_tone_stage1_asymptotic import tangent_crossing
from three_to_two_tone_stage1_global import evaluate

ROOT=Path(__file__).resolve().parent
FIG=ROOT/'three_to_two_tone_stage1_figures';FIG.mkdir(exist_ok=True)


def save(name,rows):
    with (ROOT/name).open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def fig(name):
    plt.tight_layout();plt.savefig(FIG/name,dpi=180);plt.close()


def amp_cross(ratio,n=21,d=2,b=10,phase=np.pi):
    def diff(eps):
        st.set_n(n)
        x=ex.atom(-d/n)+ratio*ex.atom(d/n)+eps*np.exp(1j*phase)*ex.atom(b/n)
        ua,ja=st.refine([-d,d],x,n)
        ub,jb=st.refine([0,b],x,n)
        return ja-jb if np.linalg.norm(ua-ub)>.2 else np.nan
    ee=np.linspace(.05,.6,24);vv=np.array([diff(v) for v in ee])
    for i in range(len(ee)-1):
        if np.isfinite(vv[i]) and np.isfinite(vv[i+1]) and vv[i]*vv[i+1]<0:
            return brentq(diff,ee[i],ee[i+1],xtol=1e-10)
    return np.nan


def main():
    n,d,b,phase=21,2.,10.,np.pi
    ec,ua0,ub0,*_=st.crossing(n,d,b,phase)
    # Exact branch continuation, with projector fit reoptimized at every epsilon.
    continuation=[]
    ee=np.linspace(ec-.075,ec+.075,151)
    for ep in ee:
        (ua,ja),(ub,jb),x=st.two_branches(n,d,b,ep,phase,ua0,ub0)
        sa,sb=st.stats(ua,x,n),st.stats(ub,x,n)
        continuation.append(dict(epsilon=ep,J_A=ja,J_B=jb,DeltaJ=ja-jb,
                                 u1_A=ua[0],u2_A=ua[1],u1_B=ub[0],u2_B=ub[1],
                                 hess_min_A=sa['hess_min'],hess_min_B=sb['hess_min'],
                                 sep_A=sa['sep_u'],sep_B=sb['sep_u']))
    save('branch_continuation.csv',continuation)
    # Finite-N tangent-space coefficient lambda_N and direct crossing sweep.
    nn_values=(17,21,25,31,41,61)
    ds=(.3,.4,.5,.6,.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4)
    tangent={nn:tangent_crossing(nn,b,phase)[0] for nn in nn_values}
    scaling=[]
    for nn in nn_values:
        for dd in ds:
            z=st.crossing(nn,dd,b,phase,lo=.025*dd*dd,hi=.115*dd*dd)
            if z is None:continue
            ep,ua,ub,ja,jb,sa,sb=z
            scaling.append(dict(N=nn,d_NDelta=dd,b_Nomega3=b,phase_rad=phase,
                                epsilon_cross=ep,lambda_numeric=ep/dd**2,
                                lambda_tangent=tangent[nn],
                                epsilon_tangent=tangent[nn]*dd**2,
                                tangent_relative_error=(tangent[nn]*dd**2-ep)/ep,
                                u1_A=ua[0],u2_A=ua[1],u1_B=ub[0],u2_B=ub[1],
                                min_hess=min(sa['hess_min'],sb['hess_min']),
                                min_branch_sep=min(sa['sep_u'],sb['sep_u'])))
            print('SCALE',nn,dd,ep,flush=True)
    # Fit p on genuinely small d, and finite-d correction c_N (empirical).
    for nn in nn_values:
        rr=[r for r in scaling if r['N']==nn]
        small=[r for r in rr if r['d_NDelta']<=1.]
        p,inter=np.polyfit(np.log([r['d_NDelta'] for r in small]),
                           np.log([r['epsilon_cross'] for r in small]),1)
        xx=np.array([r['d_NDelta']**2 for r in small])
        yy=np.array([r['lambda_numeric']-tangent[nn] for r in small])
        correction=(xx@yy)/(xx@xx)
        for r in rr:
            r['small_d_loglog_p']=p
            r['small_d_loglog_prefactor']=np.exp(inter)
            r['empirical_quartic_coefficient']=correction
            r['epsilon_tangent_plus_empirical_quartic']=tangent[nn]*r['d_NDelta']**2+correction*r['d_NDelta']**4
    save('scaling_fit.csv',scaling)
    print('SCALING COUNTS',len(scaling),flush=True)
    # Location-phase map: this is a test of the A/B seed family, not a proof
    # that no other crossing exists when the test returns not_found.
    maprows=[]
    bvals=np.arange(5.,16.,1.)
    phvals=np.linspace(0,2*np.pi,9)
    for ph in phvals:
        for bb in bvals:
            z=st.crossing(n,d,float(bb),float(ph),lo=.025,hi=.8)
            if z is None:
                status='not_found';ep=np.nan;sep=np.nan;hmin=np.nan;gap=np.nan
            else:
                ep,ua,ub,ja,jb,sa,sb=z
                sep=min(sa['sep_u'],sb['sep_u'])
                hmin=min(sa['hess_min'],sb['hess_min'])
                gap=float(np.linalg.norm(ua-ub))
                status='candidate' if (sep>1 and hmin>1e-5 and gap>1) else 'degenerate_or_unstable'
            maprows.append(dict(N=n,d_NDelta=d,b_Nomega3=float(bb),
                                r_location=(float(bb)+d)/(2*d),phase_rad=float(ph),
                                strong_amp_ratio=1.,status=status,epsilon_cross=ep,
                                min_branch_sep=sep,min_hess=hmin,branch_distance=gap,
                                independent_global_check=0))
        print('MAP PHASE',ph,flush=True)
    for ratio in (.9,.95,1.,1.05,1.1):
        ep=amp_cross(ratio)
        maprows.append(dict(N=n,d_NDelta=d,b_Nomega3=b,r_location=(b+d)/(2*d),
                            phase_rad=phase,strong_amp_ratio=ratio,
                            status='candidate' if np.isfinite(ep) else 'not_found',
                            epsilon_cross=ep,min_branch_sep=np.nan,min_hess=np.nan,
                            branch_distance=np.nan,independent_global_check=int(ratio in (.95,1.,1.05))))
    # Previously independent checks at the phase/location corners.
    for r in maprows:
        if r['b_Nomega3'] in (9.,11.) and r['phase_rad'] in (np.pi-.3,np.pi+.3):
            r['independent_global_check']=1
    save('crossing_map.csv',maprows)
    print('MAP',len(maprows),'candidates',sum(r['status']=='candidate' for r in maprows),flush=True)

    # Figure 1: full canonical half-separation domain, not a cropped frequency plot.
    mgrid=np.linspace(-np.pi,np.pi,400)
    hgrid=np.linspace(0,np.pi/2,200)
    mm,hh=np.meshgrid(mgrid,hgrid)
    cells=np.column_stack((mm.ravel(),mm.ravel(),hh.ravel(),hh.ravel()))
    x=st.x_record(n,d,b,ec,phase)
    landscape=evaluate(cells,x,n)[0].reshape(hh.shape)
    plt.figure(figsize=(9,5))
    plt.imshow(np.clip(landscape,0,8),origin='lower',aspect='auto',
               extent=(-np.pi*n,np.pi*n,0,np.pi*n/2),cmap='viridis')
    plt.colorbar(label='Profiled SSE (display clipped at 8)')
    for lab,u in [('A',ua0),('B',ub0)]:
        plt.plot((u[0]+u[1])/2,(u[1]-u[0])/2,'o',ms=7,label=lab)
    plt.xlabel('N midpoint frequency');plt.ylabel('N half separation');plt.legend()
    fig('01_full_reduced_objective.png')

    # Figure 2: terminal adaptive cells at the crossing.
    gc=list(csv.DictReader(open(ROOT/'globality_cells.csv',encoding='utf-8')))
    gc=[r for r in gc if r['epsilon_label']=='crossing']
    plt.figure(figsize=(9,5))
    for stat,color in [('excluded','#9ca3af'),('near_A','#1f77b4'),('near_B','#ff7f0e')]:
        rr=[r for r in gc if r['status']==stat]
        plt.scatter([n*(float(r['m_lo'])+float(r['m_hi']))/2 for r in rr],
                    [n*(float(r['h_lo'])+float(r['h_hi']))/2 for r in rr],
                    s=1 if stat=='excluded' else 6,c=color,label=f'{stat} ({len(rr)})',rasterized=True)
    plt.xlabel('N midpoint frequency');plt.ylabel('N half separation')
    plt.title('Adaptive terminal cells; all outside A/B boxes excluded below J*=Jbest+0.4')
    plt.legend(markerscale=3);fig('02_adaptive_globality_cells.png')

    # Figure 3: independently solved stationary equations.
    sp=list(csv.DictReader(open(ROOT/'stationary_points.csv',encoding='utf-8')))
    sp=[r for r in sp if r['epsilon_label']=='crossing']
    plt.figure(figsize=(9,5))
    for cl,mark,col in [('minimum','o','#1f77b4'),('saddle','x','#d62728'),
                        ('maximum','^','#9467bd')]:
        rr=[r for r in sp if r['classification']==cl and float(r['J'])<8]
        plt.scatter([n*float(r['m']) for r in rr],[n*float(r['h']) for r in rr],
                    marker=mark,c=col,s=12,label=f'{cl} ({len(rr)} below 8)')
    plt.xlabel('N midpoint frequency');plt.ylabel('N half separation')
    plt.legend();fig('03_stationary_points.png')

    yA=np.array([r['J_A'] for r in continuation]);yB=np.array([r['J_B'] for r in continuation])
    plt.figure(figsize=(7,5));plt.plot(ee,yA,label='A');plt.plot(ee,yB,label='B')
    mi=list(csv.DictReader(open(ROOT/'three_to_two_tone_minima.csv',encoding='utf-8')))
    third=[r for r in mi if r['rank']=='3']
    plt.scatter([float(r['epsilon']) for r in third],[float(r['J']) for r in third],
                marker='D',s=40,label='Third best found')
    plt.axvline(ec,color='k',lw=.8);plt.xlabel('epsilon');plt.ylabel('Profiled SSE')
    plt.legend();fig('04_branch_objectives_third_margin.png')
    plt.figure(figsize=(7,5));plt.plot(ee,yA-yB);plt.axhline(0,color='k',lw=.8)
    plt.axvline(ec,color='k',lw=.8);plt.xlabel('epsilon');plt.ylabel('J_A - J_B')
    fig('05_deltaJ_near_crossing.png')
    plt.figure(figsize=(7,5))
    for key,sty in [('u1_A','-'),('u2_A','--'),('u1_B','-'),('u2_B','--')]:
        plt.plot(ee,[r[key] for r in continuation],sty,label=key)
    plt.axvline(ec,color='k',lw=.8);plt.xlabel('epsilon');plt.ylabel('N fitted frequency')
    plt.legend();fig('06_fitted_frequencies.png')
    plt.figure(figsize=(7,5))
    for nn in nn_values:
        rr=[r for r in scaling if r['N']==nn]
        plt.plot([r['d_NDelta'] for r in rr],[r['epsilon_cross'] for r in rr],'-o',ms=3,label=f'N={nn}')
    plt.xlabel('N Delta');plt.ylabel('epsilon crossing');plt.legend(ncol=2)
    fig('07_epsilon_vs_NDelta.png')
    plt.figure(figsize=(7,5))
    rr=[r for r in scaling if r['N']==21]
    plt.loglog([r['d_NDelta'] for r in rr],[r['epsilon_cross'] for r in rr],'-o',label='Exact NLS crossing')
    plt.loglog([r['d_NDelta'] for r in rr],[r['epsilon_tangent'] for r in rr],'--',label='Finite-N tangent: lambda_N d²')
    plt.xlabel('N Delta');plt.ylabel('epsilon crossing');plt.legend()
    fig('08_loglog_scaling.png')
    plt.figure(figsize=(7,5))
    for ph in (0,np.pi/2,np.pi,3*np.pi/2):
        rr=[r for r in maprows if r['phase_rad']==ph and r['strong_amp_ratio']==1 and np.isfinite(r['epsilon_cross'])]
        plt.plot([r['r_location'] for r in rr],[r['epsilon_cross'] for r in rr],'-o',label=f'phase {ph/np.pi:.1f}pi')
    plt.xlabel('Weak-tone location r=(omega3-omega1)/(omega2-omega1)')
    plt.ylabel('Local A/B crossing epsilon');plt.legend()
    fig('09_epsilon_vs_weak_location.png')
    mat=np.zeros((len(phvals),len(bvals)))
    for i,ph in enumerate(phvals):
        for j,bb in enumerate(bvals):
            r=next(z for z in maprows if z['phase_rad']==ph and z['b_Nomega3']==bb and z['strong_amp_ratio']==1.)
            mat[i,j]={'not_found':0,'candidate':1,'degenerate_or_unstable':2}[r['status']]
    plt.figure(figsize=(9,5));plt.imshow(mat,origin='lower',aspect='auto',vmin=0,vmax=2,
                                        extent=(bvals[0]-.5,bvals[-1]+.5,-np.pi/8,2*np.pi+np.pi/8),cmap='viridis')
    cb=plt.colorbar(ticks=(0,1,2));cb.ax.set_yticklabels(('not found','candidate','degenerate'))
    plt.xlabel('Weak-tone normalized position b=N omega3');plt.ylabel('Weak-tone phase (rad)')
    fig('10_crossing_regime_map.png')
    plt.figure(figsize=(7,5))
    for nn in (17,21,31,61):
        rr=[r for r in scaling if r['N']==nn]
        plt.plot([r['d_NDelta'] for r in rr],[r['lambda_numeric'] for r in rr],'-o',ms=3,label=f'N={nn} exact')
        plt.axhline(tangent[nn],ls='--',lw=.8)
    plt.xlabel('N Delta');plt.ylabel('epsilon_c/(N Delta)²')
    plt.legend(ncol=2);fig('11_finite_N_theory_vs_numerics.png')


if __name__=='__main__':main()
