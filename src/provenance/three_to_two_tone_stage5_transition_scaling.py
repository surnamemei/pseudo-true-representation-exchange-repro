"""Predicted N,z width map and sparse out-of-baseline Monte Carlo checks."""
import csv
from pathlib import Path
import numpy as np
from scipy.special import ndtr, ndtri
from scipy.optimize import minimize
from three_to_two_tone_stage5_core import signal,times,refine,grid_minima

ROOT=Path(__file__).resolve().parent
rows=list(csv.DictReader((ROOT/'crossN_scaling.csv').open()))
RNG=np.random.default_rng(2052026)

def project(s,u):
    v=np.exp(1j*np.outer(s,u))
    return v@np.linalg.solve(v.conj().T@v,v.conj().T)

def fit_at(n,z,e,aseed,bseed):
    s=times(n); x=signal(n,z,10,e)
    return x,refine(aseed,x,s),refine(bseed,x,s)

out=[]; mc=[]
checks={(11,2.,30),(41,2.,30),(21,1.,30),(21,1.,40),(21,.6,40)}
for r in rows:
    if r['status']!='numerical_stationary_crossing':continue
    n=int(r['N']);z=float(r['z']);e=float(r['epsilon_cross'])
    if z not in (.4,.6,1.,1.5,2.):continue
    s=times(n);as_=np.array([float(r['u_A1']),float(r['u_A2'])]);bs_=np.array([float(r['u_B1']),float(r['u_B2'])])
    x,a,b=fit_at(n,z,e,as_,bs_)
    h=max(1e-6,e*1e-4)
    _,ap,bp=fit_at(n,z,e+h,as_,bs_)
    _,am,bm=fit_at(n,z,e-h,as_,bs_)
    slope=(ap['J']-bp['J']-am['J']+bm['J'])/(2*h)
    for snr in (30,40):
        sigma2=np.vdot(x,x).real/n*10**(-snr/10)
        q=project(s,b['u'])-project(s,a['u'])
        var=2*sigma2*np.linalg.norm(a['residual']-b['residual'])**2+sigma2**2*np.trace(q@q).real
        sd=np.sqrt(var)
        width_e=2*ndtri(.9)*sd/abs(slope)
        width_eta=width_e/e
        item=dict(N=n,z=z,snr_db=snr,epsilon_cross=e,DeltaJ_slope=slope,
            sigma_gap=sd,predicted_width_epsilon_10_90=width_e,
            predicted_width_eta_10_90=width_eta,mc_status='not_run',
            mc_width_eta_10_90='',mc_other_fraction='',mc_audit_misses='')
        if (n,z,snr) in checks:
            eta_scale=sd/(abs(slope)*e)
            etas=[-1.5*eta_scale,-.75*eta_scale,0,.75*eta_scale,1.5*eta_scale]
            counts=[];other=0;audits=0;miss=0
            for eta in etas:
                ee=e*(1+eta)
                xx,aa,bb=fit_at(n,z,ee,as_,bs_)
                ca=cb=co=0
                for k in range(200):
                    y=xx+np.sqrt(sigma2/2)*(RNG.normal(size=n)+1j*RNG.normal(size=n))
                    af=refine(aa['u'],y,s);bf=refine(bb['u'],y,s)
                    label='A' if af['J']<=bf['J'] else 'B'
                    best=af if label=='A' else bf
                    if k%10==0:
                        audits+=1
                        gr=grid_minima(y,s,grid_size=192,top=24)
                        if gr and gr[0]['J']<best['J']-1e-6:
                            miss+=1;best=gr[0]
                            if np.linalg.norm(best['u']-af['u'])<.05:label='A'
                            elif np.linalg.norm(best['u']-bf['u'])<.05:label='B'
                            else:label='other'
                    ca+=label=='A';cb+=label=='B';co+=label=='other'
                other+=co
                counts.append((eta,ca,cb,co))
                mc.append(dict(N=n,z=z,snr_db=snr,eta=eta,n_trials=200,
                               n_A=ca,n_B=cb,n_other=co,P_A=ca/200))
                print('scaling MC',n,z,snr,eta,ca,cb,co,flush=True)
            if other==0:
                eta_arr=np.array([q[0] for q in counts]);k_arr=np.array([q[1] for q in counts]);
                def loss(beta):
                    p=np.clip(ndtr(beta[0]+beta[1]*eta_arr),1e-9,1-1e-9)
                    return -sum(k_arr*np.log(p)+(200-k_arr)*np.log1p(-p))
                fit=minimize(loss,[0,-1/eta_scale],method='BFGS')
                item['mc_width_eta_10_90']=2*ndtri(.9)/abs(fit.x[1])
            item['mc_status']='sparse_200_per_point_global_audited'
            item['mc_other_fraction']=other/1000
            item['mc_audit_misses']=miss
            print('scaling summary',n,z,snr,'width',item['mc_width_eta_10_90'],
                  'pred',width_eta,'other',other,'miss',miss,'audits',audits,flush=True)
        out.append(item)

def save(name,data):
    with (ROOT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
save('transition_width_scaling.csv',out)
save('transition_scaling_trials.csv',mc)
