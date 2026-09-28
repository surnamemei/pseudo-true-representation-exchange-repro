"""Reproducible complex-AWGN branch Monte Carlo and fixed-projector gap model."""
import csv
import math
from pathlib import Path
import numpy as np
from scipy.special import ndtr, ndtri, logsumexp
from scipy.optimize import minimize
from scipy.stats import norm
from three_to_two_tone_stage5_core import times, signal, refine, grid_minima

ROOT=Path(__file__).resolve().parent
N=21; Z=2.; B=10.; EC=.2481906301722774
S=times(N)
A0=np.array([-4.772447353918614,.605170320276862])
B0=np.array([-.057653561926680,12.46341935680166])
RNG=np.random.default_rng(20260926)
TRIALS=400
FIXED=(-.2,-.1,-.05,0.,.05,.1,.2)

def project(u):
    v=np.exp(1j*np.outer(S,u))
    return v@np.linalg.solve(v.conj().T@v,v.conj().T)

def deterministic(eta):
    e=EC*(1+eta)
    x=signal(N,Z,B,e)
    a=refine(A0,x,S); b=refine(B0,x,S)
    return e,x,a,b

def sigma_gap(x,a,b,snr):
    sigma2=(np.vdot(x,x).real/N)*10**(-snr/10)
    pa=project(a['u']);pb=project(b['u'])
    q=pb-pa
    linear=2*sigma2*np.vdot(a['residual']-b['residual'],
                               a['residual']-b['residual']).real
    quadratic=sigma2**2*np.trace(q@q).real
    return sigma2,math.sqrt(max(0,linear+quadratic)),linear,quadratic

def wilson(k,n):
    p=k/n; z=1.959963984540054
    den=1+z*z/n
    center=(p+z*z/(2*n))/den
    half=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return center-half,center+half

def probit_width(rows):
    rows=[r for r in rows if r['n_other']==0]
    eta=np.array([r['eta'] for r in rows],float)
    k=np.array([r['n_A'] for r in rows],float)
    n=np.array([r['n_trials'] for r in rows],float)
    def loss(q):
        p=np.clip(ndtr(q[0]+q[1]*eta),1e-10,1-1e-10)
        return -np.sum(k*np.log(p)+(n-k)*np.log1p(-p))
    fit=minimize(loss,[0,-12],method='BFGS')
    return float((ndtri(.9)-ndtri(.1))/abs(fit.x[1])),fit.x

trial_rows=[];prob_rows=[];model_rows=[];mixture=[];width_rows=[]
for snr in (20,30,40):
    _,x0,a0,b0=deterministic(0)
    _,sd0,_,_=sigma_gap(x0,a0,b0,snr)
    slope=8.259051045764866
    sigma_eta=sd0/(slope*EC)
    adaptive=[q*sigma_eta for q in (-2,-1,-.5,.5,1,2)]
    etas=sorted(set(round(q,8) for q in (*FIXED,*adaptive)))
    snr_probs=[]
    for eta in etas:
        e,x,ad,bd=deterministic(eta)
        sigma2,sd,linvar,quadvar=sigma_gap(x,ad,bd,snr)
        detgap=ad['J']-bd['J']
        predicted=float(ndtr(-detgap/sd))
        a_ct=b_ct=other_ct=fail_ct=audit_ct=audit_miss=0
        gaps=[]
        selected=[]
        for trial in range(TRIALS):
            w=math.sqrt(sigma2/2)*(RNG.normal(size=N)+1j*RNG.normal(size=N))
            y=x+w
            af=refine(ad['u'],y,S)
            bf=refine(bd['u'],y,S)
            fits=[af,bf]
            label='A' if af['J']<=bf['J'] else 'B'
            best=fits[0 if label=='A' else 1]
            search_gap=0.
            audited=(trial%20==0)
            if audited:
                audit_ct+=1
                grid=grid_minima(y,S,grid_size=192,top=24)
                if grid and grid[0]['J']<best['J']-1e-6:
                    search_gap=best['J']-grid[0]['J']
                    audit_miss+=1
                    best=grid[0]
                    if np.linalg.norm(best['u']-af['u'])<.05: label='A'
                    elif np.linalg.norm(best['u']-bf['u'])<.05: label='B'
                    else: label='other'
            if not (af['converged'] and bf['converged']):fail_ct+=1
            a_ct+=label=='A';b_ct+=label=='B';other_ct+=label=='other'
            gap=af['J']-bf['J'];gaps.append(gap)
            selected.append((label,best['u'].copy()))
            trial_rows.append(dict(N=N,z=Z,snr_db=snr,eta=eta,epsilon=e,
              trial=trial,branch=label,u1=best['u'][0],u2=best['u'][1],
              J_global_search=best['J'],J_A=af['J'],J_B=bf['J'],
              observed_gap=gap,error_to_selected_pseudotrue=np.linalg.norm(
                  best['u']-(ad['u'] if label=='A' else bd['u'])) if label!='other' else '',
              squared_error_to_true_strong=np.sum((best['u']-np.array([-Z,Z]))**2),
              audited=int(audited),audit_search_gap=search_gap,
              both_branch_fits_stationary=int(af['converged'] and bf['converged'])))
        lower,upper=wilson(a_ct,TRIALS)
        pr=dict(N=N,z=Z,snr_db=snr,eta=eta,epsilon=e,n_trials=TRIALS,
                n_A=a_ct,n_B=b_ct,n_other=other_ct,n_unconverged=fail_ct,
                n_audited=audit_ct,n_audit_misses=audit_miss,
                P_A=a_ct/TRIALS,P_A_wilson_low=lower,P_A_wilson_high=upper,
                P_B=b_ct/TRIALS,predicted_P_A=predicted)
        prob_rows.append(pr);snr_probs.append(pr)
        model_rows.append(dict(N=N,z=Z,snr_db=snr,eta=eta,epsilon=e,
            deterministic_DeltaJ=detgap,model_sigma_DeltaJ=sd,
            model_linear_variance=linvar,model_quadratic_variance=quadvar,
            observed_DeltaJ_mean=np.mean(gaps),observed_DeltaJ_sd=np.std(gaps,ddof=1),
            P_A_predicted=predicted,P_A_observed=a_ct/TRIALS,
            absolute_probability_error=abs(predicted-a_ct/TRIALS)))
        print('noise',snr,eta,a_ct,b_ct,other_ct,'pred',round(predicted,3),
              'auditmiss',audit_miss,flush=True)
        if eta==0:
            labs=np.array([q[0] for q in selected])
            vals=np.array([q[1] for q in selected])
            pa=np.mean(labs=='A');pb=np.mean(labs=='B')
            aa=vals[labs=='A'];bb=vals[labs=='B']
            ma=aa.mean(axis=0);mb=bb.mean(axis=0)
            ca=np.cov(aa,rowvar=False);cb=np.cov(bb,rowvar=False)
            within=pa*np.trace(ca)*(len(aa)-1)/len(aa)+pb*np.trace(cb)*(len(bb)-1)/len(bb)
            between=pa*pb*np.sum((ma-mb)**2)
            total=np.trace(np.cov(vals,rowvar=False))*(len(vals)-1)/len(vals)
            direction=(mb-ma)/np.linalg.norm(mb-ma)
            t=vals@direction
            ta=aa@direction;tb=bb@direction
            ll1=np.sum(norm.logpdf(t,np.mean(t),np.std(t,ddof=0)))
            ll2=np.sum(logsumexp(np.vstack([
                np.log(pa)+norm.logpdf(t,np.mean(ta),np.std(ta,ddof=0)),
                np.log(pb)+norm.logpdf(t,np.mean(tb),np.std(tb,ddof=0))]),axis=0))
            bic1=2*np.log(len(t))-2*ll1;bic2=5*np.log(len(t))-2*ll2
            mixture.append(dict(N=N,z=Z,snr_db=snr,eta=0,n_A=len(aa),n_B=len(bb),
                mean_A_u1=ma[0],mean_A_u2=ma[1],mean_B_u1=mb[0],mean_B_u2=mb[1],
                between_mode_distance=np.linalg.norm(ma-mb),
                within_variance_trace=within,between_selection_variance_trace=between,
                total_variance_trace=total,decomposition_error=total-within-between,
                mixture_BIC_advantage_over_single_Gaussian=bic1-bic2,
                within_projected_sd_A=np.std(ta,ddof=0),
                within_projected_sd_B=np.std(tb,ddof=0)))
    empirical_width,probit_params=probit_width(snr_probs)
    width_rows.append(dict(snr_db=snr,transition_width_eta_10_90=empirical_width,
        predicted_width_eta_10_90=(ndtri(.9)-ndtri(.1))*sigma_eta,
        fitted_probit_intercept=probit_params[0],fitted_probit_slope=probit_params[1]))

def write_csv(name,rows):
    with (ROOT/name).open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
write_csv('noise_branch_trials.csv',trial_rows)
write_csv('branch_probability.csv',prob_rows)
write_csv('mixture_statistics.csv',mixture)
write_csv('noise_model_fit.csv',model_rows)
write_csv('transition_width.csv',width_rows)
