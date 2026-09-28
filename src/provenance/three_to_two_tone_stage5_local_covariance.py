"""Finite-noise local Hessian/sandwich covariance checks per stationary branch."""
import csv
from pathlib import Path
import numpy as np
import pandas as pd
from three_to_two_tone_stage5_core import times,signal,refine,objective

ROOT=Path(__file__).resolve().parent
trials=pd.read_csv(ROOT/'noise_branch_trials.csv')
N=21;S=times(N);EC=.2481906301722774
A0=[-4.772447353918614,.605170320276862]
B0=[-.057653561926680,12.46341935680166]
out=[]
for eta in (0.,.2):
    e=EC*(1+eta);x=signal(N,2,10,e)
    for label,seed in [('A',A0),('B',B0)]:
        fit=refine(seed,x,S);u=fit['u']
        h=np.empty((2,2));step=2e-4
        for k in range(2):
            du=np.eye(2)[k]*step
            h[:,k]=(objective(u+du,x,S)[1]-objective(u-du,x,S)[1])/(2*step)
        h=(h+h.T)/2
        # Gradient Jacobian with respect to the 2N real noise coordinates.
        L=np.empty((2,2*N));d=1e-5
        for k in range(N):
            basis=np.zeros(N,complex);basis[k]=d
            L[:,k]=(objective(u,x+basis,S)[1]-objective(u,x-basis,S)[1])/(2*d)
            basis[k]=1j*d
            L[:,N+k]=(objective(u,x+basis,S)[1]-objective(u,x-basis,S)[1])/(2*d)
        core=np.linalg.solve(h,L@L.T)@np.linalg.inv(h)
        for snr in (20,30,40):
            sigma2=np.vdot(x,x).real/N*10**(-snr/10)
            pred=sigma2/2*np.trace(core)
            data=trials[(trials.snr_db==snr)&(np.isclose(trials.eta,eta))&(trials.branch==label)]
            empirical=np.trace(np.cov(data[['u1','u2']].to_numpy(),rowvar=False)) if len(data)>2 else np.nan
            out.append(dict(snr_db=snr,eta=eta,branch=label,n_selected=len(data),
                            predicted_local_variance_trace=pred,
                            empirical_conditional_variance_trace=empirical,
                            conditional_over_local=empirical/pred if np.isfinite(empirical) else ''))
with (ROOT/'local_covariance_check.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(out[0]));w.writeheader();w.writerows(out)
print(pd.DataFrame(out).to_string(index=False))
