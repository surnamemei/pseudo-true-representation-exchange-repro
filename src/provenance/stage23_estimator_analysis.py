"""Local small-noise covariance and global mixture diagnostics from frozen trials.

Uses full six-real-parameter residual Hessian, including misspecification
curvature. The result is a local linearization, not a global CRB.
"""
import csv,json,time
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from three_to_two_tone_stage5_core import refine,signal,times

ROOT=Path(__file__).resolve().parent;N=21;S=times(N);EC=.2481906301722774
A=np.array([-4.772447353918614,.605170320276862]);B=np.array([-.057653561926680,12.46341935680166])
def local_cov(x,u,sigma2):
    v=np.exp(1j*np.outer(S,u));a=np.linalg.lstsq(v,x,rcond=None)[0];r=x-v@a
    jac=np.column_stack([v[:,0],1j*v[:,0],v[:,1],1j*v[:,1],1j*S*a[0]*v[:,0],1j*S*a[1]*v[:,1]])
    second=np.zeros((N,6,6),complex)
    for k in range(2):
        second[:,2*k,4+k]=second[:,4+k,2*k]=1j*S*v[:,k]
        second[:,2*k+1,4+k]=second[:,4+k,2*k+1]=-S*v[:,k]
        second[:,4+k,4+k]=-S*S*a[k]*v[:,k]
    jj=(jac.conj().T@jac).real
    mat=jj-np.einsum('n,nij->ij',np.conj(r),second).real
    inv=np.linalg.inv(mat);cov=sigma2/2*inv@jj@inv
    # Independent finite differences of the full objective gradient.
    theta=np.r_[a[0].real,a[0].imag,a[1].real,a[1].imag,u]
    def grad(t):
        aa=np.array([t[0]+1j*t[1],t[2]+1j*t[3]]);vv=np.exp(1j*np.outer(S,t[4:]));rr=x-vv@aa
        j=np.column_stack([vv[:,0],1j*vv[:,0],vv[:,1],1j*vv[:,1],1j*S*aa[0]*vv[:,0],1j*S*aa[1]*vv[:,1]])
        return -2*(j.conj().T@rr).real
    step=2e-5
    numeric=np.column_stack([(grad(theta+np.eye(6)[k]*step)-grad(theta-np.eye(6)[k]*step))/(2*step) for k in range(6)])
    discrepancy=np.max(np.abs(numeric-2*mat))/np.max(np.abs(2*mat))
    return cov[4:,4:],float(np.linalg.eigvalsh(mat)[0]),float(discrepancy)
def run():
    start=time.perf_counter();data=list(csv.DictReader((ROOT/'noise_branch_trials.csv').open()));out=[];figdata={};checks=[]
    for snr,eta in sorted(set((int(r['snr_db']),float(r['eta'])) for r in data)):
        rows=[r for r in data if int(r['snr_db'])==snr and float(r['eta'])==eta]
        e=EC*(1+eta);x=signal(N,2,10,e);fa=refine(A,x,S);fb=refine(B,x,S);sigma2=np.vdot(x,x).real/N*10**(-snr/10)
        covs=[]
        for f in (fa,fb):
            cv,eig,err=local_cov(x,f['u'],sigma2);covs.append(cv);checks.append(dict(snr=snr,eta=eta,hessian_min=eig,fd_error=err))
        u=np.array([[float(r['u1']),float(r['u2'])] for r in rows]);labels=np.array([r['branch'] for r in rows]);n=len(rows)
        refs=[fa['u'],fb['u']];probs=[np.mean(labels==j) for j in ('A','B')]
        reference=refs[0 if eta<=0 else 1] # tie at eta=0 is explicitly A.
        rec=dict(snr_db=snr,eta=eta,n=n,n_A=int(sum(labels=='A')),n_B=int(sum(labels=='B')),P_A=probs[0],
                 global_MSE_to_noiseless_winner=float(np.mean(np.sum((u-reference)**2,axis=1))),
                 global_variance_trace=float(np.var(u,axis=0).sum()),
                 weighted_local_covariance_trace=float(sum(p*np.trace(c) for p,c in zip(probs,covs))))
        means=[];within=0.;own=0.
        for j,name in enumerate(('A','B')):
            q=u[labels==name];cv=covs[j]
            mse=float(np.mean(np.sum((q-refs[j])**2,axis=1))) if len(q) else float('nan')
            mean=q.mean(axis=0) if len(q) else refs[j]
            vari=float(np.var(q,axis=0).sum()) if len(q) else 0.
            rec.update({name+'_conditional_MSE_own':mse,name+'_conditional_variance_trace':vari,
                 name+'_predicted_local_trace':float(np.trace(cv)),name+'_predicted_cov11':float(cv[0,0]),
                 name+'_predicted_cov12':float(cv[0,1]),name+'_predicted_cov22':float(cv[1,1]),
                 name+'_pseudotrue_u1':refs[j][0],name+'_pseudotrue_u2':refs[j][1]})
            means.append(mean);within+=probs[j]*vari
            if len(q):own+=probs[j]*mse
        between=probs[0]*probs[1]*np.sum((means[0]-means[1])**2)
        rec.update(weighted_conditional_MSE_own=own,within_variance_trace=within,between_variance_trace=float(between),
                   decomposition_error=rec['global_variance_trace']-within-between)
        out.append(rec)
        if eta==0:figdata[snr]=(u,labels,refs)
    with (ROOT/'estimator_mode_statistics.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=out[0]);w.writeheader();w.writerows(out)
    with (ROOT/'local_covariance_validation.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=checks[0]);w.writeheader();w.writerows(checks)
    fig,ax=plt.subplots(1,3,figsize=(10,3.05));colors={20:'#166da2',30:'#d77519',40:'#32834b'}
    for snr in (20,30,40):
        q=[r for r in out if r['snr_db']==snr];eta=[r['eta'] for r in q]
        ax[0].plot(eta,[r['P_A'] for r in q],'.-',color=colors[snr],label=f'{snr} dB')
    ax[0].set(xlabel=r'Normalized offset $\eta$',ylabel=r'$P(\widehat{j}=\mathcal{A})$',title='(a) Mode selection');ax[0].legend(frameon=False,fontsize=8)
    u,labs,refs=figdata[30];direction=(refs[1]-refs[0])/np.linalg.norm(refs[1]-refs[0]);coord=(u-refs[0])@direction
    for name,color in [('A','#166da2'),('B','#d77519')]:ax[1].hist(coord[labs==name],bins=18,alpha=.75,color=color,label=rf'$\mathcal{{{name}}}$')
    ax[1].set(xlabel='Frequency-pair coordinate along mode gap',ylabel='Trial count',title=r'(b) 30 dB, $\eta=0$');ax[1].legend(frameon=False,fontsize=8)
    q=[r for r in out if r['snr_db']==30];eta=[r['eta'] for r in q]
    ax[2].semilogy(eta,[r['global_MSE_to_noiseless_winner'] for r in q],'.-',label='Global / noiseless winner')
    ax[2].semilogy(eta,[r['weighted_conditional_MSE_own'] for r in q],'.-',label='Selected mode / own target')
    ax[2].semilogy(eta,[r['weighted_local_covariance_trace'] for r in q],'--',label='Local covariance prediction')
    ax[2].set(xlabel=r'Normalized offset $\eta$',ylabel=r'Frequency-pair squared error',title='(c) 30 dB error scales');ax[2].legend(frameon=False,fontsize=6.5)
    for a in ax:a.grid(alpha=.15)
    fig.tight_layout();fig.savefig(ROOT/'estimator_mode_mixture.pdf');plt.close(fig)
    summary=dict(source_trials=len(data),source_audits=sum(int(r['audited']) for r in data),
                 source_audit_misses=sum(float(r['audit_search_gap'])>1e-6 for r in data),
                 new_noise_draws=0,max_hessian_fd_relative_error=max(r['fd_error'] for r in checks),
                 min_full_hessian_half_eigenvalue=min(r['hessian_min'] for r in checks),
                 max_variance_identity_error=max(abs(r['decomposition_error']) for r in out),
                 seconds=time.perf_counter()-start,crossing=[r for r in out if r['eta']==0])
    (ROOT/'estimator_analysis_validation.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))
if __name__=='__main__':run()
