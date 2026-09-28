"""Ten reproducible Stage 5 scientific figures."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from three_to_two_tone_stage5_core import times,signal,refine

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'three_to_two_tone_stage5_figures';OUT.mkdir(exist_ok=True)
cross=pd.read_csv(ROOT/'crossN_scaling.csv')
asym=pd.read_csv(ROOT/'asymmetry_sweep.csv')
region=pd.read_csv(ROOT/'existence_region.csv')
prob=pd.read_csv(ROOT/'branch_probability.csv')
trials=pd.read_csv(ROOT/'noise_branch_trials.csv')
mix=pd.read_csv(ROOT/'mixture_statistics.csv')
width=pd.read_csv(ROOT/'transition_width.csv')
scaling=pd.read_csv(ROOT/'transition_width_scaling.csv')

plt.rcParams.update({'figure.dpi':150,'savefig.dpi':180,'font.size':10,
                     'axes.grid':True,'grid.alpha':.25})
def save(i,name):
    plt.tight_layout();plt.savefig(OUT/f'{i:02d}_{name}.png');plt.close()

# 1: Lambda_N approaches a fixed normalized-sampling limit.
tan=cross.drop_duplicates('N').sort_values('N')
fit=np.polyfit(1/tan.N.to_numpy()[-4:]**2,tan.lambda_N.to_numpy()[-4:],1)
plt.figure(figsize=(6,4));plt.plot(tan.N,tan.lambda_N,'o-',label=r'$\lambda_N$')
plt.axhline(fit[1],ls='--',c='k',label=fr'$1/N^2$ extrapolation {fit[1]:.7f}')
plt.xlabel('odd record length N');plt.ylabel(r'$\lambda_N$');plt.legend();save(1,'lambda_N_vs_N')

# 2: scaling collapse, numerical stationary crossings versus local law.
plt.figure(figsize=(6,4))
for n,grp in cross[cross.N<=41].groupby('N'):
    grp=grp.dropna(subset=['epsilon_cross']);plt.plot(grp.z,grp.epsilon_cross/grp.z**2,'o-',label=f'N={n}')
plt.xlabel(r'$z=N\Delta$');plt.ylabel(r'$\epsilon_c/z^2$');plt.legend(ncol=2);save(2,'scaling_collapse')

# 3 and 4: branch gaps through epsilon, with the vertical crossings.
s=times(21);a0=[-4.772447353918614,.605170320276862];b0=[-.057653561926680,12.46341935680166]
for idx,family,levels,name in [(3,'strong_amplitude_ratio',[.8,1,1.2],'amplitude_asymmetry'),
                                (4,'weak_phase',[np.pi/2,3*np.pi/4,np.pi],'phase_asymmetry')]:
    plt.figure(figsize=(6,4))
    for value in levels:
        row=asym[(asym.family==family)&(np.isclose(asym.value,value))].iloc[0]
        ee=np.linspace(max(.06,row.epsilon_c-.08),row.epsilon_c+.08,41)
        gap=[]
        for e in ee:
            x=signal(21,2,10,e,phase=value if family=='weak_phase' else np.pi,
                     ratio=value if family=='strong_amplitude_ratio' else 1)
            gap.append(refine(a0,x,s)['J']-refine(b0,x,s)['J'])
        plt.plot(ee,gap,label=f'{value:.3g}' + (' rad' if family=='weak_phase' else ''))
        plt.scatter([row.epsilon_c],[0],s=20)
    plt.axhline(0,c='k',lw=.8);plt.xlabel(r'$\epsilon$');plt.ylabel(r'$J_A-J_B$')
    plt.legend(title='phase' if family=='weak_phase' else 'A2/A1');save(idx,name)

# 5: compact measured open region; cells are numerical, not interval proofs.
pivot=region.pivot(index='ratio',columns='b',values='epsilon_c').sort_index()
plt.figure(figsize=(6,4));im=plt.imshow(pivot.to_numpy(),origin='lower',aspect='auto',cmap='viridis')
plt.xticks(range(len(pivot.columns)),[str(x) for x in pivot.columns]);
plt.yticks(range(len(pivot.index)),[str(x) for x in pivot.index]);
plt.xlabel('weak frequency b');plt.ylabel('strong amplitude ratio');
plt.colorbar(im,label=r'numerical $\epsilon_c$');save(5,'existence_region')

# 6: branch selection probability with binomial uncertainty.
plt.figure(figsize=(6,4))
for snr,grp in prob.groupby('snr_db'):
    grp=grp.sort_values('eta');plt.plot(grp.eta,grp.P_A,'o-',label=f'{snr} dB')
    plt.fill_between(grp.eta,grp.P_A_wilson_low,grp.P_A_wilson_high,alpha=.12)
plt.axvline(0,c='k',ls='--',lw=.8);plt.xlabel(r'$\eta=(\epsilon-\epsilon_c)/\epsilon_c$')
plt.ylabel('P(select A)');plt.legend();save(6,'branch_probability')

# 7: directly observed bimodality of the second fitted frequency.
plt.figure(figsize=(6,4));at=trials[(trials.snr_db==30)&(trials.eta==0)]
for branch,color in [('A','tab:blue'),('B','tab:orange')]:
    plt.hist(at.loc[at.branch==branch,'u2'],bins=35,alpha=.65,label=branch,color=color)
plt.xlabel(r'fitted second frequency $u_2=N\hat\omega_2$');plt.ylabel('trials');
plt.legend();save(7,'frequency_bimodality')

# 8: exact MSE identity needs within, between, and squared bias to a reference.
plt.figure(figsize=(6,4));xx=np.arange(len(mix));
plt.bar(xx,mix.within_variance_trace,label='within branch')
plt.bar(xx,mix.between_selection_variance_trace,bottom=mix.within_variance_trace,
        label='between branch')
plt.bar(xx,mix.bias_to_true_strong_sq,
        bottom=mix.within_variance_trace+mix.between_selection_variance_trace,
        label='squared bias to true strong pair')
plt.xticks(xx,[f'{q} dB' for q in mix.snr_db]);
plt.ylabel(r'MSE of fitted $(u_1,u_2)$ to $(-2,2)$');
plt.legend();save(8,'MSE_decomposition')

# 9: fixed-projector Gaussian cost-gap model against Monte Carlo.
plt.figure(figsize=(5,5));plt.plot([0,1],[0,1],'k--',lw=.8)
for snr,grp in prob.groupby('snr_db'):
    plt.errorbar(grp.predicted_P_A,grp.P_A,
       yerr=[np.maximum(0,grp.P_A-grp.P_A_wilson_low),
             np.maximum(0,grp.P_A_wilson_high-grp.P_A)],
       fmt='o',ms=4,label=f'{snr} dB',alpha=.8)
plt.xlabel('predicted P(A)');plt.ylabel('observed P(A)');plt.legend();save(9,'predicted_vs_observed')

# 10: SNR width, plus a distinct panel for cross-N / z predictions and sparse checks.
fig,ax=plt.subplots(1,2,figsize=(10,4))
ax[0].semilogy(width.snr_db,width.transition_width_eta_10_90,'o-',label='MC probit')
ax[0].semilogy(width.snr_db,width.predicted_width_eta_10_90,'s--',label='gap model')
ax[0].set(xlabel='SNR (dB)',ylabel='10–90% width in eta');ax[0].legend()
for n,grp in scaling[(scaling.snr_db==40)&(scaling.N.isin([11,21,41]))].groupby('N'):
    ax[1].plot(grp.z,grp.predicted_width_eta_10_90,'o-',label=f'N={n}, model')
checks=scaling[(scaling.snr_db==40)&scaling.mc_width_eta_10_90.notna()]
ax[1].scatter(checks.z,checks.mc_width_eta_10_90,marker='x',s=70,c='k',label='sparse MC')
ax[1].set(xlabel=r'$z=N\Delta$',ylabel='10–90% width in eta');ax[1].legend(fontsize=8)
fig.tight_layout();fig.savefig(OUT/'10_transition_width.png');plt.close(fig)
