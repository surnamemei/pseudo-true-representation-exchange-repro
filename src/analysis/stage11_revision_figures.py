"""Vector figures for the major scientific revision; no new Monte Carlo."""
from pathlib import Path
import json,csv
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parent;F=R/'paper/figures'
plt.rcParams.update({'font.family':'serif','font.size':8,'axes.labelsize':8,'legend.fontsize':6.5,'pdf.fonttype':42,'svg.fonttype':'none'})
blue='#225A88';orange='#B25632'
p=pd.read_csv(R/'phase_crossing_map.csv');p=p[p.crossing_found==True].sort_values('phi_over_pi')
v=json.loads((R/'stage11_phase_validation.json').read_text());res=pd.read_csv(R/'resolution_sweep.csv')
fig,ax=plt.subplots(1,2,figsize=(7.1,2.65),layout='constrained')
ax[0].plot(p.phi_over_pi,p.lambda_c,color=blue,lw=1.1,label='numerical equal-cost curve')
reps=v['representative_interval_checks'];ax[0].scatter([r['phi_over_pi'] for r in reps],[float(r['root'][2]) for r in reps],marker='s',s=23,c=orange,label='interval-verified local roots',zorder=4)
ax[0].axvline(float(v['central_collision_phi_over_pi']),color='0.45',ls='--',lw=.8,label='central-chart collision')
ax[0].scatter([float(v['cusp_phi_over_pi'])],[float(v['cusp_numerical'][1])],marker='D',s=24,facecolors='white',edgecolors=blue,zorder=5,label='numerical cusp')
ax[0].set(xlabel=r'weak-tone phase $\phi/\pi$',ylabel=r'tangent crossing $\lambda_c$',xlim=(.31,1.025),ylim=(.055,.08))
ax[0].legend(loc='lower right',frameon=False)
ax[1].plot(res.separation_DFT_bins,res.epsilon_cross,'o--',color=blue,mfc='white',ms=4,lw=1,label='numerical full-torus ordering')
ax[1].scatter([2/np.pi],[float(res.iloc[0].epsilon_cross)],s=27,marker='s',color=orange,zorder=4,label='full-domain interval certificate')
ax[1].axhline(1,color='0.5',lw=.7,ls=':');ax[1].text(.66,1.012,'unit generating amplitude',color='0.4',fontsize=6.5)
ax[1].set(xlabel=r'strong-pair separation $z/\pi$ (DFT bins)',ylabel=r'crossing amplitude $\epsilon_c$',xlim=(.58,2.07),ylim=(.18,1.13))
ax[1].legend(loc='lower right',frameon=False)
for a,l in zip(ax,'ab'):a.text(.02,.98,'('+l+')',transform=a.transAxes,va='top',weight='bold')
for ext in ('pdf','svg'):fig.savefig(F/('fig7_phase_resolution.'+ext))
plt.close(fig)
tr=pd.read_csv(R/'noise_branch_trials.csv');d=tr[(tr.N==21)&(tr.z==2)&(tr.snr_db==30)&np.isclose(tr.eta,0)]
base=json.loads((R/'stage2_reference.json').read_text())
fig,ax=plt.subplots(figsize=(3.5,2.65),layout='constrained')
for w,col,mark in [('A',blue,'o'),('B',orange,'s')]:
    g=d[d.branch==w];ax.scatter(g.u1,g.u2,s=11,alpha=.42,marker=mark,c=col,edgecolors='none',label=f'selected {w}: {len(g)}')
    ax.scatter([g.u1.mean()],[g.u2.mean()],s=52,marker='x',lw=1.4,color=col)
    q=base[w+'_at_crossing'];ax.scatter([float(q[0])],[float(q[1])],s=60,marker='+',color='black',lw=1.2,label='noiseless branch optima' if w=='A' else None)
ax.set(xlabel=r'fitted $u_1$',ylabel=r'fitted $u_2$');ax.legend(frameon=False,fontsize=7)
for ext in ('pdf','svg'):fig.savefig(F/('fig6_bimodality.'+ext))
plt.close(fig)
