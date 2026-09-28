"""Freeze Stage-7 claim/figure edits and regenerate the cross-N evidence plot."""
import csv
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parent

old=(ROOT/'paper_claim_ledger.md').read_text(encoding='utf-8')
old=old.replace('# Paper claim ledger (frozen for manuscript drafting)',
                '# Paper claim ledger (Stage 7 freeze)')
old=old.replace('A separate interval calculation certifies a **global** exchange at N=21, z=2. Other lengths, asymmetries, and noisy selections are numerical evidence.',
                'Independent interval calculations certify a **global** exchange at both N=21 and N=31, each at z=2 under matched normalized geometry. The other tested lengths, asymmetries, and noisy selections are numerical evidence.')
needle='| The N=21,z=2 crossing is transverse and not caused by fitted-tone coalescence.'
insert=('| At N=31,z=2 under the same normalized geometry, one A/B equal-cost crossing lies in the recorded ±10⁻¹⁴ epsilon interval; A/B are the only globally competitive fits throughout that bracket. | INTERVAL-CERTIFIED | [N31_global_certificate.json](N31_global_certificate.json), [N31_global_certificate.md](N31_global_certificate.md), [N31_interval_cells.csv](N31_interval_cells.csv) | Main Fig. 4 | “A second independent record length has a full-torus noncoalescent global exchange at the stated point.” | “The exchange is globally certified for every N or for a connected cross-N curve.” |\n'
        '| The N=31,z=2 crossing is transverse, with positive Hessians and finite fitted amplitudes. | INTERVAL-CERTIFIED | [N31_crossing_certificate.csv](N31_crossing_certificate.csv), [N31_global_certificate.json](N31_global_certificate.json) | Main Fig. 4 | “The interval slope lies in `[12.2910730811,12.2910731035]`, and the confluent strip is excluded.” | “The quartic small-z expansion is globally certified at z=2.” |\n')
if needle not in old:raise RuntimeError('ledger insertion anchor missing')
old=old.replace(needle,insert+needle)
old=old.replace('The A/B crossing persists at N=11,15,21,31,41 over sampled z. | NUMERICAL',
                'The A/B stationary crossing persists numerically at N=11,15,21,31,41 over sampled z; the z=2 global points at N=21 and 31 have separate interval certificates. | NUMERICAL')
old=old.replace('“Numerical equal-cost crossings and independent grid top-two checks persist across the five tested lengths.” | “A global 3→2 exchange theorem holds for every odd N.”',
                '“Numerical crossings persist across the five tested lengths; globality at z=2 is separately certified for N=21 and 31.” | “A global 3→2 exchange theorem holds for every odd N.”')
old+='\n**Stage-7 boundary:** Two different record lengths are certified at the same normalized z=2 geometry. The analytical small-z theorem is still specifically N=21; the N=31 tangent coefficients and other N values are numerical unless separately stated. No connected global interval in z between the small-z theorem and z=2 has been proved.\n'
(ROOT/'updated_claim_ledger.md').write_text(old,encoding='utf-8')

fig=(ROOT/'figure_plan.md').read_text(encoding='utf-8')
fig=fig.replace('# Main-paper figure reduction','# Main-paper figure reduction — Stage 7 freeze')
fig=fig.replace('fixed legend code throughout: solid black bracket/box = interval-certified, dashed colored line = analytical local expansion, colored marker/line = numerical, shaded band = Monte Carlo uncertainty.',
                'fixed legend code throughout: filled black square = interval-certified global point, outlined square = independently certified pointwise result, dashed colored line = proved local analytical expansion, dotted extension = unvalidated extrapolation, open colored circle = numerical point, shaded band = Monte Carlo uncertainty.')
fig=fig.replace('All cross-N/asymmetry points are **numerical**; baseline alone is interval-certified. Do not describe the 25-cell map as a certified region.',
                'At z=2, **N=21 and N=31 are separate interval-certified global points**; N=11,15,41 and off-baseline N=31 points remain numerical. Never connect the certified markers as a certified all-N curve. Asymmetry and the 25-cell map remain numerical.')
fig=fig.replace('(a) `λ_N` for N=11,15,21,31,41 (81/161/321 as small inset), or `εc/z²` across N.',
                '(a) Matched-z=2 cross-N critical amplitudes: black filled squares at N=21 and N=31, open circles at N=11,15,41; N-specific local-law quartic predictions as dashed/dotted comparison. (b) `εc(z)` for N=21 and N=31, with certified point markers, numerical markers, and the analytical local expansions separately styled. If retaining asymmetry, move it to the supplement.')
fig=fig.replace('(b) numerical `εc` versus amplitude ratio; (c) phase perturbation with φ=0,π/4 explicitly marked “no tracked bracket found.”',
                'Amplitude-ratio and phase panels, including φ=0,π/4 “no tracked bracket found,” move to the supplement.')
fig+='\n## New plotted evidence artifact\n\n[Stage-7 cross-N scaling panel](stage7_crossN_evidence.png) is generated from [crossN_scaling.csv](crossN_scaling.csv) by [stage7_manuscript_update.py](stage7_manuscript_update.py). Its filled black squares are **only** the two N=21 and N=31 z=2 global certificates. Dashed small-z laws and dotted extrapolations are not full-domain proofs. Captions must state this distinction.\n'
(ROOT/'updated_figure_plan.md').write_text(fig,encoding='utf-8')

rows=list(csv.DictReader((ROOT/'crossN_scaling.csv').open(newline='',encoding='utf-8')))
z2=sorted((r for r in rows if r['epsilon_cross'] and float(r['z'])==2 and int(r['N']) in (11,15,21,31,41)),key=lambda r:int(r['N']))
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,ax=plt.subplots(1,2,figsize=(10.2,4.1),constrained_layout=True)
for r in z2:
    n=int(r['N']);e=float(r['epsilon_cross']);q=float(r['epsilon_quartic'])
    if n in (21,31):
        ax[0].plot(n,e,'s',ms=8,color='black',zorder=4)
    else:
        ax[0].plot(n,e,'o',ms=7,mfc='white',mec='#3867a6',mew=1.5,zorder=3)
    ax[0].plot(n,q,'x',ms=6,color='#b15b3a',zorder=2)
ax[0].plot([],[],'ks',label='Certified global, z=2')
ax[0].plot([],[],'o',mfc='white',mec='#3867a6',label='Numerical crossing')
ax[0].plot([],[],'x',color='#b15b3a',label='Local quartic prediction')
ax[0].set(xlabel='Record length N',ylabel='Critical omitted-tone amplitude εc',title='Matched normalized geometry: z=2')
ax[0].legend(frameon=False,loc='lower right',fontsize=8)

for n,color in [(21,'#275d8a'),(31,'#a85c34')]:
    rr=sorted((r for r in rows if int(r['N'])==n and r['epsilon_cross']),key=lambda r:float(r['z']))
    zz=np.array([float(r['z']) for r in rr]);ee=np.array([float(r['epsilon_cross']) for r in rr])
    lam=float(rr[0]['lambda_N']);c=float(rr[0]['c_N'])
    for z,e in zip(zz,ee):
        if z==2:
            ax[1].plot(z,e,'s',ms=8,color='black',zorder=5)
        elif n==21 and z in (1.0,1.5):
            ax[1].plot(z,e,'s',ms=6,mfc='white',mec='black',mew=1.3,zorder=4)
        else:
            ax[1].plot(z,e,'o',ms=5,mfc='white',mec=color,mew=1.2,zorder=3)
    local=np.linspace(0,0.6,100)
    extrap=np.linspace(0.6,2,150)
    style='--' if n==21 else '-.'
    meaning='proved local law' if n==21 else 'numerical tangent formula'
    ax[1].plot(local,lam*local**2+c*local**4,style,color=color,lw=1.5,label=f'N={n} {meaning}')
    ax[1].plot(extrap,lam*extrap**2+c*extrap**4,':',color=color,lw=1.2)
ax[1].set(xlabel='Normalized spacing z=NΔ',ylabel='Critical amplitude εc',title='Local law vs finite-spacing crossings')
ax[1].legend(frameon=False,loc='upper left',fontsize=8)
zoom=ax[1].inset_axes([0.58,0.17,0.36,0.28])
for n,col in [(21,'#275d8a'),(31,'#a85c34')]:
    q=next(r for r in z2 if int(r['N'])==n)
    zoom.plot(2,float(q['epsilon_cross']),'s',ms=6,color='black')
    zoom.annotate(f'N={n}',(2,float(q['epsilon_cross'])),xytext=(6,3 if n==31 else -11),
                  textcoords='offset points',fontsize=7,color=col)
zoom.set_xlim(1.985,2.13);zoom.set_ylim(0.2477,0.2497)
zoom.set_xticks([2]);zoom.set_yticks([0.248,0.249])
zoom.tick_params(labelsize=7)
fig.savefig(ROOT/'stage7_crossN_evidence.png',dpi=240)
fig.savefig(ROOT/'stage7_crossN_evidence.pdf')
plt.close(fig)
