"""Ten static Stage-2 diagnostic figures; certifications reside in CSV/JSON."""
import csv,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

ROOT=Path(__file__).resolve().parent
FIG=ROOT/'three_to_two_tone_stage2_figures';FIG.mkdir(exist_ok=True)
ref=json.loads((ROOT/'stage2_reference.json').read_text())
cert=json.loads((ROOT/'validated_globality_certificate.json').read_text())
audit=json.loads((ROOT/'stage2_cell_audit.json').read_text())
ec=float(ref['epsilon_cross']);A=np.array([float(x) for x in ref['A_at_crossing']]);B=np.array([float(x) for x in ref['B_at_crossing']])

def read(name):return list(csv.DictReader((ROOT/name).open(newline='')))
def save(name):plt.tight_layout();plt.savefig(FIG/name,dpi=180);plt.close()

# Analytic finite-sum reduced objective on a plotted grid.
t=np.arange(-10,11);s=t/21
x=2*np.cos(2*s)-ec*np.exp(1j*10*s)
u1=np.linspace(-8,2,170);u2=np.linspace(-1,16,190)
U1,U2=np.meshgrid(u1,u2,indexing='ij')
D=np.zeros_like(U1,dtype=complex);h1=np.zeros_like(U1,dtype=complex);h2=np.zeros_like(U1,dtype=complex)
for tt,ss in zip(t,s):
    D+=np.exp(1j*(U2-U1)*ss)
    h1+=x[tt+10]*np.exp(-1j*U1*ss)
    h2+=x[tt+10]*np.exp(-1j*U2*ss)
den=21**2-abs(D)**2
J=np.full_like(U1,np.nan,dtype=float)
ok=(U2-U1>.3)&(den>1e-6)
J[ok]=np.vdot(x,x).real-(21*(abs(h1[ok])**2+abs(h2[ok])**2)-2*np.real(np.conj(h1[ok])*D[ok]*h2[ok]))/den[ok]
plt.figure(figsize=(7,5));plt.contourf(U1,U2,np.clip(J,0,5),levels=np.linspace(0,5,31),cmap='viridis')
plt.colorbar(label='Reduced squared error J');plt.scatter([A[0],B[0]],[A[1],B[1]],c=['tomato','cyan'],edgecolor='k',s=65)
plt.text(A[0]+.2,A[1]+.4,'A');plt.text(B[0]+.2,B[1]+.4,'B')
plt.xlabel('u1 = Nν1');plt.ylabel('u2 = Nν2');plt.title('1. Finite-record objective at the crossing')
save('01_certified_objective_landscape.png')

plt.figure(figsize=(7,5));plt.scatter([A[0],B[0]],[A[1],B[1]],c=['tomato','royalblue'],s=80)
ax=plt.gca()
for name,c,col,r in [('A',A,'tomato',.0003),('B',B,'royalblue',.01)]:
    ax.add_patch(Rectangle(c-2.00000001,4.00000002,4.00000002,fill=False,edgecolor=col,lw=2))
    plt.annotate(f'{name}: Krawczyk radius {r:g}',c,xytext=(c[0]-2,c[1]+3),arrowprops=dict(arrowstyle='->',color=col),color=col)
plt.xlim(-8,3);plt.ylim(-3,17);plt.xlabel('u1');plt.ylabel('u2');plt.title('2. Accepted frequency boxes and certified roots')
save('02_certified_AB_boxes.png')

rr=[r for r in read('interval_boxes.csv') if r['label']=='crossing']
mc=np.array([(float(r['m_lo'])+float(r['m_hi']))/2 for r in rr]);hc=np.array([(float(r['h_lo'])+float(r['h_hi']))/2 for r in rr]);
yes=np.array([r['status']!='excluded' for r in rr])
plt.figure(figsize=(8,4.5));plt.scatter(mc[~yes],hc[~yes],s=1,c='lightgray',rasterized=True,label='interval-excluded')
plt.scatter(mc[yes],hc[yes],s=8,c='red',label='inside A/B inner boxes')
plt.axhspan(0,.03,color='orange',alpha=.18,label='confluent strip')
plt.xlabel('m');plt.ylabel('h');plt.title('3. Full-torus cell partition at crossing');plt.legend(markerscale=3)
save('03_excluded_domain_globality_map.png')

cont=read('branch_continuation.csv');ee=np.array([float(r['epsilon']) for r in cont]);diff=np.array([float(r['DeltaJ']) for r in cont])
plt.figure(figsize=(7,4));plt.plot(ee,diff,label='tracked branch difference');plt.axhline(0,color='k',lw=.8);plt.axvline(ec,color='red',ls='--',label='interval-bracketed root')
plt.scatter([float(r['epsilon']) for r in ref['records']],[float(r['DeltaJ']) for r in ref['records']],color='red',s=30)
plt.xlabel('ε');plt.ylabel('J_A − J_B');plt.title('4. Transverse objective crossing');plt.legend();save('04_DeltaJ_certified_crossing.png')

h=cert['hessian_eigenvalue_bounds'];labs=[r['label'][0].upper()+r['branch'] for r in h]
mins=[float(r['lambda_min_lower'].split(',')[0].strip('[')) for r in h]
maxs=[float(r['lambda_max_upper'].split(',')[0].strip('[')) for r in h]
plt.figure(figsize=(7,4));xx=np.arange(len(h));plt.semilogy(xx,mins,'o-',label='λmin lower');plt.semilogy(xx,maxs,'s-',label='λmax upper')
plt.xticks(xx,labs);plt.ylabel('Hessian eigenvalue bound');plt.title('5. Positive-definite interval Hessians');plt.legend();save('05_Hessian_eigenvalue_bounds.png')

plt.figure(figsize=(6,4));plt.bar([z['label'] for z in audit],[float(z['coalescent_J_min_margin']) for z in audit],color='orange')
plt.axhline(0,color='k',lw=.8);plt.ylabel('J lower − global incumbent upper');plt.title('6. Exclusion of h ≤ 0.03, including confluent limit')
save('06_coalescent_strip_exclusion.png')

sc=read('finiteN_scaling.csv');nset=sorted(set(int(r['N']) for r in sc));plt.figure(figsize=(6,5))
for n in nset:
    row=[r for r in sc if int(r['N'])==n]
    plt.plot([float(r['epsilon_leading']) for r in row],[float(r['epsilon_cross']) for r in row],'.-',label=f'N={n}')
lim=plt.xlim();plt.plot(lim,lim,'k--',lw=1,label='equality');plt.xlabel('λ_N d² prediction');plt.ylabel('numerical εc');plt.title('7. Leading law vs numerical crossings');plt.legend(ncol=2,fontsize=8)
save('07_theory_vs_numerical_epsilon.png')

co=read('asymptotic_coefficients.csv');plt.figure(figsize=(6,4));plt.plot([int(r['N']) for r in co],[float(r['lambda_N']) for r in co],'o-')
plt.axhline(float(co[0]['lambda_infinity_estimate']),ls='--',color='red',label='N→∞ extrapolation')
plt.xscale('log');plt.xlabel('N');plt.ylabel('λ_N');plt.title('8. Finite-N coefficient and large-N limit');plt.legend()
save('08_lambda_N_vs_N.png')

plt.figure(figsize=(6,4))
for n in nset:
    row=[r for r in sc if int(r['N'])==n]
    plt.plot([float(r['d_NDelta'])**2 for r in row],[float(r['epsilon_cross'])/float(r['d_NDelta'])**2 for r in row],'.-',label=f'N={n}')
plt.xlabel('d² = (NΔ)²');plt.ylabel('εc / d²');plt.title('9. Collapse toward λ_N as d→0');plt.legend(ncol=2,fontsize=8)
save('09_scaling_collapse.png')

reg=read('regime_boundaries.csv');cand=[r for r in reg if r['region_status'].startswith('candidate')]
none=[r for r in reg if r['region_status'].startswith('no_')]
base=[r for r in reg if r['region_status'].startswith('certified')]
plt.figure(figsize=(7,4.5));plt.scatter([float(r['r_location']) for r in none],[float(r['phase_rad'])/np.pi for r in none],c='lightgray',s=28,label='no A/B crossing located')
plt.scatter([float(r['r_location']) for r in cand],[float(r['phase_rad'])/np.pi for r in cand],c='royalblue',s=28,label='local crossing candidate')
plt.scatter([float(r['r_location']) for r in base],[float(r['phase_rad'])/np.pi for r in base],c='red',s=85,marker='*',label='globally certified baseline')
plt.xlabel('weak-tone relative location r=(b+d)/(2d)');plt.ylabel('phase / π');plt.title('10. Targeted regime reconnaissance; other cells unclassified');plt.legend(fontsize=8)
save('10_regime_map.png')
print('saved',len(list(FIG.glob('*.png'))),'figures in',FIG)
