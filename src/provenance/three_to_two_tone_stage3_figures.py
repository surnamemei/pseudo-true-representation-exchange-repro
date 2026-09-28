"""Six theory figures. Only interval logs, not these plots, are certificates."""
import csv,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpmath import mp
from three_to_two_tone_stage3_tangent import tangent
from three_to_two_tone_stage3_full_crossings import value_grad

ROOT=Path(__file__).resolve().parent
FIG=ROOT/'three_to_two_tone_stage3_figures';FIG.mkdir(exist_ok=True)
ref=json.loads((ROOT/'stage3_tangent_reference.json').read_text())
rows=list(csv.DictReader((ROOT/'lambdaN_theory_vs_certified.csv').open(newline='')))
co=list(csv.DictReader((ROOT/'finite_sum_coefficients.csv').open(newline='')))
d=np.array([float(r['d_NDelta']) for r in rows]);e=np.array([float(r['epsilon_cross']) for r in rows])
l0=float(ref['lambda_N']);l1=float(ref['lambda_quartic']);pred2=l0*d*d;pred4=pred2+l1*d**4
def save(n):plt.tight_layout();plt.savefig(FIG/n,dpi=180);plt.close()

plt.figure(figsize=(7,4.5));plt.plot(d,e,'ko',label='multiprecision crossings')
dd=np.linspace(0,2,250);plt.plot(dd,l0*dd**2,'--',label='quadratic');plt.plot(dd,l0*dd**2+l1*dd**4,'-',label='through quartic')
plt.scatter([2],[e[-1]],s=110,marker='*',color='red',label='Stage-2 global interval crossing')
plt.xlabel('d=NΔ');plt.ylabel('critical ε');plt.title('1. Analytic critical-amplitude law');plt.legend(fontsize=8)
save('01_analytical_vs_certified_epsilon.png')

plt.figure(figsize=(6,4));plt.plot(d*d,e/d**2,'ko',label='numerical εc/d²');plt.plot(dd*dd,l0+l1*dd*dd,'-',label='λ₀+λ₁d²')
plt.axhline(l0,ls='--',color='gray',label='λ₀');plt.xlabel('d²');plt.ylabel('εc/d²');plt.title('2. Normalized scaling collapse');plt.legend()
save('02_normalized_scaling_collapse.png')

mp.dps=45
lamgrid=np.linspace(l0-.012,l0+.012,49)
ra=[];rb=[]
for x in lamgrid:
    ll=mp.mpf(repr(float(x)))
    va=mp.findroot(lambda v:tangent(v,ll)[1],mp.mpf(ref['v_A']))
    vb=mp.findroot(lambda v:tangent(v,ll)[1],mp.mpf(ref['v_B']))
    ra.append(float(tangent(va,ll)[0]));rb.append(float(tangent(vb,ll)[0]))
plt.figure(figsize=(7,4));plt.plot(lamgrid,ra,label='R_A tangent');plt.plot(lamgrid,rb,label='R_B tangent')
plt.axvline(l0,ls='--',color='k');plt.xlabel('λ=ε/d²');plt.ylabel('J/d⁴ limit');plt.title('3. Competing finite-sum tangent costs');plt.legend()
save('03_branch_cost_asymptotics.png')

fig,axs=plt.subplots(1,2,figsize=(10,4))
axs[0].plot(d*d,(e-pred2)/d**4,'ko-',label='numerical');axs[0].axhline(l1,color='red',ls='--',label='analytic λ₁')
axs[0].set_xlabel('d²');axs[0].set_ylabel('(εc−λ₀d²)/d⁴');axs[0].legend(fontsize=8)
axs[1].plot(d*d,(e-pred4)/d**6,'ko-');axs[1].set_xlabel('d²');axs[1].set_ylabel('(εc−λ₀d²−λ₁d⁴)/d⁶')
fig.suptitle('4. Remainder diagnostics (points, not uniform bounds)');save('04_residual_after_leading_law.png')

nn=np.array([int(r['N']) for r in co]);ll=np.array([float(r['lambda_N']) for r in co])
plt.figure(figsize=(6,4));plt.plot(nn,ll,'o-');plt.xlabel('finite record N');plt.ylabel('λ_N');plt.xscale('log')
plt.title('5. Finite-sum coefficients across N');save('05_lambda_N_vs_N.png')

plt.figure(figsize=(7,4.5));plt.plot(d,e,'ko-',label='numerical local-branch continuation')
plt.plot(dd,l0*dd**2+l1*dd**4,color='steelblue',ls='--',label='small-d analytic series')
plt.scatter([2],[e[-1]],marker='*',s=180,color='red',label='global interval baseline')
plt.annotate('local theorem holds for some d₀>0;\nd₀ not numerically bounded',xy=(.16,float(e[0])),xytext=(.48,.11),
             arrowprops=dict(arrowstyle='->'),fontsize=8)
plt.annotate('open global neighborhood exists;\nradius not quantified',xy=(2,e[-1]),xytext=(1.1,.21),
             arrowprops=dict(arrowstyle='->'),fontsize=8)
plt.xlabel('d=NΔ');plt.ylabel('ε');plt.title('6. Theorem region and numerical connection');plt.legend(fontsize=7,loc='upper left')
save('06_theorem_region_map.png')
print('saved',len(list(FIG.glob('*.png'))),'figures')
