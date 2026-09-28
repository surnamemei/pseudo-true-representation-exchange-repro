"""Build Stage-4 tables, certificate index, and explicitly qualified plots."""
import csv,json,hashlib,re
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpmath import mp
from three_to_two_tone_stage3_full_crossings import crossing
ROOT=Path(__file__).resolve().parent
FIG=ROOT/'three_to_two_tone_stage4_figures';FIG.mkdir(exist_ok=True)
mp.dps=75

def loadcsv(p):return list(csv.DictReader((ROOT/p).open(newline='')))
def savecsv(p,rows):
    fields=list(dict.fromkeys(k for r in rows for k in r))
    with (ROOT/p).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
def nums(s):return re.findall(r'[-+]?\d+(?:\.\d+)?(?:[Ee][-+]?\d+)?',str(s))
def savefig(name):plt.tight_layout();plt.savefig(FIG/name,dpi=190);plt.close()

cont=loadcsv('continuation_branches.csv')
boxes=json.loads((ROOT/'stage4_interval_boxes.json').read_text())
g=json.loads((ROOT/'stage4_global_interval_result.json').read_text())
probe=[r for r in loadcsv('globality_margins.csv') if r['status']=='numerical_probe_only']
repair_path=ROOT/'stage4_pointwise_repaired.json'
repaired=json.loads(repair_path.read_text()) if repair_path.exists() else []
repaired_by_z={float(r['z']):r for r in repaired}
tan=json.loads((ROOT/'stage3_tangent_reference.json').read_text())
la=float(tan['lambda_N']);lc=float(tan['lambda_quartic'])

cert=[]
for q in boxes:
    m=(mp.mpf(q['z_lo'])+mp.mpf(q['z_hi']))/2
    e=crossing(mp.nstr(m,40))['epsilon_cross']
    cert.append(dict(z_lo=q['z_lo'],z_hi=q['z_hi'],z_mid=mp.nstr(m,18),
                     epsilon_mid=e,
                     epsilon_lower=nums(q['X_bounds'][4][0])[0],
                     epsilon_upper=nums(q['X_bounds'][4][1])[-1],
                     status='connected_local_crossing_interval',
                     globality_certified=False,interval_krawczyk=q['verified'],
                     slope_lower=q['slope_lower'],A_sep_lower=q['A_sep_lower'],
                     B_sep_lower=q['B_sep_lower']))
cert.append(dict(z_lo=g['z_interval'][0],z_hi=g['z_interval'][1],
                 z_mid='1.99999985',epsilon_mid=crossing('1.99999985')['epsilon_cross'],
                 epsilon_lower=g['epsilon_crossing_box'][0],
                 epsilon_upper=g['epsilon_crossing_box'][1],
                 status='connected_global_exchange_interval',
                 globality_certified=g['full_domain_global_verified'],
                 interval_krawczyk=True,slope_lower='',A_sep_lower='',B_sep_lower=''))
savecsv('certified_crossing_curve.csv',cert)

sing=[]
for r in cont:
    sing.append(dict(z_lo=r['z_NDelta'],z_hi=r['z_NDelta'],
        A_Hessian_min=r['A_hessian_min'],B_Hessian_min=r['B_hessian_min'],
        A_Hessian_det='',B_Hessian_det='',A_separation=r['A_sep_u'],
        B_separation=r['B_sep_u'],slope=r['slope_deltaJ_epsilon'],
        status='numerical_75_digit'))
for q in boxes:
    sing.append(dict(z_lo=q['z_lo'],z_hi=q['z_hi'],
        A_Hessian_min='',B_Hessian_min='',A_Hessian_det=q['A_Hdet_lower'],
        B_Hessian_det=q['B_Hdet_lower'],A_separation=q['A_sep_lower'],
        B_separation=q['B_sep_lower'],slope=q['slope_lower'],
        status='directed_interval_lower_bounds'))
savecsv('singularity_checks.csv',sing)

for r in probe:
    r['certified_outer_margin']=''
    r['certified_z_interval']=''
    r['pointwise_global_certified']=(float(r['z_NDelta'])==2.0 or
        repaired_by_z.get(float(r['z_NDelta']),{}).get('pointwise_global_verified',False))
probe.append(dict(z_NDelta='',epsilon_cross='',branch_cost='',
    third_candidate_cost='',third_candidate_margin='',
    one_tone_boundary_cost='',one_tone_margin='',confluent_boundary_cost='',
    confluent_margin='',found_distinct_minima='',grid_size='',random_starts='',
    status='certified_exclusion_bound_not_actual_third_best',
    globality_certified=True,certified_outer_margin=g['outer_min_margin'],
    certified_z_interval='[1.9999997,2]',pointwise_global_certified=''))
savecsv('globality_margins.csv',probe)

point_path=ROOT/'stage4_pointwise_global_checks.json'
point_raw=json.loads(point_path.read_text()) if point_path.exists() else []
point=repaired if repaired else point_raw
sources=['three_to_two_tone_stage4_continuation.py',
 'three_to_two_tone_stage4_interval_curve.py','three_to_two_tone_stage4_global_interval.py',
 'three_to_two_tone_stage4_globality_probe.py','three_to_two_tone_stage4_grid_global.py',
 'three_to_two_tone_stage4_pointwise_repair.py',
 'three_to_two_tone_stage2_interval.py','three_to_two_tone_stage2_outer_global.py',
 'three_to_two_tone_stage2_inner_global.py','validated_globality_certificate.json',
 'stage4_interval_boxes.json','stage4_global_interval_result.json',
 'continuation_branches.csv','globality_margins.csv']
sources.extend(p for p in ('stage4_pointwise_global_checks.json','stage4_pointwise_repaired.json')
               if (ROOT/p).exists())
certobj=dict(model='N21 symmetric 3-to-2 undamped complex exponential LS, b10 phi_pi',
    stage2_baseline_verified=json.loads((ROOT/'validated_globality_certificate.json').read_text())['verified'],
    numerical_continuation=dict(z_interval=['0.1','2'],points=len(cont),
        secant_predictor_corrector_points=sum(r['predictor_corrector_used']=='True' for r in cont),
        tangent_seed_fallbacks=sum(r['tangent_fallback']=='True' for r in cont),
        all_hessians_positive=all(float(r['A_hessian_min'])>0 and float(r['B_hessian_min'])>0 for r in cont),
        all_separations_positive=all(float(r['A_sep_u'])>0 and float(r['B_sep_u'])>0 for r in cont),
        all_slopes_positive=all(float(r['slope_deltaJ_epsilon'])>0 for r in cont),
        status='numerical_only'),
    connected_local_interval=dict(z_interval=['1.9999','2'],boxes=len(boxes),
        all_interval_verified=all(q['verified'] for q in boxes),
        min_slope_lower=min(q['slope_lower'] for q in boxes),
        min_A_sep_lower=min(q['A_sep_lower'] for q in boxes),
        min_B_sep_lower=min(q['B_sep_lower'] for q in boxes)),
    connected_global_interval=dict(z_interval=g['z_interval'],
        globally_verified=g['full_domain_global_verified'],
        full_torus_outer_excluded_cells=g['outer_cells'],
        outer_min_margin=g['outer_min_margin'],
        inner_cells_rechecked=g['inner_cells'],
        unresolved_inner_cells=len(g['inner_remaining_failures']),
        epsilon_corridor=g['epsilon_corridor']),
    pointwise_global_checks=point,
    small_z_to_baseline_global_connection_certified=False,
    largest_explicit_connected_global_interval_in_this_run=g['z_interval'],
    optimality_of_interval_width_not_claimed=True,
    source_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
    verdict='HOLD for full small-z to z=2 global connection')
(ROOT/'stage4_certificate.json').write_text(json.dumps(certobj,indent=2))

z=np.array([float(r['z_NDelta']) for r in cont]);e=np.array([float(r['epsilon_cross']) for r in cont])
pred2=la*z*z;pred4=pred2+lc*z**4
plt.figure(figsize=(8,5));plt.plot(z,e,'o-',ms=3,label='Full local-branch crossing (numerical)')
plt.plot(z,pred4,'--',label='Quadratic + quartic asymptotic')
plt.axvspan(1.9999997,2,color='green',alpha=.6,label='Certified global band (subpixel at this scale)')
plt.xlabel(r'$z=N\Delta$');plt.ylabel(r'$\epsilon_c$');plt.legend();savefig('01_epsilon_cross_vs_NDelta.png')

plt.figure(figsize=(8,5));plt.plot(z,(pred2-e)/e*100,label='Quadratic error')
plt.plot(z,(pred4-e)/e*100,label='Through quartic error')
plt.axhline(0,color='black',lw=.7);plt.axhline(1,color='gray',lw=.7,ls=':')
plt.xlabel(r'$z=N\Delta$');plt.ylabel('Signed relative error (%)');plt.legend()
savefig('02_asymptotic_vs_full_curve.png')

plt.figure(figsize=(8,5))
for key,label in [('u_A1','A satellite'),('u_A2','A central'),('u_B1','B central'),('u_B2','B satellite')]:
    plt.plot(z,[float(r[key]) for r in cont],label=label)
plt.xlabel(r'$z=N\Delta$');plt.ylabel('Fitted normalized frequency $u=N\nu$')
plt.legend(ncol=2);savefig('03_branch_frequencies.png')

plt.figure(figsize=(8,5))
for key,label in [('A_hessian_min','A'),('B_hessian_min','B')]:
    plt.semilogy(z,[float(r[key]) for r in cont],label=label)
plt.xlabel(r'$z=N\Delta$');plt.ylabel('Smallest Hessian eigenvalue (numerical)')
plt.legend();savefig('04_hessian_eigenvalues.png')

plt.figure(figsize=(8,5))
for key,label in [('A_sep_u','A'),('B_sep_u','B')]:
    plt.plot(z,[float(r[key]) for r in cont],label=label)
plt.xlabel(r'$z=N\Delta$');plt.ylabel('Fitted frequency separation in $u$')
plt.legend();savefig('05_frequency_separation.png')

plt.figure(figsize=(8,5));xp=np.array([float(r['z_NDelta']) for r in probe[:-1]])
yp=np.array([float(r['third_candidate_margin']) for r in probe[:-1]])
plt.plot(xp,yp,'o-',label='Third local-minimum candidate (numerical search)')
plt.yscale('log');plt.xlabel(r'$z=N\Delta$');plt.ylabel('Candidate cost margin over A/B')
plt.legend();savefig('06_third_best_margin.png')

fig,ax=plt.subplots(figsize=(9,5.6))
ax.plot(z,e,'k',lw=1.6,label='Numerical crossing');ax.plot(z,pred4,'--',color='tab:orange',label='Asymptotic through $z^4$')
ax.fill_between(z,0,e,color='tab:blue',alpha=.10,label='Numerical A-lower region; globality unproved')
ax.fill_between(z,e,0.28,color='tab:red',alpha=.08,label='Numerical B-lower region; globality unproved')
ax.set_xlim(0,2.02);ax.set_ylim(0,.28);ax.set_xlabel(r'$z=N\Delta$');ax.set_ylabel(r'$\epsilon$')
ax.legend(fontsize=7,loc='upper left')
ax.text(.15,.235,'Global status unknown for most $z$',fontsize=9)
inset=ax.inset_axes([.53,.16,.43,.45])
zz=np.linspace(1.9999995,2,7)
ee=np.array([float(crossing(str(v))['epsilon_cross']) for v in zz])
inset.axvspan(0,5,color='tab:purple',alpha=.10,label='local')
inset.axvspan(0,3,color='green',alpha=.25,label='global')
inset.plot((2-zz)*1e7,(float(e[-1])-ee)*1e7,'k-',lw=1.5)
inset.set_xlim(5,0);inset.set_ylim(0,1.25)
inset.set_xlabel(r'$(2-z)/10^{-7}$',fontsize=7)
inset.set_ylabel(r'$(\epsilon_c(2)-\epsilon_c)/10^{-7}$',fontsize=7)
inset.tick_params(labelsize=6)
inset.set_title('Certified bands near $z=2$',fontsize=8)
inset.legend(fontsize=6,loc='lower left')
savefig('07_full_regime_diagram.png')
print('built stage4 outputs',len(cont),len(boxes),g['z_interval'])
