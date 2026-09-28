"""Extend the archived Stage-2 full-domain proof over a tiny z interval.

The archived directed-interval cell lower bounds are reused through the
1-Lipschitz property of distance to a fixed two-tone subspace. Non-objective
inner exclusions are rechecked with interval parameter boxes.
"""
import csv,json,re,hashlib
from pathlib import Path
import numpy as np
from mpmath import mp,iv
import three_to_two_tone_stage2_interval as si
from three_to_two_tone_stage4_interval_curve import validate

ROOT=Path(__file__).resolve().parent
mp.dps=70;iv.dps=55
ZLO='1.9999997';ZHI='2.0'
DZ=mp.mpf('3e-7');DE=mp.mpf('1.5e-7')
BASE=json.loads((ROOT/'stage2_reference.json').read_text())
E0=mp.mpf(BASE['epsilon_cross'])

def first_num(z):
    return mp.mpf(re.findall(r'[-+]?\d+(?:\.\d+)?(?:[Ee][-+]?\d+)?',z)[0])
def last_num(z):
    return mp.mpf(re.findall(r'[-+]?\d+(?:\.\d+)?(?:[Ee][-+]?\d+)?',z)[-1])
def ivx(z):return iv.mpf([str(z[0]),str(z[1])])

def stage2_point_margins(eta):
    outer=[r for r in csv.DictReader((ROOT/'interval_boxes.csv').open(newline=''))
           if r['label']=='crossing' and r['status']=='excluded']
    inner=[r for r in csv.DictReader((ROOT/'stage2_inner_cells.csv').open(newline=''))
           if r['label']=='crossing']
    U0=last_num(next(r for r in json.loads((ROOT/'stage2_outer_global_summary.json').read_text())
                     if r['label']=='crossing')['incumbent_upper'])+mp.mpf('1e-50')
    upper=(mp.sqrt(U0)+eta)**2
    out=[]
    for r in outer:
        L=first_num(r['interval_J_lower'])-mp.mpf('1e-50')
        lo=max(mp.mpf(0),mp.sqrt(max(L,0))-eta)**2
        out.append(lo-upper)
    obj=[]
    for r in inner:
        if r['status']=='objective_excluded':
            L=first_num(r['J_lower'])-mp.mpf('1e-40')
            lo=max(mp.mpf(0),mp.sqrt(max(L,0))-eta)**2
            obj.append(lo-upper)
    return outer,inner,U0,upper,min(out),min(obj)

def recheck_nonobjective(inner,e_box,d_box,uniform_upper,eta):
    si.D_STRONG=d_box
    counts={};fails=[];mins=[]
    for r in inner:
        typ=r['status'];counts[typ]=counts.get(typ,0)+1
        if typ=='inside_Krawczyk_box':continue
        if typ=='objective_excluded':
            L=first_num(r['J_lower'])-mp.mpf('1e-40')
            lo=max(mp.mpf(0),mp.sqrt(max(L,0))-eta)**2
            if lo>uniform_upper:continue
        X=[ivx((r['u1_lo'],r['u1_hi'])),ivx((r['u2_lo'],r['u2_hi']))]
        mid=[iv.mpf(str((mp.mpf(r['u1_lo'])+mp.mpf(r['u1_hi']))/2)),
             iv.mpf(str((mp.mpf(r['u2_lo'])+mp.mpf(r['u2_hi']))/2))]
        J0=si.objective_jet(mid[0],mid[1],e_box)
        JX=si.objective_jet(X[0],X[1],e_box)
        Y=[X[i]-mid[i] for i in range(2)]
        Jlower=J0.v+sum(J0.g[i]*Y[i] for i in range(2))
        Jlower+=sum(iv.mpf('0.5')*JX.H[i][j]*Y[i]*Y[j] for i in range(2) for j in range(2))
        if Jlower.a>iv.mpf(str(uniform_upper)).b:continue
        grad=[J0.g[i]+sum(JX.H[i][k]*Y[k] for k in range(2)) for i in range(2)]
        if any(g.a>0 or g.b<0 for g in grad):
            mins.append(min(float(abs(g).a) for g in grad if g.a>0 or g.b<0))
            continue
        try:
            C0=np.linalg.inv(np.array([[float(J0.H[i][j].mid) for j in range(2)] for i in range(2)]))
            C=[[iv.mpf(repr(float(C0[i,j]))) for j in range(2)] for i in range(2)]
            M=[[iv.mpf(int(i==j))-sum(C[i][k]*JX.H[k][j] for k in range(2))
               for j in range(2)] for i in range(2)]
            K=[mid[i]-sum(C[i][k]*J0.g[k] for k in range(2))
               +sum(M[i][k]*Y[k] for k in range(2)) for i in range(2)]
            if any(K[i].b<X[i].a or K[i].a>X[i].b for i in range(2)):continue
        except (np.linalg.LinAlgError,ZeroDivisionError,ValueError):pass
        fails.append(dict(branch=r['branch'],type=typ,depth=r['depth'],
                          u1_lo=r['u1_lo'],u2_lo=r['u2_lo']))
    return counts,fails,min(mins) if mins else None

def local_roots(e_box,d_box):
    si.D_STRONG=d_box
    result=[]
    for name in ('A','B'):
        center=BASE[f'{name}_at_crossing']
        radius=mp.mpf('0.0003' if name=='A' else '0.01')
        x0=[iv.mpf(v) for v in center]
        X=[iv.mpf([str(mp.mpf(center[i])-radius),str(mp.mpf(center[i])+radius)]) for i in range(2)]
        J0=si.objective_jet(x0[0],x0[1],e_box)
        JX=si.objective_jet(X[0],X[1],e_box)
        C0=np.linalg.inv(np.array([[float(J0.H[i][j].mid) for j in range(2)] for i in range(2)]))
        C=[[iv.mpf(repr(float(C0[i,j]))) for j in range(2)] for i in range(2)]
        M=[[iv.mpf(int(i==j))-sum(C[i][k]*JX.H[k][j] for k in range(2))
           for j in range(2)] for i in range(2)]
        Y=[X[i]-x0[i] for i in range(2)]
        K=[x0[i]-sum(C[i][k]*J0.g[k] for k in range(2))
           +sum(M[i][j]*Y[j] for j in range(2)) for i in range(2)]
        included=all(K[i].a>X[i].a and K[i].b<X[i].b for i in range(2))
        det=JX.H[0][0]*JX.H[1][1]-JX.H[0][1]*JX.H[1][0]
        result.append(dict(branch=name,included=bool(included),H11_lower=float(JX.H[0][0].a),
                           Hdet_lower=float(det.a),positive=bool(JX.H[0][0].a>0 and det.a>0)))
    return result

def main():
    q=validate(ZLO,ZHI,inflate=2)
    # The Krawczyk crossing box itself supplies an epsilon enclosure. DE is
    # deliberately larger; compare endpoints before using the perturbation.
    eL=E0-DE;eH=E0+DE
    actual= [first_num(q['X_bounds'][4][0]),last_num(q['X_bounds'][4][1])]
    covers=bool(q['verified'] and eL<actual[0] and eH>actual[1])
    roots_inside_local=True
    for name,indices,radius in (('A',(0,1),mp.mpf('0.0003')),
                                ('B',(2,3),mp.mpf('0.01'))):
        center=[mp.mpf(v) for v in BASE[f'{name}_at_crossing']]
        for i,j in enumerate(indices):
            roots_inside_local &= (first_num(q['X_bounds'][j][0])>center[i]-radius and
                                   last_num(q['X_bounds'][j][1])<center[i]+radius)
    S2=mp.mpf(21*21-1)/(12*21)
    eta_exact_computed=2*mp.sqrt(S2)*DZ+mp.sqrt(21)*DE
    # Rational overbound; 2 sqrt(440/252)<2.644 and sqrt(21)<4.583.
    eta=mp.mpf('1.5e-6')
    assert eta>eta_exact_computed
    outer,inner,U0,upper,margin_outer,margin_inner_obj=stage2_point_margins(eta)
    d_box=iv.mpf([ZLO,ZHI]);e_box=iv.mpf([str(eL),str(eH)])
    counts,fails,gradmargin=recheck_nonobjective(inner,e_box,d_box,upper,eta)
    local=local_roots(e_box,d_box)
    verified=bool(covers and roots_inside_local and margin_outer>0
                  and not fails and all(r['included'] and r['positive'] for r in local))
    out=dict(z_interval=[ZLO,ZHI],epsilon_corridor=[str(eL),str(eH)],
             epsilon_crossing_box=list(map(str,actual)),crossing_box_covered=covers,
             crossing_roots_inside_local_boxes=bool(roots_inside_local),
             parameter_distance_bound=str(eta),baseline_incumbent_upper=str(U0),
             uniform_candidate_upper=str(upper),outer_cells=len(outer),
             outer_min_margin=str(margin_outer),inner_cells=len(inner),
             inner_objective_min_margin=str(margin_inner_obj),
             inner_status_counts=counts,inner_remaining_failures=fails,
             inner_gradient_min_margin=None if gradmargin is None else str(gradmargin),
             local_roots=local,full_domain_global_verified=verified,
             source_hashes={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                            for name in ('interval_boxes.csv','stage2_inner_cells.csv',
                                         'stage2_crossing_global.json','stage2_reference.json')})
    (ROOT/'stage4_global_interval_result.json').write_text(json.dumps(out,indent=2))
    print(json.dumps(dict(outer_cells=len(outer),outer_min_margin=str(margin_outer),
                          inner_cells=len(inner),inner_remaining_failures=len(fails),
                          local_roots=local,full_domain_global_verified=verified),indent=2))

if __name__=='__main__':main()
