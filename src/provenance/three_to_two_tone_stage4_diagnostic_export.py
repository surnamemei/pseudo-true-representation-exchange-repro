"""Join every outer-unresolved cell with the targeted bound results."""
import csv,gzip,json,math
from decimal import Decimal
from pathlib import Path
import numpy as np
from three_to_two_tone_stage3_full_crossings import crossing

ROOT=Path(__file__).resolve().parent;OUT=ROOT/'stage4_bottleneck_cells'

def readgz(path):
    with gzip.open(path,'rt',newline='') as f:return list(csv.DictReader(f))

def zrows(z):
    p=f'z_{str(z).replace(".","p")}_'
    raw=readgz(OUT/(p+'unresolved.csv.gz'))
    taylor=readgz(OUT/(p+'taylor.csv.gz'))
    assert len(raw)==len(taylor)
    opt={int(r['cell_id']):r for r in readgz(OUT/(p+'optimized.csv.gz'))} if z==.25 else {}
    centered={int(r['cell_id']):r for r in readgz(OUT/(p+'centered_interval.csv.gz'))} if z==.25 else {}
    interval={int(r['cell_id']):r for r in readgz(OUT/(p+'interval_cells.csv.gz'))} if z==.5 else {}
    U=float(crossing(str(z))['J_A'])
    res=[]
    for a,t in zip(raw,taylor):
        cid=int(a['cell_id']);s=float(a['coalescent_sep_u'])
        if z==.25 and cid in centered:
            rigorous=centered[cid]['rigorous_lower'];gap=centered[cid]['rigorous_gap']
            status='directed_interval_objective_excluded';mechanism='centered_quadratic_plus_third_order_remainder'
        elif z==.5:
            q=interval[cid];rigorous=q['Jlower'];gap=q['margin']
            status=q['status'];mechanism='centered_second_order_interval_taylor'
        else:
            rigorous='';gap='';status='unresolved_directed_interval'
            mechanism=('float_centered_Taylor_excludes; directed_check_pending'
                       if t['Taylor_objective_excluded_float']=='True'
                       else 'float_gradient_sign_excludes; directed_check_pending')
        diagnosis=('C_numerical_no_stationary_rigorous_cost_exclusion' if z==.25 and cid in centered
                   else 'C_rigorous_cost_exclusion_stationarity_unchecked' if z==.5
                   else 'D_near_confluent_unverified' if s<.5
                   else 'C_suspected_enclosure_looseness_unverified')
        res.append(dict(z=z,cell_id=cid,
            m_lo=a['m_lo'],m_hi=a['m_hi'],h_lo=a['h_lo'],h_hi=a['h_hi'],
            u1_center=a['u1_center'],u2_center=a['u2_center'],
            J_center_float=a['J_center'],generic_lower_float=a['J_lower_float'],
            generic_gap_float=a['lower_gap_float'],
            Taylor_lower_float=t['Taylor_lower_float'],Taylor_gap_float=t['Taylor_gap_float'],
            rigorous_lower=rigorous,rigorous_gap=gap,rigorous_status=status,
            gap_to_AB_center_float=float(a['J_center'])-U,
            gram_det_center_float=a['gram_det_center'],
            gram_condition_center_float=(21+math.sqrt(max(0.,441-float(a['gram_det_center']))))/
                                       (21-math.sqrt(max(0.,441-float(a['gram_det_center'])))),
            coalescent_sep_u=s,distance_to_sep_0p1_strip_u=max(0.,s-.1),
            width_u_mean=a['width_u_mean'],width_u_halfsep=a['width_u_halfsep'],
            width_u1=float(a['width_u_mean'])+float(a['width_u_halfsep']),
            width_u2=float(a['width_u_mean'])+float(a['width_u_halfsep']),
            distance_A_u=a['distance_A_u'],distance_B_u=a['distance_B_u'],
            float_gradient_excludes=t['gradient_excluded_float'],
            enclosure_source=mechanism,diagnosis=diagnosis,
            numerical_box_min_gap=opt[cid]['minimum_gap'] if cid in opt else '',
            stationary_inside_numerical=opt[cid]['stationary_inside'] if cid in opt else ''))
    return res

def transfer_limit():
    arr=[r for r in csv.DictReader((ROOT/'interval_boxes.csv').open())
         if r['label']=='crossing' and r['status']=='excluded']
    q=min(arr,key=lambda r:Decimal(r['interval_J_lower'].strip('[]').split(',')[0]))
    m=(float(q['m_lo'])+float(q['m_hi']))/2
    h=(float(q['h_lo'])+float(q['h_hi']))/2
    u1=21*(m-h);u2=21*(m+h)
    g=1+sum(2*np.cos(k*(u2-u1)/21) for k in range(1,11))
    Jlo=Decimal(q['interval_J_lower'].strip('[]').split(',')[0])
    U=Decimal('0.8118078719175344435916832959181420930472402187285955555174407801224145')
    return dict(z_baseline=2.,cell=q,u1_center=u1,u2_center=u2,
                gram_det_center=441-g*g,
                gram_condition_center=(21+abs(g))/(21-abs(g)),coalescent_sep_u=u2-u1,
                distance_to_sep_0p1_strip_u=u2-u1-.1,
                width_m=float(q['m_hi'])-float(q['m_lo']),
                width_h=float(q['h_hi'])-float(q['h_lo']),
                J_lower_baseline=str(Jlo),J_lower_gap_baseline=str(Jlo-U),
                source='uniform z perturbation bound loses more than the small baseline cell margin')

def transfer_thresholds():
    U=0.8118078719175344436
    rows=[]
    for i,r in enumerate(csv.DictReader((ROOT/'interval_boxes.csv').open())):
        if r['label']!='crossing' or r['status']!='excluded':continue
        L=float(r['interval_J_lower'].strip('[]').split(',')[0])
        m=(float(r['m_lo'])+float(r['m_hi']))/2
        h=(float(r['h_lo'])+float(r['h_hi']))/2
        u1=21*(m-h);u2=21*(m+h)
        g=1+sum(2*np.cos(k*(u2-u1)/21) for k in range(1,11))
        # Eta=5(2-z) is the conservative Stage-4 transfer bound.
        dzcritical=(math.sqrt(L)-math.sqrt(U))/10
        rows.append(dict(source_row=i,m_lo=r['m_lo'],m_hi=r['m_hi'],
                         h_lo=r['h_lo'],h_hi=r['h_hi'],
                         u1_center=u1,u2_center=u2,gram_det_center=441-g*g,
                         gram_condition_center=(21+abs(g))/(21-abs(g)),
                         coalescent_sep_u=u2-u1,
                         distance_to_sep_0p1_strip_u=max(0.,u2-u1-.1),
                         width_m=float(r['m_hi'])-float(r['m_lo']),
                         width_h=float(r['h_hi'])-float(r['h_lo']),
                         J_lower_baseline=L,gap_baseline=L-U,
                         critical_dz_under_eta_5dz=dzcritical))
    with gzip.open(ROOT/'stage4_transfer_cell_thresholds.csv.gz','wt',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    return rows

if __name__=='__main__':
    rows=zrows(.25)+zrows(.5)
    with gzip.open(ROOT/'stage4_failure_cells.csv.gz','wt',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (ROOT/'stage4_bridge_transfer_limit.json').write_text(json.dumps(transfer_limit(),indent=2))
    transfer_thresholds()
    print('exported',len(rows),'cells')
