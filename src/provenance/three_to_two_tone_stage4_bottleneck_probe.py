"""Sparse bridge-point probe; records every terminal unresolved outer cell.

These float64 lower bounds are diagnostic only. A later directed-interval pass
must validate any exclusion used in a theorem.
"""
import csv,gzip,json,io,contextlib
from pathlib import Path
import numpy as np
from mpmath import mp
from three_to_two_tone_stage3_full_crossings import crossing
from three_to_two_tone_stage1_global import adaptive,evaluate
import three_to_two_tone_branch_study as st

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'stage4_bottleneck_cells';OUT.mkdir(exist_ok=True)
N=21

def D(v):return 1+sum(2*np.cos(k*v/N) for k in range(1,11))

def one(z,max_active=200000,maxdepth=32):
    q=crossing(str(z));e=float(q['epsilon_cross'])
    A=np.array([float(q['u_A1']),float(q['u_A2'])]);B=np.array([float(q['u_B1']),float(q['u_B2'])])
    U=float(q['J_A']);x=st.x_record(N,z,10,e,np.pi)
    log=io.StringIO()
    with contextlib.redirect_stdout(log):
        done,summary,_,un=adaptive(x,(A,B),U+1e-8,N,radius_u=2.,
                                  maxdepth=maxdepth,max_active=max_active,write_rows=False)
    (OUT/f'z_{str(z).replace(".","p")}_adaptive.log').write_text(log.getvalue())
    if len(un):
        j,lb=evaluate(un,x,N)
        m=(un[:,0]+un[:,1])/2;h=(un[:,2]+un[:,3])/2
        u1=N*(m-h);u2=N*(m+h);sep=u2-u1
        gram=np.maximum(0,N*N-np.array([D(v)**2 for v in sep]))
        width_m=N*(un[:,1]-un[:,0]);width_h=N*(un[:,3]-un[:,2])
        distA=np.hypot(u1-A[0],u2-A[1]);distB=np.hypot(u1-B[0],u2-B[1])
        rows=[]
        for i in range(len(un)):
            rows.append((i,*un[i],u1[i],u2[i],j[i],lb[i],lb[i]-U,
                         gram[i],width_m[i],width_h[i],sep[i],
                         distA[i],distB[i]))
        header=['cell_id','m_lo','m_hi','h_lo','h_hi','u1_center','u2_center',
                'J_center','J_lower_float','lower_gap_float','gram_det_center',
                'width_u_mean','width_u_halfsep','coalescent_sep_u',
                'distance_A_u','distance_B_u']
        with gzip.open(OUT/f'z_{str(z).replace(".","p")}_unresolved.csv.gz','wt',newline='') as f:
            w=csv.writer(f);w.writerow(header);w.writerows(rows)
        result=dict(z=z,complete=done,unresolved_cells=len(un),depth=summary[-1][0],
            min_gram_det=float(gram.min()),min_sep_u=float(sep.min()),
            median_sep_u=float(np.median(sep)),
            max_center_gap=float((j-U).max()),min_center_gap=float((j-U).min()),
            max_lower_gap=float((lb-U).max()),min_lower_gap=float((lb-U).min()),
            min_distance_to_A=float(distA.min()),min_distance_to_B=float(distB.min()),
            proportion_near_coalescent_sep_lt_0p1=float(np.mean(sep<.1)),
            proportion_near_A_or_B_2p5=float(np.mean(np.minimum(distA,distB)<2.5)),
            widths_u_mean=[float(width_m.min()),float(np.median(width_m)),float(width_m.max())],
            widths_u_halfsep=[float(width_h.min()),float(np.median(width_h)),float(width_h.max())],
            crossing_epsilon=e,branch_cost=U,adaptive_last=summary[-1],
            status='float_diagnostic_unresolved')
    else:
        result=dict(z=z,complete=done,unresolved_cells=0,depth=summary[-1][0],
                    crossing_epsilon=e,branch_cost=U,adaptive_last=summary[-1],
                    status='float_partition_complete_interval_validation_pending')
    return result

if __name__=='__main__':
    vals=(.25,.5,.75,1.0,1.25,1.5,1.75,2.0)
    out=[]
    for z in vals:
        q=one(z);out.append(q)
        print(z,q['complete'],q['unresolved_cells'],q['depth'],flush=True)
    (ROOT/'stage4_bridge_outer_probe.json').write_text(json.dumps(out,indent=2))
