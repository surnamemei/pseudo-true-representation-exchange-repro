"""Numerical full-torus probes only; never promote grid searches to proof."""
import csv
from pathlib import Path
import numpy as np
import three_to_two_tone_branch_study as st
from three_to_two_tone_stage3_full_crossings import crossing

ROOT=Path(__file__).resolve().parent

def main():
    rows=[]
    for z in (0.2,0.5,1.0,1.5,2.0):
        q=crossing(str(z));e=float(q['epsilon_cross'])
        ua=np.array([float(q['u_A1']),float(q['u_A2'])])
        ub=np.array([float(q['u_B1']),float(q['u_B2'])])
        x=st.x_record(21,z,10,e,np.pi)
        _,_,minima,seed_count=st.independent_minima(x,21,grid_size=128,random_count=24,seed=51)
        third=[m for m in minima if min(np.linalg.norm(m['u']-ua),np.linalg.norm(m['u']-ub))>0.025]
        best=float(q['J_A'])
        one,conf=st.boundary_cost(x,21)
        thirdj=third[0]['J'] if third else float('nan')
        rows.append(dict(z_NDelta=z,epsilon_cross=e,branch_cost=best,
                         third_candidate_cost=thirdj,
                         third_candidate_margin=thirdj-best,
                         one_tone_boundary_cost=one,one_tone_margin=one-best,
                         confluent_boundary_cost=conf,confluent_margin=conf-best,
                         found_distinct_minima=len(minima),grid_size=128,random_starts=24,
                         status='numerical_probe_only',
                         globality_certified=(z==2.0)))
        print(z,len(minima),thirdj-best,conf-best,flush=True)
    with (ROOT/'globality_margins.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

if __name__=='__main__':main()
