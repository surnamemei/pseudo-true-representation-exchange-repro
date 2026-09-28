"""Tight exterior exclusion around the tracked optima (aggregated logs)."""

import csv
from pathlib import Path

import numpy as np

import three_to_two_tone_branch_study as st
from three_to_two_tone_stage1_global import adaptive

root=Path(__file__).resolve().parent
n,d,b,phase=21,2,10,np.pi
ec=st.crossing(n,d,b,phase)[0]
rows=[]
for label,eps in [('below',ec-.035),('crossing',ec),('above',ec+.035)]:
    (ua,ja),(ub,jb),x=st.two_branches(n,d,b,eps,phase)
    threshold=min(ja,jb)-1e-4
    done,summary,_,un=adaptive(x,(ua,ub),threshold,n,radius_u=.1,
                               maxdepth=42,max_active=1200000,write_rows=False)
    for depth,cells,excluded,near,unresolved,min_lb,max_lb in summary:
        rows.append(dict(epsilon_label=label,epsilon=eps,threshold=threshold,
                         u_box_radius=.1,depth=depth,cells=cells,
                         excluded=excluded,inside_branch_boxes=near,
                         unresolved=unresolved,min_lower_bound=min_lb,
                         max_lower_bound=max_lb,complete=int(done)))
    print(label,'complete',done,'max_depth',summary[-1][0],
          'final_unresolved',len(un),flush=True)
with (root/'globality_tight_summary.csv').open('w',newline='',encoding='utf-8') as f:
    writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
