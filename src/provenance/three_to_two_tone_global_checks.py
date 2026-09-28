"""Additional independent full-circle searches around representative crossings."""

import csv
from pathlib import Path

import numpy as np

import three_to_two_tone_branch_study as st
import three_to_two_tone_explore as ex

ROOT = Path(__file__).resolve().parent
rows = []


def check(label, n, d, b, phase, eps, ratio=1.0, grid_size=512, random_count=60):
    st.set_n(n)
    x = ex.atom(-d/n) + ratio*ex.atom(d/n) + eps*np.exp(1j*phase)*ex.atom(b/n)
    ua, ja = st.refine([-d,d], x, n)
    ub, jb = st.refine([0,b], x, n)
    _, _, mins, grid_seed_count = st.independent_minima(
        x, n, grid_size=grid_size, random_count=random_count,
        seed=271828)
    sa, sb = st.stats(ua,x,n), st.stats(ub,x,n)
    best_ab = min(ja,jb)
    third = mins[2]['J'] if len(mins)>2 else np.nan
    row = dict(label=label,N=n,d_NDelta=d,b_Nomega3=b,phase_rad=phase,
               epsilon=eps,strong_amp_ratio=ratio,grid_size=grid_size,
               random_starts=random_count,grid_seed_count=grid_seed_count,
               minima_found=len(mins),J_A=ja,J_B=jb,J_best_search=mins[0]['J'],
               J_second_search=mins[1]['J'],J_third_search=third,
               best_search_minus_best_AB=mins[0]['J']-best_ab,
               top2_match_AB=int(all(any(abs(m['J']-j)<1e-8 and
                   np.linalg.norm(m['u']-u)<1e-3 for u,j in ((ua,ja),(ub,jb)))
                   for m in mins[:2])),
               branch_frequency_distance=float(np.linalg.norm(ua-ub)),
               branch_min_hess=min(sa['hess_min'],sb['hess_min']),
               branch_max_gram_cond=max(sa['gram_cond'],sb['gram_cond']))
    rows.append(row)
    print(row,flush=True)


ec = st.crossing(21,2,10,np.pi)[0]
check('primary_fine_grid',21,2,10,np.pi,ec,grid_size=2048,random_count=250)
for ph,b in ((np.pi-.3,9.0),(np.pi-.3,11.0),
             (np.pi+.3,9.0),(np.pi+.3,11.0)):
    z = st.crossing(21,2,b,ph,lo=.1,hi=.55)
    check('phase_location_corner',21,2,b,ph,z[0])
for ratio in (.95,1.05):
    def diff(eps):
        st.set_n(21)
        x=ex.atom(-2/21)+ratio*ex.atom(2/21)-eps*ex.atom(10/21)
        return st.refine([-2,2],x,21)[1]-st.refine([0,10],x,21)[1]
    from scipy.optimize import brentq
    root=brentq(diff,.2,.3)
    check('amplitude_perturbation',21,2,10,np.pi,root,ratio=ratio)

with (ROOT/'three_to_two_tone_global_checks.csv').open('w',newline='',encoding='utf-8') as f:
    writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
