"""Replay and record the two N=31 inner interval subdivisions."""
import csv
import json
from pathlib import Path

import numpy as np
from mpmath import iv

import three_to_two_tone_stage2_interval as interval
from three_to_two_tone_stage2_inner_global import process
from stage7_n31_certificate import local_check

ROOT=Path(__file__).resolve().parent
cert=json.loads((ROOT/'N31_global_certificate.json').read_text(encoding='utf-8'))
root=cert['high_precision_numerical_root']
centers=[(root['u_A1'],root['u_A2']),(root['u_B1'],root['u_B2'])]
eps=cert['epsilon_bracket']
interval.N=31
iv.dps=90
incumbent=min(local_check(centers[0],eps)[1].b,
              local_check(centers[1],eps)[1].b)
iv.dps=45
rows=[]
for which,center in zip(('A','B'),centers):
    result=process('crossing',which,eps,np.array([float(v) for v in center]),
                   incumbent,max_cells=50000,known_center_override=center,
                   local_radius_override='0.000299999999' if which=='A' else '0.009999999')
    valid,visited,cells,remaining=result
    print(which,valid,visited,remaining,flush=True)
    if not valid or remaining:
        raise RuntimeError('inner replay unresolved')
    rows.extend(cells)
with (ROOT/'N31_inner_interval_cells.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
