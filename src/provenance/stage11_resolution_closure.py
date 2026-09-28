"""Numerical one-tone/confluent closure controls at the seven reported points."""
import csv
from pathlib import Path
import numpy as np
import three_to_two_tone_branch_study as st
R=Path(__file__).resolve().parent
rows=list(csv.DictReader(open(R/'resolution_sweep.csv')))
for r in rows:
    x=st.x_record(21,float(r['z']),10,float(r['epsilon_cross']),np.pi)
    one,confl=st.boundary_cost(x,21)
    r.update(one_tone_minimum=one,confluent_minimum=confl,confluent_gap=confl-float(r['J_cross']))
    print(r['z'],r['confluent_gap'])
with (R/'resolution_sweep.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
