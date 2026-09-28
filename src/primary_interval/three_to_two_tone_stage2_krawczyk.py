"""Krawczyk enclosures for the two baseline stationary points.

This proves only local existence/uniqueness inside stated tiny boxes.
It does not by itself prove global optimality.
"""

import csv
import json
from pathlib import Path

import numpy as np
from mpmath import iv

import three_to_two_tone_branch_study as st
from three_to_two_tone_stage2_interval import objective_jet,interval_str

iv.dps=90
ROOT=Path(__file__).resolve().parent
ref=json.loads((ROOT/'stage2_reference.json').read_text(encoding='utf-8'))


def sup_abs(x):
    return abs(x).b


def check(label,e_value,which,radius='1e-10',return_taylor=False):
    center=ref['A_at_crossing'] if which=='A' else ref['B_at_crossing']
    # For epsilon +/-0.035 use the corresponding high-precision branch centers.
    if label in ('below','above'):
        center=next(r[which] for r in ref['records'] if r['label']==label)
    e=iv.mpf(e_value)+iv.mpf(['-1e-95','1e-95'])
    x0=[iv.mpf(center[0]),iv.mpf(center[1])]
    rr=iv.mpf(radius)
    X=[x0[0]+iv.mpf([-rr.b,rr.b]),x0[1]+iv.mpf([-rr.b,rr.b])]
    j0=objective_jet(x0[0],x0[1],e)
    jX=objective_jet(X[0],X[1],e)
    H=jX.H
    # A floating approximation to the inverse is allowed in Krawczyk;
    # all subsequent products use interval arithmetic.
    hf=np.array([[float(j0.H[i][j].mid) for j in range(2)] for i in range(2)])
    C=np.linalg.inv(hf)
    cc=[[iv.mpf(repr(float(C[i,j]))) for j in range(2)] for i in range(2)]
    M=[[iv.mpf(int(i==j))-sum(cc[i][k]*H[k][j] for k in range(2))
        for j in range(2)] for i in range(2)]
    Y=[X[i]-x0[i] for i in range(2)]
    K=[x0[i]-sum(cc[i][k]*j0.g[k] for k in range(2))
       +sum(M[i][k]*Y[k] for k in range(2)) for i in range(2)]
    included=all(K[i].a>X[i].a and K[i].b<X[i].b for i in range(2))
    contraction=max(sum(sup_abs(M[i][j]) for j in range(2)) for i in range(2))
    strict=H[0][0].a>0 and (H[0][0]*H[1][1]-H[0][1]*H[1][0]).a>0
    # Taylor enclosure for the objective at the unique root.
    taylor=j0.v
    for i in range(2):
        taylor+=j0.g[i]*Y[i]
    for i in range(2):
        for k in range(2):
            taylor+=iv.mpf('0.5')*H[i][k]*Y[i]*Y[k]
    output=dict(label=label,branch=which,epsilon=str(e_value),radius_u=radius,
                x0=center,X1=interval_str(X[0]),X2=interval_str(X[1]),
                K1=interval_str(K[0]),K2=interval_str(K[1]),
                krawczyk_inclusion=bool(included),
                contraction_upper=str(contraction.b),
                contraction_lt_one=bool(contraction<1),
                H11_lower=str(H[0][0].a),
                determinant_lower=str((H[0][0]*H[1][1]-H[0][1]*H[1][0]).a),
                positive_definite=bool(strict),
                J_root_enclosure=interval_str(taylor),
                J_lower=str(taylor.a),J_upper=str(taylor.b))
    return (output,taylor) if return_taylor else output


if __name__=='__main__':
    records=[]
    for r in ref['records']:
        for which in ('A','B'):
            z=check(r['label'],r['epsilon'],which)
            records.append(z)
            print(r['label'],which,z['krawczyk_inclusion'],
                  z['contraction_lt_one'],z['positive_definite'],
                  z['H11_lower'][:25],z['determinant_lower'][:25],flush=True)
    (ROOT/'stage2_local_interval_checks.json').write_text(json.dumps(records,indent=2),encoding='utf-8')
