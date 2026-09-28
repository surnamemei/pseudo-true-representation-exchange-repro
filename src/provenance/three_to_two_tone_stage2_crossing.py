"""Interval bracketing of the transverse A/B crossing."""

import csv
import json
from pathlib import Path

import mpmath as mp
from mpmath import iv

from three_to_two_tone_stage2_krawczyk import check,ref

ROOT=Path(__file__).resolve().parent
iv.dps=90;mp.mp.dps=110


def D(u):
    return iv.mpf(1)+sum(2*iv.cos(iv.mpf(k)*u/21) for k in range(1,11))


def dJ_de(u1,u2,e):
    g=D(u2-u1);den=iv.mpf(441)-g*g
    p1=D(u1+2)+D(u1-2);p2=D(u2+2)+D(u2-2)
    z1=D(u1-10);z2=D(u2-10)
    h1=p1-e*z1;h2=p2-e*z2
    norm_deriv=iv.mpf(0)
    for t in range(-10,11):
        s=iv.mpf(t)/21
        norm_deriv+=2*e-4*iv.cos(2*s)*iv.cos(10*s)
    cap_deriv=-2*21*(h1*z1+h2*z2)+2*g*(z1*h2+h1*z2)
    return norm_deriv-cap_deriv/den


ec=mp.mpf(ref['epsilon_cross'])
rad=mp.mpf('1e-15')
eL=mp.nstr(ec-rad,105);eU=mp.nstr(ec+rad,105)
rows=[]
for name,e in [('lower',eL),('upper',eU)]:
    A,JA=check('crossing',e,'A',return_taylor=True)
    B,JB=check('crossing',e,'B',return_taylor=True)
    difference=JA-JB
    rows.append(dict(endpoint=name,epsilon=e,DeltaJ_interval=str(difference),
                     sign='negative' if difference.b<0 else ('positive' if difference.a>0 else 'undetermined'),
                     A_krawczyk=A['krawczyk_inclusion'],B_krawczyk=B['krawczyk_inclusion'],
                     A_H_positive=A['positive_definite'],B_H_positive=B['positive_definite'],
                     A_objective=A['J_root_enclosure'],B_objective=B['J_root_enclosure']))

ei=iv.mpf([eL,eU])
A,_=check('crossing',ei,'A',return_taylor=True)
B,_=check('crossing',ei,'B',return_taylor=True)
uA=[iv.mpf(v)+iv.mpf(['-1e-10','1e-10']) for v in ref['A_at_crossing']]
uB=[iv.mpf(v)+iv.mpf(['-1e-10','1e-10']) for v in ref['B_at_crossing']]
slope=dJ_de(*uA,ei)-dJ_de(*uB,ei)
output=dict(precision_dps=iv.dps,epsilon_lower=eL,epsilon_upper=eU,
            endpoints=rows,
            A_krawczyk_entire_interval=A['krawczyk_inclusion'],
            B_krawczyk_entire_interval=B['krawczyk_inclusion'],
            A_H_positive_entire_interval=A['positive_definite'],
            B_H_positive_entire_interval=B['positive_definite'],
            DeltaJ_derivative_interval=str(slope),
            transverse_positive=bool(slope.a>0),
            bracket_certified=bool(rows[0]['sign']=='negative' and rows[1]['sign']=='positive'
                                   and slope.a>0 and A['krawczyk_inclusion'] and B['krawczyk_inclusion']))
(ROOT/'stage2_crossing_interval.json').write_text(json.dumps(output,indent=2),encoding='utf-8')
with (ROOT/'crossing_certificate.csv').open('w',newline='',encoding='utf-8') as f:
    writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
print('bracket',eL,eU,flush=True)
print('signs',rows[0]['sign'],rows[1]['sign'],'slope',slope,flush=True)
print('certified',output['bracket_certified'],flush=True)
