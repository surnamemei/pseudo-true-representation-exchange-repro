"""Annotate predicate-wise tested floors; no minimality inference is made."""
import csv
from pathlib import Path

p=Path(__file__).with_name('largeN_condition_thresholds.csv')
with p.open(newline='',encoding='utf-8') as f:
    rows=list(csv.DictReader(f))
mapping={
    'C1/C4 joint root inclusion and equality':(5001,'fails at 1501; passes at 5001'),
    'C1 tangent Gram factors':(1501,'passes at 1501'),
    'C2 curvature':(1501,'passes at 1501'),
    'C3 satellite coefficient':(1501,'passes at 1501'),
    'C5 transversality':(1501,'passes at 1501'),
    'C6 regular cost':(10001,'full partition fails at 5001; passes after targeted refinement at 10001'),
    'C6 derivative sign':(10001,'full C6 partition fails at 5001; passes after targeted refinement at 10001'),
    'C6 monotone cell':(10001,'full C6 partition fails at 5001; passes after targeted refinement at 10001'),
    'C7 coalescent':(1501,'passes at 1501'),
    'C8 finite-period tail':(1501,'passes at 1501'),
}
assert {r['condition'] for r in rows}==set(mapping)
for row in rows:
    n,evidence=mapping[row['condition']]
    row['smallest_tested_predicate_sufficient_odd_floor']=n
    row['lower_trial_evidence']=evidence
fields=list(rows[0])
with p.open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=fields)
    w.writeheader();w.writerows(rows)
