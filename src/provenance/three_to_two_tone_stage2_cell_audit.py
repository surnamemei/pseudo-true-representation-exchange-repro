"""Independent coverage audit of the saved outer interval partition."""
import csv,json,re
from decimal import Decimal,getcontext
from pathlib import Path
import numpy as np
import three_to_two_tone_branch_study as st

ROOT=Path(__file__).resolve().parent
getcontext().prec=90
EPS=Decimal('1e-12')
RADIUS=Decimal('2.00000001')
STRIP=Decimal('0.03')
NUM=re.compile(r'[-+]?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?')

def first_number(s):
    m=NUM.search(s)
    if m is None:raise ValueError(s)
    return Decimal(m.group())

def audit():
    ref=json.loads((ROOT/'stage2_reference.json').read_text())
    rows=list(csv.DictReader((ROOT/'interval_boxes.csv').open(newline='')))
    out=[]
    for rec in ref['records']:
        label=rec['label'];rs=[r for r in rows if r['label']==label]
        (ua,_),(ub,_),_=st.two_branches(21,2,10,float(rec['epsilon']),np.pi)
        centers={'near_A':tuple(Decimal(repr(float(x))) for x in ua),
                 'near_B':tuple(Decimal(repr(float(x))) for x in ub)}
        incumbent=first_number(next(r for r in rs if r['status']=='excluded')['interval_J_lower'])
        # The incumbent is sourced from the saved outer summary, not a cell.
        summ=next(z for z in json.loads((ROOT/'stage2_outer_global_summary.json').read_text()) if z['label']==label)
        incumbent=first_number(summ['incumbent_upper'])
        accepted_margin=None;exclusion_margin=None;strip_margin=None;strip_count=0
        failures=[]
        for r in rs:
            mlo,mhi,hlo,hhi=[Decimal(r[k]) for k in ('m_lo','m_hi','h_lo','h_hi')]
            if r['status']=='excluded':
                if r['verified']!='True':failures.append(('unverified_exclusion',r))
                margin=first_number(r['interval_J_lower'])-incumbent
                exclusion_margin=margin if exclusion_margin is None else min(exclusion_margin,margin)
                if hlo- EPS<=STRIP:
                    strip_count+=1
                    strip_margin=margin if strip_margin is None else min(strip_margin,margin)
            else:
                if r['status'] not in centers:failures.append(('unresolved',r));continue
                c1,c2=centers[r['status']]
                deviations=[]
                for m in (mlo-EPS,mhi+EPS):
                    for h in (hlo-EPS,hhi+EPS):
                        deviations.extend((abs(Decimal(21)*(m-h)-c1),abs(Decimal(21)*(m+h)-c2)))
                margin=RADIUS-max(deviations)
                accepted_margin=margin if accepted_margin is None else min(accepted_margin,margin)
                if margin<=0:failures.append(('outside_inner_box',r))
                if hlo-EPS<=STRIP:failures.append(('accepted_intersects_strip',r))
        out.append(dict(label=label,total_cells=len(rs),excluded_cells=sum(r['status']=='excluded' for r in rs),
                        accepted_cells=sum(r['status']!='excluded' for r in rs),
                        inner_box_centers_u={k:[str(v[0]),str(v[1])] for k,v in centers.items()},
                        accepted_box_min_margin_u=str(accepted_margin),
                        excluded_J_min_margin=str(exclusion_margin),
                        coalescent_h_max=str(STRIP),coalescent_separation_u_max=str(Decimal(42)*STRIP),
                        coalescent_intersected_cells=strip_count,
                        coalescent_J_min_margin=str(strip_margin),
                        failures=len(failures),complete=len(failures)==0))
    (ROOT/'stage2_cell_audit.json').write_text(json.dumps(out,indent=2))
    print(json.dumps(out,indent=2))

if __name__=='__main__':audit()
