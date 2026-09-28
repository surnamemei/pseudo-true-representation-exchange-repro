"""Check that the two-global-branch result persists over the root bracket."""
import csv,json
from pathlib import Path
from mpmath import iv
import numpy as np
import three_to_two_tone_branch_study as st
from three_to_two_tone_stage2_krawczyk import check
from three_to_two_tone_stage2_inner_global import process
from three_to_two_tone_stage2_outer_global import verify_cell,data_energy

ROOT=Path(__file__).resolve().parent
ref=json.loads((ROOT/'stage2_reference.json').read_text())
cross=json.loads((ROOT/'stage2_crossing_interval.json').read_text())
er=[cross['epsilon_lower'],cross['epsilon_upper']]

def main():
    ja,tja=check('crossing',er,'A',return_taylor=True)
    jb,tjb=check('crossing',er,'B',return_taylor=True)
    best_upper=min(tja.b,tjb.b)
    iv.dps=65
    e=iv.mpf(er)
    E=data_energy(e)
    rows=[r for r in csv.DictReader((ROOT/'interval_boxes.csv').open(newline=''))
          if r['label']=='crossing' and r['status']=='excluded']
    failed=[];min_margin=None
    for r in rows:
        z=(int(r['depth']),float(r['m_lo']),float(r['m_hi']),
           float(r['h_lo']),float(r['h_hi']),0.,0.,'excluded')
        ok,bb,why=verify_cell(z,e,E,best_upper)
        if not ok:failed.append((r,why))
        elif bb is not None:
            margin=bb.a-best_upper
            min_margin=margin if min_margin is None else min(min_margin,margin)
    (ua,_),(ub,_),_=st.two_branches(21,2,10,float(ref['epsilon_cross']),np.pi)
    iv.dps=45
    inner=[]
    for which,center in [('A',ua),('B',ub)]:
        complete,visited,cells,remaining=process('crossing',which,er,center,best_upper,max_cells=10000)
        inner.append(dict(branch=which,complete=complete,visited=visited,
                          terminal=len(cells),remaining=remaining))
    out=dict(epsilon_bracket=er,outer_excluded_cells=len(rows),outer_failed=len(failed),
             outer_min_margin=str(min_margin),inner=inner,
             local_inclusion=[ja['krawczyk_inclusion'],jb['krawczyk_inclusion']],
             full_bracket_global=not failed and all(x['complete'] for x in inner)
                and ja['krawczyk_inclusion'] and jb['krawczyk_inclusion'])
    (ROOT/'stage2_crossing_global.json').write_text(json.dumps(out,indent=2,default=str))
    print(json.dumps(out,indent=2,default=str),flush=True)

if __name__=='__main__':main()
