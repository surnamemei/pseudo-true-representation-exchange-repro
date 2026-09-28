"""Tight local boxes for archived Stage-4 pointwise outer exclusions."""
import json
from pathlib import Path
import numpy as np
from mpmath import iv,mp
import three_to_two_tone_stage2_interval as si
from three_to_two_tone_stage2_inner_global import process
from three_to_two_tone_stage4_interval_curve import validate,root
from three_to_two_tone_stage4_grid_global import local_box

ROOT=Path(__file__).resolve().parent
iv.dps=55;mp.dps=75

def repair(q):
    z=q['z'];si.D_STRONG=iv.mpf(z)
    proof=validate(z,z,base_radius='1e-48')
    a1,a2,b1,b2,ec=root(mp.mpf(z))
    er=[mp.nstr(ec-mp.mpf('1e-47'),65),mp.nstr(ec+mp.mpf('1e-47'),65)]
    e=iv.mpf(er)
    JA=si.objective_jet(iv.mpf(mp.nstr(a1,65)),iv.mpf(mp.nstr(a2,65)),e)
    JB=si.objective_jet(iv.mpf(mp.nstr(b1,65)),iv.mpf(mp.nstr(b2,65)),e)
    U=min(JA.v.b,JB.v.b)
    inner=[]
    for lab,center in [('A',(a1,a2)),('B',(b1,b2))]:
        done,n,rows,remain=process('crossing',lab,er,tuple(map(float,center)),U,
            max_cells=30000,known_center_override=[mp.nstr(v,65) for v in center],
            local_radius_override='9e-6')
        loc=local_box(center,e,'1e-5')
        inner.append(dict(branch=lab,complete=done,visited=n,remaining=remain,local=loc))
        print(z,lab,done,n,loc,flush=True)
    # The original outer interval was narrower in epsilon. Transfer by the
    # 1-Lipschitz residual bound: data perturbation <5e-47. Here ||x||<11,
    # so the lower-minus-upper squared-cost gap changes by <1e-40.
    outer_transfer=bool(q['outer_proposal_complete'] and q['outer_failures']==0
        and q['outer_min_margin'] is not None and q['outer_min_margin']>1e-40)
    done=bool(proof['verified'] and outer_transfer and
        all(t['complete'] and t['local']['inclusion'] and t['local']['positive'] for t in inner))
    return dict(z=z,root_krawczyk_tight=proof['verified'],
                outer_transfer_valid=outer_transfer,inner=inner,
                pointwise_global_verified=done,
                outer_transfer_gap_loss_bound='1e-40',
                note='isolated z node only; does not certify z intervals')

if __name__=='__main__':
    raw=json.loads((ROOT/'stage4_pointwise_global_checks.json').read_text())
    out=[]
    for q in raw:out.append(repair(q))
    (ROOT/'stage4_pointwise_repaired.json').write_text(json.dumps(out,indent=2))
