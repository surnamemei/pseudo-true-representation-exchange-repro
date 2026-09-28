"""Tight pointwise root/inner repair after a completed outer interval partition."""
import json,sys
from pathlib import Path
from mpmath import iv,mp
import three_to_two_tone_stage2_interval as si
from three_to_two_tone_stage2_inner_global import process
from three_to_two_tone_stage4_interval_curve import validate,root
from three_to_two_tone_stage4_grid_global import local_box

ROOT=Path(__file__).resolve().parent
iv.dps=55;mp.dps=75

def repair(z,outer_cells,outer_failures,outer_min_margin):
    zs=str(z);si.D_STRONG=iv.mpf(zs)
    proof=validate(zs,zs,base_radius='1e-20')
    a1,a2,b1,b2,ec=root(mp.mpf(zs))
    er=[mp.nstr(ec-mp.mpf('1e-19'),65),mp.nstr(ec+mp.mpf('1e-19'),65)]
    e=iv.mpf(er)
    JA=si.objective_jet(iv.mpf(mp.nstr(a1,65)),iv.mpf(mp.nstr(a2,65)),e)
    JB=si.objective_jet(iv.mpf(mp.nstr(b1,65)),iv.mpf(mp.nstr(b2,65)),e)
    U=min(JA.v.b,JB.v.b)
    inner=[]
    radius='1e-6' if mp.mpf(zs)<=mp.mpf('.75') else '1e-5'
    accepted_radius='9e-7' if mp.mpf(zs)<=mp.mpf('.75') else '9e-6'
    for lab,center in [('A',(a1,a2)),('B',(b1,b2))]:
        done,n,rows,remain=process('crossing',lab,er,tuple(map(float,center)),U,
            max_cells=30000,known_center_override=[mp.nstr(v,65) for v in center],
            local_radius_override=accepted_radius)
        loc=local_box(center,e,radius)
        inner.append(dict(branch=lab,complete=done,visited=n,remaining=remain,local=loc))
        print(zs,lab,done,n,loc,flush=True)
    # Expanding epsilon from 1e-50 to 1e-19 changes ||x|| by <5e-19.
    # Squared residual changes by <1e-17; 1e-16 is a safe rational overbound.
    outer_ok=bool(outer_failures==0 and outer_min_margin>1e-16)
    verified=bool(proof['verified'] and outer_ok and
                  all(v['complete'] and v['local']['inclusion'] and v['local']['positive']
                      for v in inner))
    return dict(z=zs,outer_cells=outer_cells,outer_failures=outer_failures,
                outer_min_margin=outer_min_margin,root_verified=proof['verified'],
                epsilon_bracket=er,outer_transfer_loss_bound='1e-16',
                outer_transfer_valid=outer_ok,inner=inner,
                pointwise_global_verified=verified,
                note='isolated z node only')

if __name__=='__main__':
    z=sys.argv[1];nc=int(sys.argv[2]);fail=int(sys.argv[3]);margin=float(sys.argv[4])
    q=repair(z,nc,fail,margin)
    fn=ROOT/f'stage4_bridge_pointwise_{z.replace(".","p")}.json'
    fn.write_text(json.dumps(q,indent=2))
    print(json.dumps(dict(z=z,verified=q['pointwise_global_verified'],
                          inner=[(r['branch'],r['complete'],r['visited']) for r in q['inner']]),indent=2))
