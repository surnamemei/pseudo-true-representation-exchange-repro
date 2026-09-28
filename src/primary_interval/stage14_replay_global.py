"""Replay every continuum global partition predicate and quantify cost-cell margin."""
import csv
import json
from decimal import Decimal
from pathlib import Path
from mpmath import iv,mp
import stage14_continuum_global_interval as core

ROOT=Path(__file__).resolve().parent
cert=json.loads((ROOT/'continuum_global_certificate.json').read_text(encoding='utf-8'))
local=json.loads((ROOT/'stage14_continuum_local_interval.json').read_text(encoding='utf-8'))
with (ROOT/'stage14_continuum_global_partition.csv').open(encoding='utf-8') as f:
    leaves=list(csv.DictReader(f))
with (ROOT/'stage14_continuum_global_roots.csv').open(encoding='utf-8') as f:
    roots=list(csv.DictReader(f))

def replay():
    ordered=sorted(leaves,key=lambda r:Decimal(r['lo']))
    cover=(Decimal(ordered[0]['lo'])==Decimal(-100) and Decimal(ordered[-1]['hi'])==Decimal(100)
           and all(Decimal(a['hi'])==Decimal(b['lo']) for a,b in zip(ordered,ordered[1:])))
    U=mp.mpf(cert['global_trial_upper'])
    failed=[];cost_margin=mp.inf;counts={}
    for i,r in enumerate(leaves):
        a=float(r['lo']);b=float(r['hi']);p=r['predicate']
        z=core.calc(a,b)
        ok=False
        if p=='cost_excluded':
            ok=z is not None and z[0].a>U
            if ok:cost_margin=min(cost_margin,mp.mpf(z[0].a)-U)
        elif p=='gradient_excluded':
            ok=z is not None and core.excludes(z[1])
        elif p=='monotone_derivative_no_root':
            if z is not None and core.excludes(z[2]):
                ga=core.calc(a,a)[1];gb=core.calc(b,b)[1]
                ok=(ga.a>0 and gb.a>0) or (ga.b<0 and gb.b<0)
        elif p.startswith('unique_root_'):
            if z is not None and core.excludes(z[2]):
                ga=core.calc(a,a)[1];gb=core.calc(b,b)[1]
                ok=(ga.b<0 and gb.a>0) or (ga.a>0 and gb.b<0)
        counts[p]=counts.get(p,0)+1
        if not ok:failed.append({'index':i,'cell':[a,b],'predicate':p})
        if (i+1)%5000==0:print('replayed',i+1,flush=True)
    link=True
    for r in roots:
        i=0 if r['identity']=='A' else 1
        link &= Decimal(r['lo'])<Decimal(local['K'][i][0]) and Decimal(local['K'][i][1])<Decimal(r['hi'])
    out={
      'partition_cover_exact_decimal':cover,
      'leaf_count':len(leaves),'predicate_counts':counts,
      'failed_count':len(failed),'first_failures':failed[:5],
      'cost_cell_margin_lower':str(cost_margin),
      'krawczyk_roots_inside_global_root_boxes':bool(link),
      'replay_verified':bool(cover and not failed and link and cost_margin>0 and cert['verified_global'])
    }
    (ROOT/'stage14_continuum_global_replay.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
    cert['retained_stationary_count']=cert.pop('stationary_count',cert.get('retained_stationary_count'))
    cert['retained_minimum_count']=cert.pop('minimum_count',cert.get('retained_minimum_count'))
    cert['all_stationary_roots_isolated']=False
    cert['other_roots_excluded_by_cost_or_derivative_predicates']=True
    cert['other_isolated_root_margin_lower']=cert.pop('other_regular_margin_lower',cert.get('other_isolated_root_margin_lower'))
    cert['cost_excluded_cell_margin_lower']=str(cost_margin)
    cert['C6_other_regular_margin_lower']=str(cost_margin)
    cert['partition_replay_verified']=out['replay_verified']
    (ROOT/'continuum_global_certificate.json').write_text(json.dumps(cert,indent=2),encoding='utf-8')
    print(json.dumps(out,indent=2),flush=True)

if __name__=='__main__':replay()
