"""Certify adjacent epsilon slabs share the same A/B stationary roots."""
import json
import sys
from decimal import Decimal
from pathlib import Path

from mpmath import iv,mp
from three_to_two_tone_stage3_full_crossings import value_grad
import stage7_n31_certificate as lc
import stage17_parametric_jet as p
import three_to_two_tone_stage2_interval as it

ROOT=Path(__file__).resolve().parent
mp.dps=90;iv.dps=85
it.N=21;lc.interval.N=21;lc.interval.objective_jet=p.objective
suffix=sys.argv[1] if len(sys.argv)>1 else ''
root=json.loads((ROOT/f'stage17_epsilon_root_tiling{suffix}.json').read_text())
rows=root['records'];assert root['all_pass'] and len(rows)==root['slabs']
def endpoint(s):return Decimal(s.strip('[]').split(',')[0].strip())
def D(z):
    return iv.mpf(1)+sum(2*iv.cos(iv.mpf(k)*z/21) for k in range(1,11))
small=Decimal('1e-20')
gram_min=None
for row in rows:
    for branch in ('A','B'):
        X=row['branches'][branch]['X']
        a=iv.mpf([str(endpoint(X[0][0])-small),str(endpoint(X[0][1])+small)])
        b=iv.mpf([str(endpoint(X[1][0])-small),str(endpoint(X[1][1])+small)])
        gram=iv.mpf(441)-D(b-a)**2
        assert gram.a>0
        v=mp.mpf(gram.a)
        gram_min=v if gram_min is None else min(v,gram_min)
joints=[];minimum_containment=None
for i in range(len(rows)-1):
    left,right=rows[i],rows[i+1]
    e=left['epsilon_interval'][1]
    assert e==right['epsilon_interval'][0]
    for branch in ('A','B'):
        seed=tuple(mp.mpf(v) for v in left['branches'][branch]['center'])
        root_point=mp.findroot(lambda u,v:value_grad(u,v,mp.mpf(e),2)[1:],
                               seed,tol=mp.mpf('1e-65'),maxsteps=40)
        check,_,Xsmall=lc.local_check([str(v) for v in root_point],[e,e],'1e-8')
        assert check['inclusion'] and check['positive_Hessian']
        gaps=[]
        for side in (left,right):
            Xlarge=side['branches'][branch]['X']
            for k in range(2):
                low=Decimal(str(mp.mpf(Xsmall[k].a)))-endpoint(Xlarge[k][0])
                high=endpoint(Xlarge[k][1])-Decimal(str(mp.mpf(Xsmall[k].b)))
                gaps.extend((low,high))
        mg=min(gaps)
        assert mg>small,(i,branch,mg)
        minimum_containment=mg if minimum_containment is None else min(mg,minimum_containment)
        joints.append(dict(left_slab=i,right_slab=i+1,branch=branch,
                           epsilon=e,center=[str(v) for v in root_point],
                           tiny_inclusion=True,min_box_containment=str(mg)))
    if (i+1)%32==0:print('joints',i+1,flush=True)
out=dict(slabs=len(rows),joint_boundaries=len(rows)-1,checks=len(joints),
         all_joints_pass=True,minimum_box_containment=str(minimum_containment),
         minimum_gram_determinant_lower=str(gram_min),records=joints)
(ROOT/f'stage17_epsilon_joints{suffix}.json').write_text(json.dumps(out,indent=2))
print(json.dumps({k:v for k,v in out.items() if k!='records'},indent=2))
