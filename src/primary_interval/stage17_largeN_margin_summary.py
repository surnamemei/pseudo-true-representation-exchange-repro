"""Final predicate margins for the N0=10001 uniform cover with four refinements."""
import json
import sys
from pathlib import Path
sys.argv=['stage17_uniform_trial.py','10001']
import stage17_uniform_trial as s

trial=json.loads((s.ROOT/'stage17_interval_trial_N10001.json').read_text())
refined=json.loads((s.ROOT/'stage17_refinement_N10001.json').read_text())
assert refined['complete']
failed={v['index'] for v in trial['first_failures']}
calc,_,_,_=s.make_core(s.local_check()[0][2])
U=s.iv.mpf(trial['trial_cost_upper'])
mins={}
def update(p,m,index):
    v=s.mp.mpf(m.a)
    if p not in mins or v<s.mp.mpf(mins[p]['margin_lower']):
        mins[p]=dict(margin_lower=str(v),index=index)
for i,row in enumerate(s.leaves):
    if i in failed:continue
    lo,hi=float(row['lo']),float(row['hi'])
    z=calc(lo,hi);p=row['predicate']
    if p=='cost_excluded':
        m=z[0].a-U.b
    elif p=='gradient_excluded':
        m=z[1].a if z[1].a>0 else -z[1].b
    elif p=='monotone_derivative_no_root':
        ga,gb=calc(lo,lo)[1],calc(hi,hi)[1]
        cur=z[2].a if z[2].a>0 else -z[2].b
        aa=ga.a if ga.a>0 else -ga.b
        bb=gb.a if gb.a>0 else -gb.b
        m=min(cur.a,aa.a,bb.a)
    else:
        m=z[2].a
    assert m>0,(i,p,m)
    update(p,m,i)
    if (i+1)%5000==0:print('checked',i+1,flush=True)
for j,row in enumerate(refined['leaves']):
    update(row['predicate'],s.iv.mpf(row['margin_lower']),f"refined_{j}")
out=dict(N0=10001,original_cells=len(s.leaves),replaced_parents=len(failed),
         refined_leaves=len(refined['leaves']),predicate_margins=mins,
         local=trial['local'],coalescent_margin_lower=trial['coalescent_margin_lower'],
         tail_margin_lower=trial['tail_margin_lower'])
(s.ROOT/'stage17_largeN_margins.json').write_text(json.dumps(out,indent=2))
print(json.dumps({k:v for k,v in out.items() if k!='local'},indent=2))
