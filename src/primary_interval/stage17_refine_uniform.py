"""Targeted interval refinement of failed uniform N0 partition cells."""
import json
import sys
from pathlib import Path

sys.argv=['stage17_uniform_trial.py',sys.argv[1] if len(sys.argv)>1 else '10001']
import stage17_uniform_trial as s

N0=s.N0
trial=json.loads((s.ROOT/f'stage17_interval_trial_N{N0}.json').read_text())
assert trial['local']['inclusion'] and trial['local']['signs']
calc,_,_,_=s.make_core(s.local_check()[0][2])
U=s.iv.mpf(trial['trial_cost_upper'])

def classify(lo,hi):
    z=calc(lo,hi)
    if z is None:return None,None
    if z[0].a>U.b:return 'cost_excluded',str(s.mp.mpf(z[0].a-U.b))
    if s.inf.excludes(z[1]):
        m=max(s.mp.mpf(z[1].a),-s.mp.mpf(z[1].b))
        return 'gradient_excluded',str(m)
    if s.inf.excludes(z[2]):
        ga,gb=calc(lo,lo)[1],calc(hi,hi)[1]
        if ga.a>0 and gb.a>0:return 'monotone_no_root',str(min(s.mp.mpf(ga.a),s.mp.mpf(gb.a)))
        if ga.b<0 and gb.b<0:return 'monotone_no_root',str(min(-s.mp.mpf(ga.b),-s.mp.mpf(gb.b)))
    return None,None

leaves=[];unresolved=[];visited=0
def verify(lo,hi,depth,parent):
    global visited
    visited+=1
    p,m=classify(lo,hi)
    if p:
        leaves.append(dict(parent_index=parent,lo=str(lo),hi=str(hi),depth=depth,predicate=p,margin_lower=m))
        return
    if depth>=30 or visited>=100000:
        unresolved.append(dict(parent_index=parent,lo=str(lo),hi=str(hi),depth=depth))
        return
    mid=(lo+hi)/2
    verify(lo,mid,depth+1,parent)
    verify(mid,hi,depth+1,parent)

for row in trial['first_failures']:
    verify(float(row['cell'][0]),float(row['cell'][1]),0,row['index'])
out=dict(N0=N0,original_cells=len(s.leaves),original_failures=trial['failed_total'],
         refined_parents=len(trial['first_failures']),visited=visited,
         terminal_leaves=len(leaves),unresolved=unresolved,leaves=leaves,
         complete=trial['failed_total']==len(trial['first_failures']) and not unresolved and trial['root_boxes_linked'])
(s.ROOT/f'stage17_refinement_N{N0}.json').write_text(json.dumps(out,indent=2))
print(json.dumps({k:v for k,v in out.items() if k!='leaves'},indent=2))
