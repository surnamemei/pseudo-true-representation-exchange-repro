"""Independent exact-decimal volume and containment audit of Stage-17 covers."""
import csv
import json
import re
import sys
from collections import defaultdict
from decimal import Decimal,getcontext
from pathlib import Path

R=Path(__file__).resolve().parent
getcontext().prec=130
radius=Decimal(sys.argv[1]) if len(sys.argv)>1 else Decimal('0.0025')
suffix='' if radius==Decimal('0.0025') else '_'+str(radius).replace('.','p')
outer_suffix=str(radius).replace('.','p')
outer=json.loads((R/f'stage17_epsilon_outer_refinement_{outer_suffix}_14287.json').read_text())
inner=json.loads((R/f'stage17_epsilon_inner_tiled{suffix}.json').read_text())
root=json.loads((R/f'stage17_epsilon_root_tiling{suffix}.json').read_text())
joint=json.loads((R/f'stage17_epsilon_joints{suffix}.json').read_text())
assert outer['complete_outer'] and inner['complete'] and root['all_pass'] and joint['all_joints_pass']

def dec(s):return Decimal(str(s))
def area(cell,names):
    return (dec(cell[names[1]])-dec(cell[names[0]]))*(dec(cell[names[3]])-dec(cell[names[2]]))
def inside(cell,parent,names):
    return (dec(cell[names[0]])>=dec(parent[names[0]]) and
            dec(cell[names[1]])<=dec(parent[names[1]]) and
            dec(cell[names[2]])>=dec(parent[names[2]]) and
            dec(cell[names[3]])<=dec(parent[names[3]]))

outer_names=('m_lo','m_hi','h_lo','h_hi')
src=[r for r in csv.DictReader((R/'interval_boxes.csv').open()) if r['label']=='crossing' and r['status']=='excluded']
def key(r):return dec(re.findall(r'[-+]?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?',r['interval_J_lower'])[0])
src.sort(key=key)
assert len(src)==outer['parents']==14287
outer_sum=defaultdict(Decimal);outer_count=defaultdict(int)
for leaf in outer['refined_leaves']:
    k=leaf['parent_index'];c=leaf['cell'];assert 0<=k<len(src)
    assert inside(c,src[k],outer_names)
    outer_sum[k]+=area(c,outer_names)
    outer_count[k]+=1
assert len(outer_sum)==len(src)
for k,parent in enumerate(src):assert outer_sum[k]==area(parent,outer_names),(k,outer_sum[k])

inner_names=('u1_lo','u1_hi','u2_lo','u2_hi')
orig=[r for r in csv.DictReader((R/'stage2_inner_cells.csv').open()) if r['label']=='crossing']
assert len(orig)==819
inner_sum=defaultdict(Decimal);inner_count=defaultdict(int)
for leaf in inner['records']:
    k=leaf['parent_index'];c=leaf['cell'];i,j=leaf['epsilon_slabs']
    assert 0<=k<len(orig) and 0<=i<j<=root['slabs']
    assert inside(c,orig[k],inner_names)
    inner_sum[k]+=area(c,inner_names)*(j-i)
    inner_count[k]+=1
assert len(inner_sum)==len(orig)
for k,parent in enumerate(orig):
    assert inner_sum[k]==area(parent,inner_names)*root['slabs'],(k,inner_sum[k])

for i in range(root['slabs']-1):
    assert root['records'][i]['epsilon_interval'][1]==root['records'][i+1]['epsilon_interval'][0]
final=[str(Decimal('0.248190630172278')-radius),
       str(Decimal('0.248190630172277')+radius)]
outer_interval=[dec(v) for v in outer['epsilon_interval']]
root_interval=[dec(v) for v in root['interval']]
assert outer_interval[0]<dec(final[0])<dec(final[1])<outer_interval[1]
assert root_interval[0]<dec(final[0])<dec(final[1])<root_interval[1]
bracket=json.loads((R/'stage2_crossing_global.json').read_text())['epsilon_bracket']
assert dec(final[0])<dec(bracket[0])<dec(bracket[1])<dec(final[1])
out=dict(final_closed_interval=final,final_width=str(dec(final[1])-dec(final[0])),
         outer_parent_cells=len(src),outer_terminal_leaves=len(outer['refined_leaves']),
         outer_refined_parent_count=sum(v>1 for v in outer_count.values()),
         inner_parent_cells=len(orig),inner_terminal_leaves=len(inner['records']),
         inner_unique_parent_count=len(inner_sum),epsilon_slabs=root['slabs'],
         root_joints=joint['checks'],all_exact_cover_identities_pass=True)
(R/f'stage17_epsilon_coverage_audit{suffix}.json').write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
