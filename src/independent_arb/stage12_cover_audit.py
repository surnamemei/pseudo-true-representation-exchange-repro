import csv,json
from decimal import Decimal
from pathlib import Path
from collections import defaultdict
from stage12_arb_backend import I,arb
R=Path(__file__).resolve().parent

def cover(rows,keys):
 rects=[tuple(Decimal(r[k]) for k in keys) for r in rows]
 ys=sorted(set(y for a,b,c,d in rects for y in (c,d)));ix={v:i for i,v in enumerate(ys)};n=len(ys)-1
 lo=min(a for a,b,c,d in rects);hi=max(b for a,b,c,d in rects);ev=defaultdict(list)
 for a,b,c,d in rects:ev[a].append((ix[c],ix[d]-1,1));ev[b].append((ix[c],ix[d]-1,-1))
 mn=[0]*(4*n+8);lazy=[0]*(4*n+8)
 def add(i,l,r,a,b,v):
  if a<=l and r<=b:mn[i]+=v;lazy[i]+=v;return
  m=(l+r)//2
  if a<=m:add(2*i,l,m,a,b,v)
  if b>m:add(2*i+1,m+1,r,a,b,v)
  mn[i]=lazy[i]+min(mn[2*i],mn[2*i+1])
 missing=[];mincover=10**9
 for x in sorted(ev):
  for a,b,v in ev[x]:add(1,0,n-1,a,b,v)
  if x<hi:
   mincover=min(mincover,mn[1])
   if mn[1]<1:missing.append(str(x))
 return dict(cells=len(rows),domain=list(map(str,[lo,hi,ys[0],ys[-1]])),minimum_cover_count=mincover,uncovered_slabs=missing,complete=not missing)

out=[]
audit=json.loads((R/'stage2_cell_audit.json').read_text())
c31=json.loads((R/'N31_global_certificate.json').read_text())['high_precision_numerical_root']
for N in (21,31):
 outer=list(csv.DictReader(open('interval_boxes.csv' if N==21 else 'N31_interval_cells.csv')))
 inner=list(csv.DictReader(open('stage2_inner_cells.csv' if N==21 else 'N31_inner_interval_cells.csv')))
 for label in sorted(set(r.get('label','crossing') for r in outer)):
  rr=cover([r for r in outer if r.get('label','crossing')==label],('m_lo','m_hi','h_lo','h_hi'))
  a,b,c,d=[I(v) for v in rr['domain']];pi=I(arb.pi());inflate=I('1e-12')
  rr.update(N=N,case=label,kind='outer',true_torus_covered=bool((a-inflate).hi < -pi.hi and (b+inflate).lo>pi.hi and c.hi<=0 and (d+inflate).lo>(pi/2).hi));out.append(rr)
  for branch in ('A','B'):
   rr=cover([r for r in inner if r.get('label','crossing')==label and r['branch']==branch],('u1_lo','u1_hi','u2_lo','u2_hi'))
   center=(next(a for a in audit if a['label']==label)['inner_box_centers_u']['near_'+branch] if N==21 else [repr(float(c31['u_'+branch+'1'])),repr(float(c31['u_'+branch+'2']))])
   ca,cb=map(Decimal,center);rad=Decimal('2.00000001')
   rr['expected_neighborhood_matches']=list(map(Decimal,rr['domain']))==[ca-rad,ca+rad,cb-rad,cb+rad]
   rr.update(N=N,case=label,kind='inner_'+branch);out.append(rr)
Path('stage12_replay/coverage_audit.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2));assert all(r['complete'] and r.get('true_torus_covered',True) and r.get('expected_neighborhood_matches',True) for r in out)
refined=[]
for N in (21,31):
 rows=list(csv.DictReader((R/f'stage12_replay/N{N}_refinement_leaves.csv').open(newline='')))
 groups=defaultdict(list)
 for row in rows:groups[(row['case'],row['branch'],row['parent'])].append(row)
 for (case,branch,parent),rs in groups.items():
  q=cover(rs,('u1_lo','u1_hi','u2_lo','u2_hi'))
  q.update(N=N,case=case,branch=branch,parent=parent,parent_domain_matches=list(map(Decimal,q['domain']))==list(map(Decimal,parent.split(','))))
  assert q['complete'] and q['parent_domain_matches'];refined.append(q)
Path('stage12_replay/refinement_coverage_audit.json').write_text(json.dumps(refined,indent=2))
print('Verified refinement parent covers:',len(refined))
