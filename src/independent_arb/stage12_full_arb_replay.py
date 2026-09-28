from stage12_arb_backend import *
from functools import lru_cache
import stage12_arb_backend as backend
backend.energy=lru_cache(None)(backend.energy)
energy=backend.energy
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'stage12_replay';OUT.mkdir(exist_ok=True)
REFINEMENT_LEAVES=[]

def root_check(N,center,e,r='1e-10'):
 x0=[I(v) for v in center];X=[box(v,r) for v in center]
 j0=objective(*x0,e,N);jx=objective(*X,e,N);Y=[X[i]-x0[i] for i in range(2)]
 K,M=krawczyk(x0,Y,j0,jx)
 det=jx.H[0][0]*jx.H[1][1]-jx.H[0][1]*jx.H[1][0]
 gram=N*N-d0(X[1]-X[0],N)**2
 J=taylor(j0,jx,Y)
 contraction=max(sum(M[i][k].absmax() for k in range(2)) for i in range(2))
 return dict(X=X,J=J,K=K,included=all(K[i].lo>X[i].lo and K[i].hi<X[i].hi for i in range(2)),contraction=contraction,H11=jx.H[0][0],det=det,gram=gram)
def krawczyk(mid,Y,j0,jx):
 h=[[j0.H[i][j].mid() for j in range(2)] for i in range(2)]
 d=h[0][0]*h[1][1]-h[0][1]*h[1][0]
 C=[[I((h[1][1]/d).mid()),I((-h[0][1]/d).mid())],[I((-h[1][0]/d).mid()),I((h[0][0]/d).mid())]]
 M=[[I(int(i==j))-sum(C[i][k]*jx.H[k][j] for k in range(2)) for j in range(2)] for i in range(2)]
 K=[mid[i]-sum(C[i][k]*j0.g[k] for k in range(2))+sum(M[i][k]*Y[k] for k in range(2)) for i in range(2)]
 return K,M
def taylor(j0,jx,Y):
 return j0.v+sum(j0.g[i]*Y[i] for i in range(2))+sum(I('0.5')*jx.H[i][k]*Y[i]*Y[k] for i in range(2) for k in range(2))

def center_projected(m,h,e,E,N):
 tmax=(N-1)//2;h1=I(2)-e;h2=I(0);c1sq=I(1);c2sq=I(0);w1sq=I(0);w2sq=I(0)
 for t in range(1,tmax+1):
  z=h*t;co,si=z.cos(),z.sin();sc=si/z;sp=(z*co-si)/(z*z)
  p=(I(2)*t/N-m*t).cos()+(-I(2)*t/N-m*t).cos()-e*(I(10)*t/N-m*t).cos()
  q=(I(2)*t/N-m*t).sin()+(-I(2)*t/N-m*t).sin()-e*(I(10)*t/N-m*t).sin()
  h1+=2*co*p;h2+=2*t*sc*q;c1sq+=2*co**2;c2sq+=2*t*t*sc**2;w1sq+=2*t*t*si**2;w2sq+=2*t**4*sp**2
 return E-h1**2/c1sq-h2**2/c2sq,c1sq.sqrt(),c2sq.sqrt(),w1sq.sqrt(),w2sq.sqrt()
def outer_bound(row,e,N):
 a,b,c,d=[Decimal(row[k]) for k in ('m_lo','m_hi','h_lo','h_hi')];inflate=Decimal('1e-12')
 m=I(str((a+b)/2));h=I(str((c+d)/2));rm=I(str((b-a)/2+inflate));rh=I(str((d-c)/2+inflate))
 E=energy(e,N);Jc,c1,c2,w1,w2=center_projected(m,h,e,E,N);tm=(N-1)//2
 t2=I(sum(t**4 for t in range(-tm,tm+1))).sqrt();t3=I(sum(t**6 for t in range(-tm,tm+1))).sqrt()
 d1=w1*rh+t2*rh**2/2;d2=w2*rh+t3*rh**2/6
 if not(c1.lo>d1.hi and c2.lo>d2.hi and Jc.lo>0):return I(0),I(0)
 theta=I(2).sqrt()*tm*rm+((d1/(c1-d1))**2+(d2/(c2-d2))**2).sqrt()
 residual=Jc.sqrt()-E.sqrt()*theta
 return (residual**2 if residual.lo>0 else I(0)),theta

def inner(row,e,N,U,center):
 ds=[Decimal(row[k]) for k in ('u1_lo','u1_hi','u2_lo','u2_hi')];a,b,c,d=ds
 if row['status']=='inside_Krawczyk_box':
  r=Decimal('0.0003' if row['branch']=='A' else '0.01');c1,c2=map(Decimal,center)
  return 'inside_Krawczyk_box',a>=c1-r and b<=c1+r and c>=c2-r and d<=c2+r,''
 X=[I(str(a),str(b)),I(str(c),str(d))];mid=[I(str((a+b)/2)),I(str((c+d)/2))]
 j0=objective(*mid,e,N);jx=objective(*X,e,N);Y=[X[i]-mid[i] for i in range(2)]
 L=taylor(j0,jx,Y);g=[j0.g[i]+sum(jx.H[i][k]*Y[k] for k in range(2)) for i in range(2)]
 proofs={'objective_excluded':bool(L.lo>U),'gradient_excluded':any(q.lo>0 or q.hi<0 for q in g)}
 if row['status']=='krawczyk_excluded' or not any(proofs.values()):
  K,_=krawczyk(mid,Y,j0,jx);proofs['krawczyk_excluded']=any(K[i].hi<X[i].lo or K[i].lo>X[i].hi for i in range(2))
 preferred=row['status'];actual=preferred if proofs.get(preferred,False) else next((k for k,v in proofs.items() if v),'NONE')
 return actual,actual!='NONE',bounds(L)

def refined_inner(row,e,N,U,center):
 actual,ok,val=inner(row,e,N,U,center)
 if ok:return actual,ok,val,0
 stack=[(row,0)];leaves=0;visited=0
 while stack:
  cell,depth=stack.pop();visited+=1
  status,passed,_=inner(cell,e,N,U,center)
  if passed:
   leaves+=1
   REFINEMENT_LEAVES.append(dict(N=N,case=row.get('label','crossing'),branch=row['branch'],parent=','.join(row[k] for k in ('u1_lo','u1_hi','u2_lo','u2_hi')),u1_lo=cell['u1_lo'],u1_hi=cell['u1_hi'],u2_lo=cell['u2_lo'],u2_hi=cell['u2_hi'],depth=depth,predicate=status))
   continue
  if depth>=12:return 'NONE',False,val,visited
  a,b,c,d=[Decimal(cell[k]) for k in ('u1_lo','u1_hi','u2_lo','u2_hi')]
  lk,hk=('u1_lo','u1_hi') if b-a>=d-c else ('u2_lo','u2_hi')
  mid=(Decimal(cell[lk])+Decimal(cell[hk]))/2
  left=dict(cell);right=dict(cell);left[hk]=str(mid);right[lk]=str(mid)
  stack.extend([(left,depth+1),(right,depth+1)])
 return 'Arb_refined_cover',True,val,visited

def run(N):
 started=time.time();checks=[];fails=[];allrows=[]
 def add(case,kind,branch,idx,ok,value='',original='',actual=''):
  r=dict(N=N,case=case,kind=kind,branch=branch,index=idx,agreement='PASS' if ok else 'FAIL',original_predicate=original,arb_predicate=actual,arb_value=value)
  allrows.append(r)
  if not ok:fails.append(r);print('FAIL',r,flush=True)
 c21=json.loads((ROOT/'validated_globality_certificate.json').read_text());c31=json.loads((ROOT/'N31_global_certificate.json').read_text())
 ref=json.loads((ROOT/'stage2_reference.json').read_text());audit=json.loads((ROOT/'stage2_cell_audit.json').read_text())
 if N==21:
  epsbr=[c21['crossing']['epsilon_lower'],c21['crossing']['epsilon_upper']]
  cases={r['label']:dict(e=(I(*epsbr) if r['label']=='crossing' else I(r['epsilon']))+box(0,'1e-95'),centers={k:r[k] for k in ('A','B')},floatcenters={k:next(a for a in audit if a['label']==r['label'])['inner_box_centers_u']['near_'+k] for k in ('A','B')}) for r in ref['records']}
 else:
  epsbr=c31['epsilon_bracket'];roots=c31['high_precision_numerical_root'];centers={k:[roots['u_'+k+'1'],roots['u_'+k+'2']] for k in ('A','B')}
  cases={'crossing':dict(e=I(*epsbr)+box(0,'1e-95'),centers=centers,floatcenters={k:[repr(float(v)) for v in centers[k]] for k in centers})}
 for label,c in cases.items():
  rr={}
  for branch in ('A','B'):
   rr[branch]=root_check(N,c['centers'][branch],c['e'])
   for radius,rrr in [('1e-10',rr[branch]),('0.0003' if branch=='A' else '0.01',root_check(N,c['centers'][branch],c['e'],'0.0003' if branch=='A' else '0.01'))]:
    for kind,ok,v in [('root_inclusion',rrr['included'],';'.join(bounds(v) for v in rrr['K'])),('contraction',rrr['contraction']<1,str(rrr['contraction'])),('H11',rrr['H11'].lo>0,bounds(rrr['H11'])),('Hessian_det',rrr['det'].lo>0,bounds(rrr['det'])),('Gram',rrr['gram'].lo>0,bounds(rrr['gram']))]:add(label,kind,branch,radius,bool(ok),v)
  c['U']=min(q['J'].hi for q in rr.values());c['roots']=rr
  # Fresh incumbent is a feasible POINT fit, which is valid even without a root proof.
  c['U']=min(c['U'],min(objective(*[I(v) for v in c['centers'][k]],c['e'],N).v.hi for k in ('A','B')))
 cross=cases['crossing'];r=cross['roots'];s=slope_branch(*r['A']['X'],cross['e'],N)-slope_branch(*r['B']['X'],cross['e'],N)
 add('crossing','cost_gap_slope','A-B','',bool(s.lo>0),bounds(s))
 for label,e in zip(('lower','upper'),epsbr):
  endpoint=I(e)+box(0,'1e-95')
  ra=root_check(N,cross['centers']['A'],endpoint);rb=root_check(N,cross['centers']['B'],endpoint);diff=ra['J']-rb['J']
  add(label,'endpoint_sign','A-B','',bool(diff.hi<0 if label=='lower' else diff.lo>0),bounds(diff))
 outer=list(csv.DictReader((ROOT/('interval_boxes.csv' if N==21 else 'N31_interval_cells.csv')).open(newline='')))
 counts={};margins={}
 for i,row in enumerate(outer):
  label=row.get('label','crossing');c=cases[label];status=row['status'];val=''
  if status=='excluded':
   bound,theta=outer_bound(row,c['e'],N);ok=bool(bound.lo>c['U']);val=bounds(bound)+'; U='+str(c['U'])+'; theta='+bounds(theta)
   margin=float(bound.lo-c['U']);margins[label]=min(margins.get(label,float('inf')),margin)
  else:
   a,b,cc,d=[Decimal(row[k]) for k in ('m_lo','m_hi','h_lo','h_hi')];infl=Decimal('1e-12');u1=(N*(a-d-2*infl),N*(b-cc+2*infl));u2=(N*(a+cc-2*infl),N*(b+d+2*infl));rad=Decimal('2.00000001')
   ok=any(u1[0]>=Decimal(v[0])-rad and u1[1]<=Decimal(v[0])+rad and u2[0]>=Decimal(v[1])-rad and u2[1]<=Decimal(v[1])+rad for v in c['floatcenters'].values())
  add(label,'outer', '',i,ok,val,status,status);counts['outer_'+status]=counts.get('outer_'+status,0)+1
  if i%2500==0:print('N',N,'outer',i,'/',len(outer),'elapsed',round(time.time()-started),flush=True)
 innerrows=list(csv.DictReader((ROOT/('stage2_inner_cells.csv' if N==21 else 'N31_inner_interval_cells.csv')).open(newline='')))
 refinements=0;refinement_nodes=0
 for i,row in enumerate(innerrows):
  label=row.get('label','crossing');c=cases[label];actual,ok,val,nref=refined_inner(row,c['e'],N,c['U'],c['centers'][row['branch']])
  refinements+=int(nref>0);refinement_nodes+=nref
  add(label,'inner',row['branch'],i,ok,val,row['status'],actual);counts['inner_'+row['status']]=counts.get('inner_'+row['status'],0)+1
  if i%200==0:print('N',N,'inner',i,'/',len(innerrows),flush=True)
 with (OUT/f'N{N}_predicates.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=allrows[0]);w.writeheader();w.writerows(allrows)
 with (OUT/f'N{N}_refinement_leaves.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=REFINEMENT_LEAVES[0]);w.writeheader();w.writerows(REFINEMENT_LEAVES)
 result=dict(N=N,backend='python-flint/Arb',python_flint=flint.__version__,FLINT=flint.__FLINT_VERSION__,precision_bits=ctx.prec,predicates=len(allrows),passed=len(allrows)-len(fails),failures=fails,counts=counts,outer_min_margin=margins,elapsed_seconds=time.time()-started,refined_archived_inner_cells=refinements,refinement_nodes=refinement_nodes,alternate_inner_proofs=sum(r['kind']=='inner' and r['original_predicate']!=r['arb_predicate'] and r['agreement']=='PASS' for r in allrows))
 (OUT/f'N{N}_summary.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__':run(int(sys.argv[1]))
