"""Subdivide only the 81 archived inner cells failing the wider z transfer."""
import csv,json
from decimal import Decimal,getcontext
from pathlib import Path
import numpy as np
from mpmath import iv,mp
import three_to_two_tone_stage2_interval as si
from three_to_two_tone_stage4_interval_curve import root
from three_to_two_tone_stage4_targeted_transfer import e0,U

ROOT=Path(__file__).resolve().parent
iv.dps=45;mp.dps=75;getcontext().prec=110

def main():
    failures=json.loads((ROOT/'stage4_targeted_transfer_inner_tight.json').read_text())['failures']
    source=[r for r in csv.DictReader((ROOT/'stage2_inner_cells.csv').open()) if r['label']=='crossing']
    ref={(r['branch'],r['status'],r['depth'],r['u1_lo'],r['u2_lo']):r for r in source}
    db=iv.mpf(['1.9999','2']);si.D_STRONG=db
    eb=iv.mpf([str(e0-mp.mpf('2.5e-5')),str(e0+mp.mpf('1e-6'))])
    upper=(mp.sqrt(U)+mp.mpf('4e-4'))**2
    a,b,c,d,_=root(mp.mpf('1.99995'))
    centers={'A':(Decimal(mp.nstr(a,65)),Decimal(mp.nstr(b,65))),
             'B':(Decimal(mp.nstr(c,65)),Decimal(mp.nstr(d,65)))}
    # Both are strictly inside larger parameter-uniform positive-Hessian
    # boxes of radii .0003 (A) and .008 (B).
    radii={'A':Decimal('.00029'),'B':Decimal('.0079')}
    counts={};nvisited=0;remaining=[];per=[]
    for i,f in enumerate(failures):
        r=ref[(f['branch'],f['type'],f['depth'],f['u1_lo'],f['u2_lo'])]
        stack=[(Decimal(r['u1_lo']),Decimal(r['u1_hi']),
                Decimal(r['u2_lo']),Decimal(r['u2_hi']),0)]
        start=nvisited;localremaining=0
        while stack:
            xlo,xhi,ylo,yhi,dep=stack.pop();nvisited+=1
            if nvisited%1000==0:print('inside',i,'visited',nvisited,'stack',len(stack),'depth',dep,flush=True)
            lab=f['branch'];(cx,cy)=centers[lab];rad=radii[lab]
            if xlo>=cx-rad and xhi<=cx+rad and ylo>=cy-rad and yhi<=cy+rad:
                status='inside_convex_core'
            else:
                X=[iv.mpf([str(xlo),str(xhi)]),iv.mpf([str(ylo),str(yhi)])]
                mid=[iv.mpf(str((xlo+xhi)/2)),iv.mpf(str((ylo+yhi)/2))]
                J0=si.objective_jet(mid[0],mid[1],eb)
                JX=si.objective_jet(X[0],X[1],eb)
                Y=[X[k]-mid[k] for k in range(2)]
                L=J0.v+sum(J0.g[k]*Y[k] for k in range(2))
                L+=sum(iv.mpf('.5')*JX.H[k][j]*Y[k]*Y[j] for k in range(2) for j in range(2))
                if L.a>iv.mpf(str(upper)).b:
                    status='objective_excluded'
                else:
                    G=[J0.g[k]+sum(JX.H[k][j]*Y[j] for j in range(2)) for k in range(2)]
                    if any(g.a>0 or g.b<0 for g in G):status='gradient_excluded'
                    else:
                        status=''
                        try:
                            H=np.array([[float(J0.H[k][j].mid) for j in range(2)] for k in range(2)])
                            C0=np.linalg.inv(H)
                            C=[[iv.mpf(repr(float(C0[k,j]))) for j in range(2)] for k in range(2)]
                            M=[[iv.mpf(int(k==j))-sum(C[k][l]*JX.H[l][j] for l in range(2))
                                for j in range(2)] for k in range(2)]
                            K=[mid[k]-sum(C[k][j]*J0.g[j] for j in range(2))
                               +sum(M[k][j]*Y[j] for j in range(2)) for k in range(2)]
                            if any(K[k].b<X[k].a or K[k].a>X[k].b for k in range(2)):
                                status='krawczyk_excluded'
                        except (ValueError,ZeroDivisionError,np.linalg.LinAlgError):pass
            if status:
                counts[status]=counts.get(status,0)+1
                continue
            if dep>=28 or nvisited>=100000:
                localremaining+=1;remaining.append(dict(original=i,branch=lab,
                    u1_lo=str(xlo),u1_hi=str(xhi),u2_lo=str(ylo),u2_hi=str(yhi),depth=dep))
                continue
            if xhi-xlo>=yhi-ylo:
                mid=(xlo+xhi)/2
                stack.extend([(xlo,mid,ylo,yhi,dep+1),(mid,xhi,ylo,yhi,dep+1)])
            else:
                mid=(ylo+yhi)/2
                stack.extend([(xlo,xhi,ylo,mid,dep+1),(xlo,xhi,mid,yhi,dep+1)])
        per.append(dict(original=i,branch=f['branch'],visited=nvisited-start,remaining=localremaining))
        if i%10==0:print(i,'of',len(failures),'visited',nvisited,'remaining',len(remaining),flush=True)
    out=dict(original_failure_cells=len(failures),visited=nvisited,counts=counts,
             remaining=remaining,per_original_cell=per,
             verified=not remaining,
             note='covers only the 81 original inner failure cells; other archived inner cells use prior recheck')
    (ROOT/'stage4_targeted_transfer_inner_refined.json').write_text(json.dumps(out,indent=2))
    print(json.dumps({k:out[k] for k in ('original_failure_cells','visited','counts','verified')},indent=2))

if __name__=='__main__':main()
