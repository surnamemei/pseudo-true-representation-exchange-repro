"""Directed-interval Taylor/gradient/interval-Newton checks for selected cells."""
import csv,gzip,json
from pathlib import Path
import numpy as np
from mpmath import iv,mp
import three_to_two_tone_stage2_interval as si
from three_to_two_tone_stage4_interval_curve import root,validate

ROOT=Path(__file__).resolve().parent;OUT=ROOT/'stage4_bottleneck_cells';N=21
iv.dps=35;mp.dps=75

def ivstr(x):return str(x.a),str(x.b)

def check_cell(r,z,e,U,do_newton=True):
    a,b,c,d=[mp.mpf(r[k]) for k in ('m_lo','m_hi','h_lo','h_hi')]
    u1=(N*(a-d),N*(b-c));u2=(N*(a+c),N*(b+d))
    X=[iv.mpf([str(u1[0]),str(u1[1])]),iv.mpf([str(u2[0]),str(u2[1])])]
    mid=[iv.mpf(str((u1[0]+u1[1])/2)),iv.mpf(str((u2[0]+u2[1])/2))]
    Y=[X[i]-mid[i] for i in range(2)]
    try:
        J0=si.objective_jet(mid[0],mid[1],e)
        JX=si.objective_jet(X[0],X[1],e)
    except (ZeroDivisionError,ValueError):
        return dict(cell_id=int(r['cell_id']),status='undefined_interval_gram')
    L=J0.v+sum(J0.g[i]*Y[i] for i in range(2))
    L+=sum(iv.mpf('0.5')*JX.H[i][j]*Y[i]*Y[j] for i in range(2) for j in range(2))
    if L.a>U:
        return dict(cell_id=int(r['cell_id']),status='objective_excluded',
                    Jlower=str(L.a),margin=str(L.a-U))
    G=[J0.g[i]+sum(JX.H[i][k]*Y[k] for k in range(2)) for i in range(2)]
    if any(g.a>0 or g.b<0 for g in G):
        return dict(cell_id=int(r['cell_id']),status='gradient_excluded',
                    Jlower=str(L.a),margin=str(L.a-U),
                    grad1=ivstr(G[0]),grad2=ivstr(G[1]))
    if do_newton:
        try:
            C0=np.linalg.inv(np.array([[float(J0.H[i][j].mid) for j in range(2)] for i in range(2)]))
            C=[[iv.mpf(repr(float(C0[i,j]))) for j in range(2)] for i in range(2)]
            M=[[iv.mpf(int(i==j))-sum(C[i][k]*JX.H[k][j] for k in range(2))
                for j in range(2)] for i in range(2)]
            K=[mid[i]-sum(C[i][k]*J0.g[k] for k in range(2))
               +sum(M[i][k]*Y[k] for k in range(2)) for i in range(2)]
            if any(K[i].b<X[i].a or K[i].a>X[i].b for i in range(2)):
                return dict(cell_id=int(r['cell_id']),status='newton_excluded',
                            Jlower=str(L.a),margin=str(L.a-U))
        except (np.linalg.LinAlgError,ZeroDivisionError,ValueError):pass
    return dict(cell_id=int(r['cell_id']),status='unresolved',
                Jlower=str(L.a),margin=str(L.a-U),
                grad1=ivstr(G[0]),grad2=ivstr(G[1]),
                u1_bounds=list(map(str,u1)),u2_bounds=list(map(str,u2)))

def setup(z):
    si.D_STRONG=iv.mpf(str(z))
    a,b,c,d,ec=root(mp.mpf(str(z)))
    # At small z the 55-digit interval evaluation has a rounding floor above
    # 1e-48.  The 1e-20 Krawczyk box is still far narrower than any cell.
    q=validate(str(z),str(z),base_radius='1e-20')
    e=iv.mpf([mp.nstr(ec-mp.mpf('1e-19'),65),mp.nstr(ec+mp.mpf('1e-19'),65)])
    A=si.objective_jet(iv.mpf(mp.nstr(a,65)),iv.mpf(mp.nstr(b,65)),e)
    B=si.objective_jet(iv.mpf(mp.nstr(c,65)),iv.mpf(mp.nstr(d,65)),e)
    U=min(A.v.b,B.v.b)
    return e,U,q['verified']

def sample(z,n=200):
    p=OUT/f'z_{str(z).replace(".","p")}_unresolved.csv.gz'
    with gzip.open(p,'rt',newline='') as f:rows=list(csv.DictReader(f))
    t=OUT/f'z_{str(z).replace(".","p")}_taylor.csv.gz'
    with gzip.open(t,'rt',newline='') as f:diag=list(csv.DictReader(f))
    margins=np.array([float(r['Taylor_gap_float']) for r in diag])
    rng=np.random.default_rng(412)
    ids=list(rng.choice(len(rows),size=min(n//2,len(rows)),replace=False))
    ids+=list(np.argsort(margins)[:min(n//2,len(rows))])
    ids=list(dict.fromkeys(ids))
    e,U,root_ok=setup(z)
    checks=[check_cell(rows[i],z,e,U) for i in ids]
    counts={k:sum(q['status']==k for q in checks) for k in set(q['status'] for q in checks)}
    out=dict(z=z,samples=len(ids),counts=counts,root_verified=root_ok,
             cells=checks,status='directed_interval_sample')
    (OUT/f'z_{str(z).replace(".","p")}_interval_sample.json').write_text(json.dumps(out,indent=2))
    return {k:out[k] for k in ('z','samples','counts','root_verified')}

def batch(z):
    p=OUT/f'z_{str(z).replace(".","p")}_unresolved.csv.gz'
    with gzip.open(p,'rt',newline='') as f:rows=list(csv.DictReader(f))
    e,U,root_ok=setup(z)
    dest=OUT/f'z_{str(z).replace(".","p")}_interval_cells.csv.gz'
    counts={};minmargin=None
    with gzip.open(dest,'wt',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['cell_id','status','Jlower','margin','grad1','grad2','u1_bounds','u2_bounds'],extrasaction='ignore')
        w.writeheader()
        for i,r in enumerate(rows):
            q=check_cell(r,z,e,U)
            w.writerow(q)
            counts[q['status']]=counts.get(q['status'],0)+1
            if 'margin' in q:
                val=float(q['margin'].strip('[]').split(',')[0])
                minmargin=val if minmargin is None else min(minmargin,val)
            if i%500==0:print(z,i,'of',len(rows),counts,flush=True)
    out=dict(z=z,total=len(rows),counts=counts,root_verified=root_ok,
             min_recorded_margin=minmargin,
             status='directed_interval_checks_on_float_unresolved_partition')
    (OUT/f'z_{str(z).replace(".","p")}_interval_summary.json').write_text(json.dumps(out,indent=2))
    return out

if __name__=='__main__':
    for z in (.25,.5):print(sample(z),flush=True)
