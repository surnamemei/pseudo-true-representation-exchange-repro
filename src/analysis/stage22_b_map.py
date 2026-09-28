"""High-precision continuum b continuation with numerical competitor diagnostics."""
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from mpmath import mp
from scipy.optimize import minimize_scalar

ROOT = Path(__file__).resolve().parent
mp.dps = 65
M2 = mp.mpf(1)/12
M4 = mp.mpf(1)/80
V2 = M4-M2*M2


def d(x, k=0):
    x = mp.mpf(x)
    if abs(x) < 1:
        return mp.fsum((-1)**m*x**(2*m-k)/
                       (mp.mpf(2*m+1)*mp.factorial(2*m-k)*2**(2*m))
                       for m in range((k+1)//2, 52))
    out = 2*mp.sin(x/2)/x
    for j in range(1, k+1):
        out = (2*mp.sin(x/2+j*mp.pi/2)/2**j-j*out)/x
    return out


def parts(v, lam, b):
    dv, dp = d(v), d(v, 1)
    q0 = -M2-lam*d(b)
    q1 = 2*lam*d(b, 1)
    gram = 1-dv*dv-dp*dp/M2
    numer = d(v, 2)-lam*d(b-v)-dv*q0+dp*q1/(2*M2)
    base = M4-2*lam*d(b, 2)+lam*lam-q0*q0-q1*q1/(4*M2)
    return gram, numer, base


def loss(v, lam, b):
    a, bb, cc = parts(v, lam, b)
    return cc-bb*bb/a


def coal(lam, b):
    q0 = -M2-lam*d(b)
    q2 = -M4+lam*d(b,2)-M2*q0
    return parts(mp.mpf(2),lam,b)[2]-q2*q2/V2


def solve(b, seed):
    b = mp.mpf(b)
    f1 = lambda va, vb, lam: mp.diff(lambda v:loss(v,lam,b), va)
    f2 = lambda va, vb, lam: mp.diff(lambda v:loss(v,lam,b), vb)
    f3 = lambda va, vb, lam: loss(va,lam,b)-loss(vb,lam,b)
    va, vb, lam = mp.findroot((f1,f2,f3), seed, tol=mp.mpf('1e-52'),maxsteps=70)
    ca, cb = loss(va,lam,b),loss(vb,lam,b)
    ga,ba,_ = parts(va,lam,b)
    gb,bb,_ = parts(vb,lam,b)
    curv_a=mp.diff(lambda v:loss(v,lam,b),va,2)
    curv_b=mp.diff(lambda v:loss(v,lam,b),vb,2)
    slope=mp.diff(lambda L:loss(va,L,b)-loss(vb,L,b),lam)
    coal_gap=coal(lam,b)-ca
    out = dict(b=float(b),lambda_inf=float(lam),v_A=float(va),v_B=float(vb),
               R_star=float((ca+cb)/2),beta_A=float(ba/ga),beta_B=float(bb/gb),
               curvature_A=float(curv_a),curvature_B=float(curv_b),
               slope=float(slope),coalescent_margin=float(coal_gap),
               lambda_inf_hp=mp.nstr(lam,55),v_A_hp=mp.nstr(va,55),
               v_B_hp=mp.nstr(vb,55),R_star_hp=mp.nstr((ca+cb)/2,55),
               beta_A_hp=mp.nstr(ba/ga,55),beta_B_hp=mp.nstr(bb/gb,55),
               curvature_A_hp=mp.nstr(curv_a,55),curvature_B_hp=mp.nstr(curv_b,55),
               slope_hp=mp.nstr(slope,55),
               coalescent_margin_hp=mp.nstr(coal_gap,55))
    out['_mp']=(va,vb,lam)
    return out


def d_np(x, k):
    x=np.asarray(x,float)
    safe=np.where(np.abs(x)<.12,1.,x)
    if k==0:
        y=2*np.sin(safe/2)/safe
        t=1-x*x/24+x**4/1920-x**6/322560
    elif k==1:
        y=np.cos(safe/2)/safe-2*np.sin(safe/2)/safe**2
        t=-x/12+x**3/480-x**5/53760
    else:
        y=-np.sin(safe/2)/(2*safe)-2*np.cos(safe/2)/safe**2+4*np.sin(safe/2)/safe**3
        t=-1/12+x*x/160-x**4/10752
    return np.where(np.abs(x)<.12,t,y)


def loss_np(v, lam, b):
    dv,dp=d_np(v,0),d_np(v,1)
    q0=-1/12-lam*d_np(b,0)
    q1=2*lam*d_np(b,1)
    a=1-dv*dv-12*dp*dp
    bb=d_np(v,2)-lam*d_np(b-v,0)-dv*q0+6*dp*q1
    cc=1/80-2*lam*d_np(b,2)+lam*lam-q0*q0-3*q1*q1
    with np.errstate(divide='ignore',invalid='ignore'):
        return np.where(np.abs(v)<.2,np.nan,cc-bb*bb/a)


def tail_bound(lam,b):
    v=100.
    dmax=2/v
    dpmax=1/v+2/v**2
    ddmax=.5/v+2/v**2+4/v**3
    q0=-1/12-lam*float(d(b))
    q1=2*lam*float(d(b,1))
    bmax=ddmax+abs(lam)*2/(v-b)+dmax*abs(q0)+dpmax*abs(q1)*6
    amin=1-dmax*dmax-12*dpmax*dpmax
    cc=float(parts(mp.mpf(2),mp.mpf(lam),mp.mpf(b))[2])
    return cc-bmax*bmax/amin


def regular_scan(out):
    b,lam=out['b'],out['lambda_inf']
    grid=np.linspace(-100,100,8001)
    vals=loss_np(grid,lam,b)
    ix=np.where(np.isfinite(vals[1:-1]) & (vals[1:-1]<vals[:-2]) &
                (vals[1:-1]<vals[2:]))[0]+1
    minima=[]
    for k in ix:
        lo,hi=grid[k-1],grid[k+1]
        res=minimize_scalar(lambda x:float(loss(mp.mpf(x),mp.mpf(lam),mp.mpf(b))),
                            bounds=(lo,hi),method='bounded',
                            options={'xatol':1e-12})
        if res.success:
            minima.append((res.x,res.fun))
    # Explicit narrow endpoint probes so a near-zero regular well is not silently missed.
    for lo,hi in [(-.2,-.025),(.025,.2)]:
        res=minimize_scalar(lambda x:float(loss(mp.mpf(x),mp.mpf(lam),mp.mpf(b))),
                            bounds=(lo,hi),method='bounded',options={'xatol':1e-12})
        if res.success and abs(res.x-lo)>.00001 and abs(res.x-hi)>.00001:
            minima.append((res.x,res.fun))
    minima.sort(key=lambda r:r[0])
    dedup=[]
    for x,y in minima:
        if not dedup or abs(x-dedup[-1][0])>.03:
            dedup.append((x,y))
    third=[(x,y) for x,y in dedup if abs(x-out['v_A'])>.05 and
           abs(x-out['v_B'])>.05]
    best_third=min(third,key=lambda r:r[1]) if third else (float('nan'),float('nan'))
    out['regular_minima_scanned']=len(dedup)
    out['nearest_regular_competitor_v']=best_third[0]
    out['nearest_regular_competitor_margin']=best_third[1]-out['R_star']
    out['tail_lower_margin']=tail_bound(lam,b)-out['R_star']
    return out


def main():
    anchor=(mp.mpf('-4.060030216180106892789'),mp.mpf('12.42390445039620145'),
            mp.mpf('0.0665069703923955354559'))
    result={}
    for direction in ((10,9,8,7,6),(11,12,13,14)):
        seed=anchor if direction[0]==10 else result[10]['_mp']
        for b in direction:
            try:
                row=solve(b,seed)
                regular_scan(row)
                seed=row['_mp']
                row['status']=('A_clean_numerical' if row['lambda_inf']>0 and
                               row['curvature_A']>0 and row['curvature_B']>0 and
                               abs(row['beta_A'])>1e-10 and abs(row['beta_B'])>1e-10 and
                               abs(row['slope'])>1e-10 and
                               row['coalescent_margin']>0 and
                               row['nearest_regular_competitor_margin']>0 and
                               row['tail_lower_margin']>0 else 'requires_review')
                if row['status']=='A_clean_numerical' and (abs(row['v_A'])<.1 or
                    row['coalescent_margin']<1e-6 or abs(row['beta_A'])>1e3):
                    row['status']='A_near_coalescent_numerical'
                result[b]=row
                print(b,row['lambda_inf'],row['v_A'],row['v_B'],
                      row['coalescent_margin'],row['nearest_regular_competitor_margin'],
                      row['status'],flush=True)
            except Exception as exc:
                result[b]=dict(b=b,status='E_no_candidate',error=str(exc))
                print(b,'FAILED',exc,flush=True)
    rows=[]
    for b in range(6,15):
        row=dict(result[b])
        row.pop('_mp',None)
        rows.append(row)
    fields=list(dict.fromkeys(k for row in rows for k in row))
    with (ROOT/'lambda_inf_vs_b.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=fields)
        w.writeheader();w.writerows(rows)
    ok=[row for row in rows if 'lambda_inf' in row and row['b']<=13]
    fig,ax=plt.subplots(figsize=(5.8,3.4))
    ax.plot([r['b'] for r in ok],[r['lambda_inf'] for r in ok],
            '-o',color='#116ba2',lw=1.8,markersize=4,
            label='Numerical continuum A/B crossing')
    ax.scatter([10],[result[10]['lambda_inf']],marker='s',s=52,color='#622b88',
               zorder=4,label='Independently certified $b=10$')
    if 14 in result and 'lambda_inf' in result[14]:
        ax.scatter([14],[result[14]['lambda_inf']],marker='^',s=58,
                   edgecolor='#b35b12',facecolor='none',linewidth=1.5,zorder=4,
                   label='Near central chart (numerical)')
    ax.set(xlabel=r'Weak-tone location $b$',ylabel=r'Critical coefficient $\lambda_\infty(b)$')
    ax.set_xticks(range(6,15))
    ax.grid(alpha=.2)
    ax.legend(frameon=False,fontsize=8)
    fig.tight_layout()
    fig.savefig(ROOT/'lambda_inf_vs_b.pdf')
    plt.close(fig)


if __name__=='__main__':main()
