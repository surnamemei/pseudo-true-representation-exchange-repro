"""Noiseless, finite-record 3-to-2 tone projection study.

Run from F:/Research with: python three_to_two_tone_branch_study.py
This is numerical evidence, not a computer-assisted proof of globality.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import minimum_filter
from scipy.optimize import brentq, minimize, minimize_scalar

import three_to_two_tone_explore as ex

ROOT = Path(__file__).resolve().parent
FIG = ROOT / 'three_to_two_tone_figures'
FIG.mkdir(exist_ok=True)


def set_n(n):
    assert n % 2 == 1
    ex.N = n
    ex.t = np.arange(n) - (n - 1) / 2


def x_record(n, d, b, eps, phase):
    set_n(n)
    return ex.record(d, b, eps, phase)


def wrap_u(u, n):
    return np.sort((np.asarray(u) + np.pi * n) % (2 * np.pi * n) - np.pi * n)


def value(u, x, n):
    return ex.cost(np.asarray(u) / n, x)


def refine(u0, x, n):
    r = minimize(value, np.asarray(u0, float), args=(x, n), method='Nelder-Mead',
                 options={'xatol': 2e-10, 'fatol': 2e-12, 'maxiter': 1400})
    u = wrap_u(r.x, n)
    return u, value(u, x, n)


def hessian(u, x, n, h=0.002):
    e1 = np.array([h, 0.0]); e2 = np.array([0.0, h])
    f0 = value(u, x, n)
    h11 = (value(u+e1,x,n)-2*f0+value(u-e1,x,n))/h**2
    h22 = (value(u+e2,x,n)-2*f0+value(u-e2,x,n))/h**2
    h12 = (value(u+e1+e2,x,n)-value(u+e1-e2,x,n)
           -value(u-e1+e2,x,n)+value(u-e1-e2,x,n))/(4*h**2)
    return np.linalg.eigvalsh([[h11,h12],[h12,h22]])


def stats(u, x, n):
    set_n(n)
    v = np.column_stack((ex.atom(u[0]/n),ex.atom(u[1]/n)))
    amps = np.linalg.lstsq(v,x,rcond=None)[0]
    gram_cond = np.linalg.cond(v.conj().T@v)
    eig = hessian(u,x,n)
    waveform = v@amps
    return {'sep_u': float(u[1]-u[0]),
            'max_abs_amp': float(np.max(np.abs(amps))),
            'gram_cond': float(gram_cond),
            'hess_min': float(eig[0]),
            'hess_max': float(eig[1]),
            'waveform': waveform}


def two_branches(n,d,b,eps,phase,seed_a=None,seed_b=None):
    x = x_record(n,d,b,eps,phase)
    ua,ja = refine([-d,d] if seed_a is None else seed_a,x,n)
    ub,jb = refine([0,b] if seed_b is None else seed_b,x,n)
    return (ua,ja),(ub,jb),x


def independent_minima(x,n,grid_size=512,random_count=0,seed=13):
    set_n(n)
    grid=np.linspace(-np.pi,np.pi,grid_size,endpoint=False)
    arr=ex.grid_cost(x,grid)
    filt=minimum_filter(arr,size=3,mode='wrap')
    ii=np.argwhere(np.isfinite(arr)&(arr<=filt+1e-11))
    ii=sorted(ii,key=lambda k:arr[tuple(k)])
    starts=[np.sort(grid[k]*n) for k in ii[:200]]
    rng=np.random.default_rng(seed)
    starts.extend(np.sort(rng.uniform(-np.pi*n,np.pi*n,size=(random_count,2)),axis=1))
    found=[]
    for st in starts:
        u,j=refine(st,x,n)
        if not any(np.linalg.norm(u-r['u'])<1e-3 for r in found):
            found.append({'u':u,'J':j,'grid_seed':len(starts)})
    found.sort(key=lambda r:r['J'])
    return arr,grid,found,len(ii)


def boundary_cost(x,n):
    set_n(n)
    def one_cost(u):
        v=ex.atom(u/n)
        return float(np.vdot(x,x).real-abs(np.vdot(v,x))**2/n)
    def confluent_cost(u):
        v=np.column_stack((ex.atom(u/n),1j*ex.t*ex.atom(u/n)))
        c=np.linalg.lstsq(v,x,rcond=None)[0]
        e=x-v@c
        return float(np.vdot(e,e).real)
    grid=np.linspace(-np.pi*n,np.pi*n,1001)
    out=[]
    for f in (one_cost,confluent_cost):
        vals=np.array([f(v) for v in grid])
        ix=np.where((vals<=np.roll(vals,1))&(vals<=np.roll(vals,-1)))[0]
        opts=[minimize_scalar(f,bounds=(grid[max(0,i-1)],grid[min(len(grid)-1,i+1)]),
                              method='bounded',options={'xatol':1e-11}) for i in ix]
        out.append(min([z.fun for z in opts]+[float(np.min(vals))]))
    return out


def crossing(n,d,b,phase,lo=.12,hi=.5):
    def diff(eps):
        (ua,ja),(ub,jb),_=two_branches(n,d,b,eps,phase)
        if np.linalg.norm(ua-ub)<.1:
            return np.nan
        return ja-jb
    ee=np.linspace(lo,hi,20)
    vv=np.array([diff(e) for e in ee])
    for i in range(len(ee)-1):
        if np.isfinite(vv[i]) and np.isfinite(vv[i+1]) and vv[i]*vv[i+1]<0:
            root=brentq(diff,ee[i],ee[i+1],xtol=1e-11)
            (ua,ja),(ub,jb),x=two_branches(n,d,b,root,phase)
            sa=stats(ua,x,n); sb=stats(ub,x,n)
            return root,ua,ub,ja,jb,sa,sb
    return None


def save_csv(path,rows):
    if not rows:
        return
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0].keys()))
        w.writeheader();w.writerows(rows)


def figure(name):
    plt.tight_layout();plt.savefig(FIG/name,dpi=180);plt.close()


def main():
    n,d,b,phase=21,2.0,10.0,np.pi
    cross=crossing(n,d,b,phase)
    if cross is None:
        raise RuntimeError('No crossing in primary geometry')
    ec,ua_c,ub_c,ja_c,jb_c,sa_c,sb_c=cross
    print('PRIMARY CROSS',ec,ua_c,ub_c,ja_c,jb_c,'seps',sa_c['sep_u'],sb_c['sep_u'],flush=True)
    eps_values=np.linspace(max(.12,ec-.11),ec+.11,91)
    data=[]
    for eps in eps_values:
        (ua,ja),(ub,jb),x=two_branches(n,d,b,float(eps),phase,ua_c,ub_c)
        sa=stats(ua,x,n);sb=stats(ub,x,n)
        wave_dist=float(np.linalg.norm(sa['waveform']-sb['waveform']))
        for label,u,j,s in [('A',ua,ja,sa),('B',ub,jb,sb)]:
            data.append({'N':n,'d_NDelta':d,'b_Nomega3':b,'phase_rad':phase,
                         'epsilon':float(eps),'branch':label,'J':j,'u1_Nnu1':u[0],
                         'u2_Nnu2':u[1],'sep_Nrad':s['sep_u'],
                         'max_abs_amp':s['max_abs_amp'],'gram_cond':s['gram_cond'],
                         'hess_min_u':s['hess_min'],'hess_max_u':s['hess_max'],
                         'waveform_branch_distance':wave_dist,
                         'selected_of_AB':int((label=='A' and ja<jb) or (label=='B' and jb<ja))})
    save_csv(ROOT/'three_to_two_tone_branch_data.csv',data)
    minima_rows=[];checks=[]
    for eps in (ec-.035,ec,ec+.035):
        x=x_record(n,d,b,eps,phase)
        arr,grid,mins,grid_count=independent_minima(x,n,grid_size=1024,random_count=120)
        co,conf=boundary_cost(x,n)
        checks.append({'eps':eps,'best':mins[0]['J'],'second':mins[1]['J'],
                       'one':co,'confluent':conf,'minima_count':len(mins),
                       'grid_seed_count':grid_count})
        for rank,z in enumerate(mins):
            u=z['u'];s=stats(u,x,n)
            minima_rows.append({'N':n,'d_NDelta':d,'b_Nomega3':b,'phase_rad':phase,
                                'epsilon':eps,'rank':rank+1,'J':z['J'],
                                'u1_Nnu1':u[0],'u2_Nnu2':u[1],
                                'sep_Nrad':s['sep_u'],'max_abs_amp':s['max_abs_amp'],
                                'gram_cond':s['gram_cond'],'hess_min_u':s['hess_min'],
                                'one_tone_boundary_J':co,'confluent_boundary_J':conf,
                                'grid_size':1024,'random_starts':120})
        # Plot the informative portion while searching the full alias-free circle.
        gu=grid*n;mask=(gu>=-19)&(gu<=20)
        shown=np.clip(arr[np.ix_(mask,mask)],0,4)
        plt.figure(figsize=(7,6))
        plt.imshow(shown,origin='lower',extent=(gu[mask][0],gu[mask][-1],gu[mask][0],gu[mask][-1]),
                   aspect='equal',cmap='viridis',vmin=0,vmax=4)
        plt.colorbar(label='Profiled SSE (clipped at 4)')
        (ua_plot,_),(ub_plot,_),_=two_branches(n,d,b,eps,phase,ua_c,ub_c)
        for lab,u in [('A',ua_plot),('B',ub_plot)]:
            plt.plot(u[1],u[0],'o',ms=8,label=lab)
        plt.xlabel('u2 = N nu2');plt.ylabel('u1 = N nu1');plt.title(f'Objective landscape, epsilon={eps:.6f}')
        plt.legend();figure(f'01_landscape_epsilon_{eps:.6f}.png')
        print('GLOBAL CHECK',checks[-1],flush=True)
    save_csv(ROOT/'three_to_two_tone_minima.csv',minima_rows)
    eps=np.array(eps_values)
    by={q:[r for r in data if r['branch']==q] for q in ('A','B')}
    def ys(q,key): return np.array([r[key] for r in by[q]])
    plt.figure(figsize=(8,5))
    for q in ('A','B'):
        for k,sty in [('u1_Nnu1','-'),('u2_Nnu2','--')]:plt.plot(eps,ys(q,k),sty,label=f'{q}: {k.split("_")[0]}')
    plt.axvline(ec,color='k',lw=1);plt.xlabel('epsilon');plt.ylabel('Fitted normalized frequency u=N nu');plt.legend(ncol=2)
    figure('02_fitted_frequencies_vs_epsilon.png')
    plt.figure(figsize=(7,5))
    for q in ('A','B'):plt.plot(eps,ys(q,'J'),label=f'Branch {q}')
    plt.axvline(ec,color='k',lw=1);plt.xlabel('epsilon');plt.ylabel('Profiled SSE');plt.legend()
    figure('03_branch_objectives_vs_epsilon.png')
    plt.figure(figsize=(7,5));plt.plot(eps,ys('A','J')-ys('B','J'))
    plt.axhline(0,color='k',lw=.8);plt.axvline(ec,color='k',lw=.8)
    plt.xlabel('epsilon');plt.ylabel('J_A - J_B');figure('04_deltaJ_vs_epsilon.png')
    plt.figure(figsize=(7,5))
    for q in ('A','B'):plt.plot(eps,ys(q,'hess_min_u'),label=f'Branch {q}')
    plt.axvline(ec,color='k',lw=.8);plt.xlabel('epsilon');plt.ylabel('Smallest profiled Hessian eigenvalue (u units)');plt.legend()
    figure('05_hessian_minimum_eigenvalue.png')
    plt.figure(figsize=(7,5))
    for q in ('A','B'):plt.plot(eps,ys(q,'sep_Nrad'),label=f'Branch {q}')
    plt.axvline(ec,color='k',lw=.8);plt.xlabel('epsilon');plt.ylabel('N |nu2 - nu1|');plt.legend()
    figure('06_fitted_frequency_separation.png')
    boundary=[]
    for nn in (17,21,25,31):
        for dd in (1.6,1.8,2.0,2.2,2.4):
            z=crossing(nn,dd,b,phase,lo=.06,hi=.75)
            if z:
                e,ua,ub,ja,jb,sa,sb=z
                boundary.append({'N':nn,'d_NDelta':dd,'b_Nomega3':b,'phase_rad':phase,
                                 'epsilon_cross':e,'u_distance':float(np.linalg.norm(ua-ub)),
                                 'min_hess':min(sa['hess_min'],sb['hess_min']),
                                 'max_gram_cond':max(sa['gram_cond'],sb['gram_cond'])})
                print('BOUNDARY',boundary[-1],flush=True)
    plt.figure(figsize=(7,5))
    for nn in sorted(set(r['N'] for r in boundary)):
        rr=[r for r in boundary if r['N']==nn]
        plt.plot([r['d_NDelta'] for r in rr],[r['epsilon_cross'] for r in rr],'-o',label=f'N={nn}')
    plt.xlabel('N Delta');plt.ylabel('epsilon at local-branch crossing');plt.legend()
    figure('07_crossing_boundary_vs_NDelta.png')
    robust=[]
    for ph in np.pi+np.array([-.3,-.15,0,.15,.3]):
        for bb in (9,9.5,10,10.5,11):
            z=crossing(n,d,bb,float(ph),lo=.10,hi=.55)
            if z:
                e,ua,ub,ja,jb,sa,sb=z
                robust.append({'N':n,'d_NDelta':d,'b_Nomega3':bb,'phase_rad':ph,
                               'epsilon_cross':e,'u_distance':float(np.linalg.norm(ua-ub)),
                               'min_hess':min(sa['hess_min'],sb['hess_min']),
                               'max_gram_cond':max(sa['gram_cond'],sb['gram_cond'])})
    mat=np.full((5,5),np.nan)
    for r in robust:
        i=np.argmin(abs((np.pi+np.array([-.3,-.15,0,.15,.3]))-r['phase_rad']))
        j=np.argmin(abs(np.array([9,9.5,10,10.5,11])-r['b_Nomega3']))
        mat[i,j]=r['epsilon_cross']
    plt.figure(figsize=(7,5));plt.imshow(mat,origin='lower',aspect='auto',extent=[8.75,11.25,-.375,.375],cmap='viridis')
    plt.colorbar(label='Crossing epsilon');plt.xlabel('Weak-tone position b=N omega3');plt.ylabel('Phase minus pi (rad)')
    figure('08_phase_location_robustness.png')
    # Separate existence checks for a small perturbation of the strong amplitude ratio.
    amp_tests=[]
    for ratio in (.95,1.0,1.05):
        def rec(e):
            set_n(n)
            return ex.atom(-d/n)+ratio*ex.atom(d/n)+e*np.exp(1j*phase)*ex.atom(b/n)
        def diff(e):
            x=rec(e);ua,ja=refine([-d,d],x,n);ub,jb=refine([0,b],x,n)
            return ja-jb if np.linalg.norm(ua-ub)>.1 else np.nan
        ee=np.linspace(.15,.4,20);vv=np.array([diff(e) for e in ee])
        roots=[]
        for i in range(len(ee)-1):
            if np.isfinite(vv[i]) and np.isfinite(vv[i+1]) and vv[i]*vv[i+1]<0:
                roots.append(brentq(diff,ee[i],ee[i+1],xtol=1e-10))
        amp_tests.append({'strong_amp_ratio':ratio,'epsilon_cross':roots[0] if roots else np.nan})
    save_csv(ROOT/'three_to_two_tone_robustness.csv',boundary+robust)
    print('ROBUST COUNT',len(robust),'AMP TESTS',amp_tests,flush=True)
    print('CHECKS',checks,flush=True)
    print('OUTPUT',FIG,flush=True)


if __name__=='__main__':
    main()
