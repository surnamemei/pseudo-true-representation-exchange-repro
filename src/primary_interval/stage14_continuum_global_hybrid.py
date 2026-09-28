"""Full-line continuum global exclusion using cost, derivative, and root predicates."""
import csv
import json
from pathlib import Path
from mpmath import iv,mp
import stage14_continuum_global_interval as core

ROOT=Path(__file__).resolve().parent
iv.dps=core.iv.dps;mp.dps=core.mp.dps
ref=json.loads((ROOT/'stage14_continuum_candidate.json').read_text(encoding='utf-8'))
local=json.loads((ROOT/'stage14_continuum_local_interval.json').read_text(encoding='utf-8'))
va=float(ref['v_A']);vb=float(ref['v_B'])
U=core.calc(va,va)[0].b

def main():
    stack=[(k/4,(k+1)/4) for k in range(-400,400)]
    leaves=[];roots=[];un=[];vis=0
    while stack:
        a,b=stack.pop();vis+=1
        z=core.calc(a,b)
        if z is not None and z[0].a>U:
            leaves.append((a,b,'cost_excluded'));continue
        if z is not None and core.excludes(z[1]):
            leaves.append((a,b,'gradient_excluded'));continue
        if z is not None and core.excludes(z[2]):
            ga=core.calc(a,a)[1];gb=core.calc(b,b)[1]
            if (ga.a>0 and gb.a>0) or (ga.b<0 and gb.b<0):
                leaves.append((a,b,'monotone_derivative_no_root'));continue
            if (ga.b<0 and gb.a>0) or (ga.a>0 and gb.b<0):
                lo,hi=a,b
                for _ in range(35):
                    mid=(lo+hi)/2;gm=core.calc(mid,mid)[1]
                    if not core.excludes(gm):break
                    if (gm.a>0)==(ga.a>0):lo=mid
                    else:hi=mid
                rv=core.calc(lo,hi)
                kind='min' if z[2].a>0 else 'max'
                identity='A' if abs((lo+hi)/2-va)<.25 else ('B' if abs((lo+hi)/2-vb)<.25 else 'other')
                if kind=='min' and identity=='other' and rv[0].a<=U:
                    # This should only happen from interval overwidth.
                    m=(a+b)/2;stack.extend([(a,m),(m,b)]);continue
                roots.append({'lo':lo,'hi':hi,'v':(lo+hi)/2,'kind':kind,'identity':identity,
                              'R_lo':str(mp.mpf(rv[0].a)),'R_hi':str(mp.mpf(rv[0].b)),
                              'H_lo':str(mp.mpf(rv[2].a)),'H_hi':str(mp.mpf(rv[2].b))})
                leaves.append((a,b,'unique_root_'+identity+'_'+kind));continue
        if b-a<1e-10:
            un.append([a,b]);continue
        m=(a+b)/2;stack.extend([(a,m),(m,b)])
        if vis%1000==0:print('visited',vis,'pending',len(stack),flush=True)
    roots.sort(key=lambda z:z['v'])
    mins=[r for r in roots if r['kind']=='min']
    other=[r for r in mins if r['identity']=='other']
    other_margin=min((mp.mpf(r['R_lo'])-mp.mpf(U) for r in other),default=mp.inf)
    q2=-core.M4+core.L*core.D(iv.mpf(10),2)-core.M2*core.q0
    coal=core.C-q2*q2/(core.M4-core.M2*core.M2)
    V=iv.mpf(100)
    dd=2/V;dpb=(1+dd)/V;d2b=1/(2*V)+2*dpb/V;dbw=2/(V-10)
    AA=1-dd*dd-dpb*dpb/core.M2
    BB=d2b+core.L.b*dbw+dd*abs(core.q0).b+dpb*abs(core.q1).b/(2*core.M2.a)
    tail=core.C-BB*BB/AA
    out={
      'domain':[-100,100],'lambda_box':core.bounds(core.L),'precision_dps':iv.dps,
      'visited':vis,'unresolved':un,'retained_stationary_count':len(roots),
      'retained_minimum_count':len(mins),'all_stationary_roots_isolated':False,
      'other_roots_excluded_by_cost_or_derivative_predicates':True,
      'root_identities':[r['identity'] for r in mins],
      'global_trial_upper':str(mp.mpf(U)),
      'other_isolated_root_margin_lower':str(other_margin),
      'coalescent_interval':core.bounds(coal),
      'coalescent_margin_lower':str(mp.mpf(coal.a)-mp.mpf(U)),
      'tail_lower_interval':core.bounds(tail),
      'tail_margin_lower':str(mp.mpf(tail.a)-mp.mpf(U)),
      'verified_global':bool(not un and len([r for r in mins if r['identity']=='A'])==1
                             and len([r for r in mins if r['identity']=='B'])==1
                             and other_margin>0 and coal.a>U and tail.a>U and local['verified_local'])
    }
    (ROOT/'continuum_global_certificate.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
    with (ROOT/'stage14_continuum_global_roots.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(roots[0]));w.writeheader();w.writerows(roots)
    with (ROOT/'stage14_continuum_global_partition.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.writer(f);w.writerow(['lo','hi','predicate']);w.writerows(leaves)
    print(json.dumps(out,indent=2),flush=True)

if __name__=='__main__':main()
