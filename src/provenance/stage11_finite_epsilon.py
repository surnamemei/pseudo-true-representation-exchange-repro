"""Finite-width parameter proof attempt and historical margin diagnosis."""
import csv,json,re
from pathlib import Path
from mpmath import iv,mp
import numpy as np
import three_to_two_tone_stage2_interval as it
import three_to_two_tone_stage2_outer_global as outer
from three_to_two_tone_stage3_full_crossings import value_grad
from stage7_n31_certificate import local_check
import three_to_two_tone_branch_study as st

ROOT=Path(__file__).resolve().parent
it.N=21;iv.dps=60;mp.dps=65
ref=json.loads((ROOT/'stage2_reference.json').read_text())
def num(s):return mp.mpf(re.findall(r'[-+]?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?',s)[0])
def ib(x):return [str(mp.mpf(x.a)),str(mp.mpf(x.b))]
rows=[r for r in csv.DictReader(open(ROOT/'interval_boxes.csv')) if r['label']=='crossing' and r['status']=='excluded']
rows.sort(key=lambda r:num(r['interval_J_lower']))
ec=mp.mpf(json.loads((ROOT/'stage2_crossing_global.json').read_text())['epsilon_bracket'][0])+mp.mpf('1e-15')

def main():
    smallest=rows[0]
    cell=(int(smallest['depth']),*[float(smallest[k]) for k in ('m_lo','m_hi','h_lo','h_hi')],0,0,'excluded')
    m=(cell[1]+cell[2])/2;h=(cell[3]+cell[4])/2
    u=np.array([21*(m-h),21*(m+h)])
    hh=mp.findroot(lambda a,b:value_grad(a,b,ec,2)[1:],tuple(u),tol=mp.mpf('1e-55'),maxsteps=80)
    stationary=value_grad(*hh,ec,2)
    coords=dict(original_cell=smallest,center_u=list(u),u1_interval=[21*(cell[1]-cell[4]),21*(cell[2]-cell[3])],u2_interval=[21*(cell[1]+cell[3]),21*(cell[2]+cell[4])],center_objective=str(value_grad(*u,ec,2)[0]),stationary_from_cell=[str(v) for v in hh],stationary_cost=str(stationary[0]),stationary_precision=mp.dps)
    for w in ('A','B'):
        q=np.array([float(v) for v in ref[w+'_at_crossing']])
        coords['center_Linf_distance_to_'+w]=float(np.max(abs(u-q)))
        coords['center_distance_outside_radius2_'+w]=float(max(0,np.max(abs(u-q))-2))
    # Independent stationary exclusion over the whole minimum-margin cell.
    stack=[(*coords['u1_interval'],*coords['u2_interval'],0)]
    visits=0;unresolved=[]
    gradient_epsilon=iv.mpf(json.loads((ROOT/'stage2_crossing_global.json').read_text())['epsilon_bracket'])
    while stack:
        a,b,c,d,dep=stack.pop();visits+=1
        X=[iv.mpf([a,b]),iv.mpf([c,d])]
        mid=[iv.mpf((a+b)/2),iv.mpf((c+d)/2)]
        JX=it.objective_jet(*X,gradient_epsilon)
        J0=it.objective_jet(*mid,gradient_epsilon)
        grad=[J0.g[i]+sum(JX.H[i][k]*(X[k]-mid[k]) for k in range(2)) for i in range(2)]
        if any(v.a>0 or v.b<0 for v in grad):continue
        if dep>=12:unresolved.append([a,b,c,d]);continue
        if b-a>=d-c:
            m=(a+b)/2;stack.extend([(a,m,c,d,dep+1),(m,b,c,d,dep+1)])
        else:
            m=(c+d)/2;stack.extend([(a,b,c,m,dep+1),(a,b,m,d,dep+1)])
    coords['gradient_exclusion_subdivision']=dict(visited=visits,unresolved=unresolved,no_stationary_point=not unresolved,epsilon_interval=ib(gradient_epsilon))
    out=[]
    # Twelve noninfinitesimal slabs. Test all local predicates and the 40
    # weakest historical spatial cells. Failures are retained, never promoted.
    for idx in range(12):
        lo=ec-mp.mpf('.03')+idx*mp.mpf('.005');hi=lo+mp.mpf('.005');mid=(lo+hi)/2
        ep=[str(lo),str(hi)];E=it.data_energy(iv.mpf(ep))
        rr=[];Js=[]
        for w in ('A','B'):
            q=mp.findroot(lambda a,b:value_grad(a,b,mid,2)[1:],tuple(map(mp.mpf,ref[w+'_at_crossing'])),tol=mp.mpf('1e-55'))
            point=it.objective_jet(iv.mpf(str(q[0])),iv.mpf(str(q[1])),iv.mpf(ep))
            Js.append(point.v.b) # a fixed feasible fit, uniform incumbent
            choices=[]
            for radius in ('.005','.01','.025','.05','.1'):
                cert,_,_=local_check([str(v) for v in q],ep,radius)
                choices.append(dict(radius=radius,inclusion=cert['inclusion'],positive_Hessian=cert['positive_Hessian'],contraction=cert['contraction_upper']))
            rr.append(dict(branch=w,center=[str(v) for v in q],attempts=choices))
        incumbent=min(Js);fail=[];passed=0
        for k,row in enumerate(rows[:40]):
            c=(int(row['depth']),*[float(row[n]) for n in ('m_lo','m_hi','h_lo','h_hi')],0,0,'excluded')
            ok,lb,why=outer.verify_cell(c,iv.mpf(ep),E,incumbent)
            if ok:passed+=1
            else:fail.append(dict(rank=k,cell={n:row[n] for n in ('m_lo','m_hi','h_lo','h_hi')},reason=why,lower=None if lb is None else ib(lb)))
        out.append(dict(epsilon_interval=ep,local=rr,uniform_incumbent=str(incumbent),tested_outer_cells=40,passed=passed,failures=fail))
        print('epsilon slab',idx,'local',[[q['inclusion'] and q['positive_Hessian'] for q in c['attempts']] for c in rr],'outer passed',passed,flush=True)
    result=dict(requested_interval=[str(ec-mp.mpf('.03')),str(ec+mp.mpf('.03'))],historical_isolated_cases_certified=['epsilon_c-0.035','epsilon_c','epsilon_c+0.035'],historical_outer_summaries=json.loads((ROOT/'stage2_outer_global_summary.json').read_text()),minimum_outer_cell=coords,attempts=out,finite_width_global_certified=False,reason='parameter dependency in fixed-frequency root boxes and reused outer bounds; failed cells explicitly retained; no complete finite-width partition',status='CERTIFICATION-LIMITED')
    (ROOT/'finite_epsilon_certificate.json').write_text(json.dumps(result,indent=2))
if __name__=='__main__':main()
