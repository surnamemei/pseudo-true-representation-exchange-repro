"""Parametric root, Gram, and slope checks for the finite epsilon interval."""
import json
from pathlib import Path
from mpmath import iv, mp
import stage7_n31_certificate as local
import stage17_parametric_jet as param
import three_to_two_tone_stage2_interval as it

ROOT=Path(__file__).resolve().parent
it.N=21
local.interval.N=21
local.interval.objective_jet=param.objective
iv.dps=90;mp.dps=100
ref=json.loads((ROOT/'stage2_reference.json').read_text())
cross=json.loads((ROOT/'stage2_crossing_global.json').read_text())
ec=(mp.mpf(cross['epsilon_bracket'][0])+mp.mpf(cross['epsilon_bracket'][1]))/2
r=mp.mpf('0.00001')
e=[str(ec-r),str(ec+r)]
checks={}
def gram_bounds(X,eps):
    def D(v):return iv.mpf(1)+sum(2*iv.cos(iv.mpf(k)*v/21) for k in range(1,11))
    u1,u2=X;e=iv.mpf(eps)
    g=D(u2-u1);den=iv.mpf(441)-g*g
    h1=D(u1+2)+D(u1-2)-e*D(u1-10)
    h2=D(u2+2)+D(u2-2)-e*D(u2-10)
    a1=(21*h1-g*h2)/den;a2=(21*h2-g*h1)/den
    return dict(Gram_determinant=[str(den.a),str(den.b)],Gram_positive=bool(den.a>0),
                amplitudes=[[str(v.a),str(v.b)] for v in (a1,a2)])
for w,rad in [('A','0.0003'),('B','0.01')]:
    cert,_,X=local.local_check(ref[w+'_at_crossing'],e,rad)
    gram=gram_bounds(X,e)
    checks[w]=dict(krawczyk=cert,gram=gram)
# local_check output X entries are interval_str strings; use the returned
# frequency boxes from a fresh call to avoid parsing formatted intervals.
_,_,XA=local.local_check(ref['A_at_crossing'],e,'0.0003')
_,_,XB=local.local_check(ref['B_at_crossing'],e,'0.01')
ca=param.coefficients(*XA)
cb=param.coefficients(*XB)
da=2*ca[0].v*iv.mpf(e)+ca[1].v
db=2*cb[0].v*iv.mpf(e)+cb[1].v
slope=da-db
out=dict(epsilon_interval=e,crossing_bracket=cross['epsilon_bracket'],
         branches=checks,slope_interval=[str(slope.a),str(slope.b)],
         slope_positive=bool(slope.a>0),
         all_checks_pass=bool(all(v['krawczyk']['inclusion'] and
                                  v['krawczyk']['contraction_lt_one'] and
                                  v['krawczyk']['positive_Hessian'] and
                                  v['gram']['Gram_positive'] for v in checks.values())
                              and slope.a>0))
(ROOT/'stage17_epsilon_local.json').write_text(json.dumps(out,indent=2))
print(json.dumps(dict(epsilon_interval=e,
                      A={k:checks['A']['krawczyk'][k] for k in ('inclusion','positive_Hessian','contraction_upper')},
                      B={k:checks['B']['krawczyk'][k] for k in ('inclusion','positive_Hessian','contraction_upper')},
                      gram_A=checks['A']['gram']['Gram_positive'],
                      gram_B=checks['B']['gram']['Gram_positive'],
                      slope=out['slope_interval'],all_checks_pass=out['all_checks_pass']),indent=2))
