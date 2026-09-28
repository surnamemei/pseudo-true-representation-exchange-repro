"""Independent five-column least-squares check of interval scalar costs."""
import json
from pathlib import Path
from mpmath import mp
import stage11_tangent_global as g
from three_to_two_tone_stage3_tangent import tangent
mp.dps=75
R=Path(__file__).resolve().parent
r=json.loads((R/'smallz_global_certificate.json').read_text());tests=[]
for v in ['-.9','-.4','-.01','.01','.4','.9']+[str(x['v']) for x in r['stationary_points']]:
    q=tangent(mp.mpf(v),mp.mpf(g.ref['lambda_N']))[0]
    z=g.calc(float(v),float(v))[0]
    tests.append(dict(v=v,interval_contains_independent_5column_cost=mp.mpf(z.a)<=q<=mp.mpf(z.b),cost=str(q)))
(R/'stage11_tangent_crosschecks.json').write_text(json.dumps(tests,indent=2))
assert all(t['interval_contains_independent_5column_cost'] for t in tests)
print('Independent tangent cost checks:',len(tests),'passed')
