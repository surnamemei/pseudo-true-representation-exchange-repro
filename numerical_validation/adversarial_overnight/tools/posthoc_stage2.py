"""POST_HOC_SENSITIVITY: fixed-delta scaling using every global switch into the B-like fit,
regardless of the pre-registered C2 regularity of the lower-epsilon participant.
Does not replace the pre-registered (INCONCLUSIVE) Stage-2 classification."""
import json, glob, math, sys
import numpy as np
sys.path.insert(0, '.')
from s2_phase import fits, cfgB, DELTAS, ZS_B, ZFIT
out = {"tag": "POST_HOC_SENSITIVITY", "note": "switch into the B-like (satellite within 2 of 10.3) fit; P may be near-coalescent"}
for d in DELTAS:
    pts, rows = [], []
    for z in ZS_B:
        a = json.load(open(f"stage2_phase/analysis/{cfgB(d, z)['cid']}.json"))
        ex = [s for s in a['switches'] if s.get('type') == 'exchange']
        def sat(p): return p['u1'] if abs(p['u1']) > abs(p['u2']) else p['u2']
        cand = [s for s in ex if abs(sat(s['Q']) - 10.3) < 2.0]
        s = cand[0] if cand else None
        r = dict(z=z, n_exchanges=len(ex))
        if s:
            r.update(eps_c=s['eps_c'], eps_over_z=s['eps_c'] / z, status=s['status'], P_sep=s['P']['sep'],
                     P_cls=s['P']['cls'], P_amp=s['P']['amp_max'], Q_sat=sat(s['Q']), C5=s['C5'],
                     third_margin_rel=s['third_margin_rel'])
            if z in ZFIT and s['C5']:
                pts.append((z, s['eps_c']))
        rows.append(r)
    f = {}
    if len(pts) == len(ZFIT):
        zz, ee = zip(*pts)
        f = dict(full=fits(zz, ee), drop_smallest=fits(zz[1:], ee[1:]), drop_largest=fits(zz[:-1], ee[:-1]))
    out[str(d)] = dict(rows=rows, fits=f, tangent_eps_over_z=0.6274474 * math.sin(d * math.pi))
json.dump(out, open('stage2_phase/POSTHOC_fixed_delta_all_switches.json', 'w'), indent=1)
for d in DELTAS:
    o = out[str(d)]
    print(d, 'tangent', round(o['tangent_eps_over_z'], 5), [(r['z'], round(r.get('eps_over_z', float('nan')), 5), r.get('P_cls'), round(r.get('P_sep', float('nan')), 3), r.get('C5')) for r in o['rows']])
    if o['fits']:
        F = o['fits']
        print('   alpha', round(F['full']['alpha'], 4), 'rel rms z1', round(F['full']['relres_rms_z1'], 4), 'z2', round(F['full']['relres_rms_z2'], 4), '| drop small', round(F['drop_smallest']['alpha'], 4), 'drop large', round(F['drop_largest']['alpha'], 4))
