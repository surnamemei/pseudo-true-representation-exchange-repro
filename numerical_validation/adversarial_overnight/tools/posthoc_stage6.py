"""POST_HOC_SENSITIVITY for Stage 6: nature of OUT labels and nearest-root probit widths.
Does not replace the pre-registered Stage-6 labels or hypotheses."""
import csv, json, math
import numpy as np
from scipy.special import ndtr, ndtri
from scipy.optimize import minimize
rows = list(csv.DictReader(open('stage6_noise/trials.csv')))
out = {"tag": "POST_HOC_SENSITIVITY"}
outs = [r for r in rows if r['label'] == 'OUT']
dA = np.array([float(r['dA']) for r in outs]); dB = np.array([float(r['dB']) for r in outs])
out['n_OUT'] = len(outs)
out['OUT_by_snr'] = {s: sum(r['snr'] == s for r in outs) for s in ('20', '30', '40')}
out['OUT_nearest_A_fraction'] = float(np.mean(dA < dB)) if len(outs) else None
out['OUT_dmin_quantiles'] = np.quantile(np.minimum(dA, dB), [0, .25, .5, .75, .9, 1]).tolist() if len(outs) else None
out['OUT_cls_counts'] = {c: sum(r['best_cls'] == c for r in outs) for c in ('regular', 'coalescent')}
etas = sorted(set(float(r['eta']) for r in rows))
widths = {}
for s in ('20', '30', '40'):
    k = []; n = []
    for e in etas:
        rr = [r for r in rows if r['snr'] == s and abs(float(r['eta']) - e) < 1e-12]
        near_A = sum(float(r['dA']) < float(r['dB']) for r in rr)
        k.append(near_A); n.append(len(rr))
    k = np.array(k, float); n = np.array(n, float); eta = np.array(etas)
    def nll(q):
        p = np.clip(ndtr(q[0] + q[1] * eta), 1e-12, 1 - 1e-12)
        return -np.sum(k * np.log(p) + (n - k) * np.log1p(-p))
    f = minimize(nll, [0, -10.0], method='Nelder-Mead', options=dict(xatol=1e-10, fatol=1e-12, maxiter=4000))
    widths[s] = dict(width_nearest_root_labels=float((ndtri(.9) - ndtri(.1)) / abs(f.x[1])),
                     P_nearestA=list(map(float, k / n)))
out['probit_width_nearest_root'] = widths
json.dump(out, open('stage6_noise/POSTHOC_labels_and_widths.json', 'w'), indent=1)
print(json.dumps({k: v for k, v in out.items() if k != 'probit_width_nearest_root'}, indent=1))
for s, w in widths.items(): print(s, round(w['width_nearest_root_labels'], 4), [round(x, 3) for x in w['P_nearestA']])
