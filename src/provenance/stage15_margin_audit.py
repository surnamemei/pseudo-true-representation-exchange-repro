"""Read-only replay of Stage-14 partition margins relevant to an effective N0."""
import csv
import json
from pathlib import Path
from mpmath import mp
import stage14_continuum_global_interval as core

ROOT = Path(__file__).resolve().parent
with (ROOT / 'stage14_continuum_global_partition.csv').open(encoding='utf-8') as f:
    cells = list(csv.DictReader(f))
U = mp.mpf(json.loads((ROOT / 'continuum_global_certificate.json').read_text())['global_trial_upper'])

out = {}
for i, cell in enumerate(cells):
    lo, hi = float(cell['lo']), float(cell['hi'])
    p = cell['predicate']
    z = core.calc(lo, hi)
    if z is None:
        continue
    if p == 'cost_excluded':
        margin = mp.mpf(z[0].a) - U
    elif p == 'gradient_excluded':
        margin = max(mp.mpf(z[1].a), -mp.mpf(z[1].b))
    elif p == 'monotone_derivative_no_root':
        ga, gb = core.calc(lo, lo)[1], core.calc(hi, hi)[1]
        margin = min(max(mp.mpf(ga.a), -mp.mpf(ga.b)),
                     max(mp.mpf(gb.a), -mp.mpf(gb.b)),
                     max(mp.mpf(z[2].a), -mp.mpf(z[2].b)))
    else:
        margin = max(mp.mpf(z[2].a), -mp.mpf(z[2].b))
    if p not in out or margin < mp.mpf(out[p]['minimum_margin']):
        out[p] = {'minimum_margin': str(margin), 'cell': [cell['lo'], cell['hi']],
                  'width': str(hi-lo), 'index': i}
    if (i + 1) % 5000 == 0:
        print('audited', i + 1, flush=True)
(ROOT / 'stage15_partition_margin_audit.json').write_text(json.dumps(out, indent=2), encoding='utf-8')
print(json.dumps(out, indent=2))
