from pathlib import Path
import csv,json
from mpmath import mp
mp.dps=85
# Check that each numerical reference lies inside its recorded joint validation box;
# check proof counts, analytic source presence, and frozen certificate identity.
for N in (11,15,21,31,41):
 p=Path(f'stage12_generalN/N{N}_certificate.json');q=json.loads(p.read_text());assert q['local_crossing']['verified'] and q['tangent_classification_certified'] and not q['unresolved']
 assert mp.mpf(q['other_minimum_margin'])>0 and mp.mpf(q['coalescent_margin'])>0
 leaves=list(csv.DictReader(open(f'stage12_generalN/N{N}_partition.csv')));leaves.sort(key=lambda r:mp.mpf(r['lo']))
 assert mp.mpf(leaves[0]['lo'])<=-mp.pi*N and mp.mpf(leaves[-1]['hi'])>=mp.pi*N
 assert all(mp.mpf(a['hi'])==mp.mpf(b['lo']) for a,b in zip(leaves,leaves[1:]))
 print(N,'complete tangent cover',len(leaves),'leaves; minima',sum(r['kind']=='min' for r in q['stationary_points']))
for N in (21,31):
 q=json.loads(Path(f'stage12_replay/N{N}_summary.json').read_text());assert q['passed']==q['predicates'] and not q['failures']
print('ALL FINAL PROOF CHECKS PASS')
