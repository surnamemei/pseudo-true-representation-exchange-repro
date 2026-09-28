"""Hash active Stage-17 proof and submission artifacts."""
import csv
import hashlib
import json
from pathlib import Path

root=Path(__file__).resolve().parent
epsilon=json.loads((root/'finite_width_epsilon_certificate.json').read_text())
fixed={
    'main PDF':['paper/main_tsp_reviewfriendly.pdf'],
    'supplement PDF':['supplement/supplement.pdf'],
    'active manuscript source':[
        'paper/main_tsp_reviewfriendly.tex','paper/sections/body_reviewfriendly.tex',
        'paper/sections/continuum_stage14.tex','paper/sections/evidence_table.tex',
        'paper/sections/main_certificate_table_10page.tex',
        'paper/figures/fig_stage17_tangent_landscape.pdf',
        'supplement/supplement.tex','supplement/S1b_continuum_largeN.tex',
        'supplement/S2_interval_certificate_protocol.tex',
        'supplement/S4_extended_generality.tex','supplement/S5_independent_replay.tex'],
    'Stage 17 report':[
        'REPRODUCIBILITY_README.md','paper/CURRENT_VERSION.md',
        'inner_cell_accounting.md','inner_cell_accounting.csv',
        'tightened_N0_derivation.md','tightened_N0_certificate.json',
        'largeN_condition_thresholds.csv',
        'finite_width_epsilon_certificate.md','finite_width_epsilon_certificate.json',
        'epsilon_interval_obligations.csv','structural_stability_corollary.md',
        'final_hostile_review.md','final_submission_blockers.md',
        'stage17_epsilon_coeff_crosscheck.json',
        'stage17_epsilon_coeff_crosscheck.py',
        'stage17_annotate_condition_thresholds.py',
        'stage17_interval_trial_N1501.json','stage17_interval_trial_N5001.json'],
    'Stage 17 interval proof log':list(epsilon['sha256']),
}
seen=set();rows=[]
for role,paths in fixed.items():
    for rel in paths:
        if rel in seen:continue
        seen.add(rel)
        p=root/rel
        assert p.is_file(),p
        h=hashlib.sha256()
        with p.open('rb') as f:
            for block in iter(lambda:f.read(4*1024*1024),b''):
                h.update(block)
        rows.append(dict(role=role,path=rel,bytes=p.stat().st_size,sha256=h.hexdigest()))
out=root/'stage17_artifact_manifest.csv'
with out.open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=['role','path','bytes','sha256'])
    w.writeheader();w.writerows(rows)
print(len(rows),'manifest entries')
