"""N31 data-flow audit with poisoned N21 metadata and inner-root defaults."""
import csv,json,hashlib
from pathlib import Path
import numpy as np
from mpmath import mp,iv
ROOT=Path(__file__).resolve().parent
def main():
    import three_to_two_tone_stage3_full_crossings as fc
    import three_to_two_tone_stage2_inner_global as inner
    import three_to_two_tone_stage2_interval as it
    from stage7_n31_certificate import local_check
    cert=json.loads((ROOT/'N31_global_certificate.json').read_text());r=cert['high_precision_numerical_root'];ep=cert['epsilon_bracket']
    mp.dps=80;seed=tuple(mp.mpf(r[k]) for k in ('u_A1','u_A2','u_B1','u_B2','epsilon_cross'))
    a=fc.crossing(2,N=31,seed=seed)
    old=fc.ref.copy();fc.ref.update(lambda_N='123',lambda_quartic='456',v_A='999',v_B='998',A_coeff=['0','0','888'],B_coeff=['0','0','777'])
    b=fc.crossing(2,N=31,seed=seed);fc.ref=old
    changes=[k for k in a if a[k]!=b[k]]
    invariant=['epsilon_cross','u_A1','u_A2','u_B1','u_B2','J_A','J_B','Rscaled_A','Rscaled_B']
    assert all(a[k]==b[k] for k in invariant)
    it.N=31;iv.dps=90
    centers=[[r['u_A1'],r['u_A2']],[r['u_B1'],r['u_B2']]]
    U=min(local_check(centers[0],ep)[1].b,local_check(centers[1],ep)[1].b)
    inner.ref={'A_at_crossing':['999','998'],'B_at_crossing':['997','996']}
    iv.dps=45;logs=[];tests=[]
    for w,c in zip(('A','B'),centers):
        done,n,rows,remaining=inner.process('crossing',w,ep,np.array(list(map(float,c))),U,max_cells=50000,known_center_override=c,local_radius_override='0.000299999999' if w=='A' else '0.009999999')
        logs.extend(rows);tests.append(dict(branch=w,complete=done,visited=n,remaining=remaining))
        print('poisoned defaults',tests[-1],flush=True)
    original=list(csv.DictReader(open(ROOT/'N31_inner_interval_cells.csv')))
    # csv's string serialization is the archived comparison convention.
    equal=len(logs)==len(original) and all({k:str(v) for k,v in row.items()}==ref for row,ref in zip(logs,original))
    result=dict(poisoned_quartic_changed_fields=changes,invariant_root_fields=invariant,root_fields_exactly_equal=True,poisoned_inner_defaults=tests,inner_rows_identical=equal,inner_row_count=len(logs),frozen_sha256=hashlib.sha256((ROOT/'N31_global_certificate.json').read_bytes()).hexdigest())
    (ROOT/'stage11_N31_poison_test.json').write_text(json.dumps(result,indent=2));print(result,flush=True)
if __name__=='__main__':main()
