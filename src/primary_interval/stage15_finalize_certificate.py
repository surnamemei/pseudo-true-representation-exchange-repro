"""Assemble the effective-N0 certificate from a completed all-t interval replay."""
import hashlib
import json
import csv
from decimal import Decimal, getcontext
from pathlib import Path
from mpmath import mp, iv

ROOT=Path(__file__).resolve().parent
getcontext().prec=70
N0=35377
trial=json.loads((ROOT/f'stage15_interval_trial_N{N0}.json').read_text())
assert trial['verified'] and trial['first_failures']==[] and trial['root_boxes_linked']
old=json.loads((ROOT/'continuum_global_certificate.json').read_text())
replay=json.loads((ROOT/'stage14_continuum_global_replay.json').read_text())
local=json.loads((ROOT/'stage14_continuum_local_interval.json').read_text())
assert old['verified_global'] and replay['replay_verified'] and local['verified_local']
with (ROOT/'stage14_continuum_global_partition.csv').open(newline='',encoding='utf-8') as f:
    cells=list(csv.DictReader(f))
ordered=sorted(cells,key=lambda r:Decimal(r['lo']))
cover=(Decimal(ordered[0]['lo'])==Decimal(-100) and Decimal(ordered[-1]['hi'])==Decimal(100)
       and all(Decimal(a['hi'])==Decimal(b['lo']) for a,b in zip(ordered,ordered[1:])))
root_cells=[r for r in cells if r['predicate'].startswith('unique_root_')]
K=trial['local']['K']
linked_exact=(len(root_cells)==2 and all(
    Decimal(r['lo'])<Decimal(K[0 if '_A_' in r['predicate'] else 1][0]) and
    Decimal(K[0 if '_A_' in r['predicate'] else 1][1])<Decimal(r['hi'])
    for r in root_cells))
assert cover and linked_exact

mp.dps=80
iv.dps=80
eta=mp.sinh(mp.mpf('0.5'))/mp.mpf('0.5')-1
def rem(n):
    z=(mp.mpf(165)/n)**2
    return mp.exp(mp.mpf('27.5'))/(1-eta)*z**9/(1-z)
assert rem(N0)<mp.mpf('1e-30') and rem(N0-2)>mp.mpf('1e-30')
half=iv.mpf(1)/2
eta_iv=((iv.exp(half)-iv.exp(-half))/2)/half-1
def rem_iv(n):
    z=(iv.mpf(165)/n)**2
    return iv.exp(iv.mpf(55)/2)/(1-eta_iv)*z**9/(1-z)
Bn=rem_iv(N0);Bprev=rem_iv(N0-2)
assert Bn.b<iv.mpf('1e-30').a and Bprev.a>iv.mpf('1e-30').b

cm={
 'curvature':min(Decimal(local['Rvv_A'][0]),Decimal(local['Rvv_B'][0])),
 'satellite_amplitude':min(Decimal(local['beta_A'][0]),-Decimal(local['beta_B'][1])),
 'transversality':Decimal(local['crossing_slope'][0]),
 'regular_competitor_cost':Decimal(old['C6_other_regular_margin_lower']),
 'coalescent_cost':Decimal(old['coalescent_margin_lower'])
}
tm={
 'curvature':min(Decimal(trial['local']['curvature_A_lower']),Decimal(trial['local']['curvature_B_lower'])),
 'satellite_amplitude':min(Decimal(trial['local']['beta_A_lower']),Decimal(trial['local']['beta_B_abs_lower'])),
 'transversality':Decimal(trial['local']['slope_lower']),
 'regular_competitor_cost':Decimal(trial['minimum_predicate_margins']['cost_excluded']),
 'coalescent_cost':Decimal(trial['coalescent_margin_lower'])
}
erosion={k:max(Decimal(0),cm[k]-tm[k]) for k in cm}
assert all(tm[k]>cm[k]/2 for k in cm)

def sha(name):
    return hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

out={
 'status':'verified_interval_assisted_effective_threshold',
 'family':{'b':10,'phi':'pi','N_parity':'odd'},
 'N0_analytic':None,
 'N0_analytic_reason':'Closed-form global rational-loss sup-norm constants were not reduced to a usable stand-alone threshold.',
 'kernel_remainder_first_allowed_odd_N':N0,
 'kernel_remainder_at_N0':mp.nstr(rem(N0),35),
 'kernel_remainder_at_previous_odd_N':mp.nstr(rem(N0-2),35),
 'kernel_remainder_at_N0_directed_interval':[str(mp.mpf(Bn.a)),str(mp.mpf(Bn.b))],
 'kernel_remainder_previous_odd_directed_interval':[str(mp.mpf(Bprev.a)),str(mp.mpf(Bprev.b))],
 'N0_interval':N0,
 'minimality_scope':'Smallest odd N satisfying the selected analytic 1e-30 kernel remainder inequality; its all-t interval predicates pass at that floor. Not the smallest true N with branch exchange.',
 't_interval':['0',trial['t_range'][1]],
 'Krawczyk_inclusion':trial['local']['inclusion'],
 'Krawczyk_contraction_upper':trial['local']['contraction_upper'],
 'partition_cells':sum(trial['predicate_counts'].values()),
 'partition_predicate_counts':trial['predicate_counts'],
 'partition_all_passed':trial['failed_count_capped']==0,
 'partition_cover_from_frozen_replay':replay['partition_cover_exact_decimal'],
 'partition_cover_rechecked_exact_decimal':cover,
 'root_cells_linked_rechecked_exact_decimal':linked_exact,
 'continuum_margin_lower':{k:str(v) for k,v in cm.items()},
 'uniform_finite_margin_lower':{k:str(v) for k,v in tm.items()},
 'E_i_uniform_erosion_bound_for_all_N_at_least_N0':{k:str(v) for k,v in erosion.items()},
 'safety_factor_two_passes':all(tm[k]>cm[k]/2 for k in cm),
 'finite_tail_margin_lower':trial['tail_margin_lower'],
 'finite_coalescent_margin_lower':trial['coalescent_margin_lower'],
 'frozen_partition_sha256':sha('stage14_continuum_global_partition.csv'),
 'checker_sha256':sha('stage15_uniform_interval.py'),
 'trial_file':f'stage15_interval_trial_N{N0}.json',
 'coverage':'Every odd N >= N0; only a possibly N-dependent sufficiently small z neighborhood.',
 'not_covered':'No connected proof to z=2 and no assertion for every odd 3 <= N < N0.'
}
(ROOT/'explicit_N0_certificate.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps({k:out[k] for k in ('status','N0_analytic','N0_interval','partition_cells','safety_factor_two_passes')},indent=2))
