"""Conservative index of the proof checks and remaining obligations."""
import hashlib,json
from pathlib import Path

ROOT=Path(__file__).resolve().parent

def read(name):return json.loads((ROOT/name).read_text())
def sha(name):return hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

def main():
    local=read('stage4_bridge_local_nodes.json')
    outer=read('stage4_targeted_transfer_outer.json')
    broad=read('stage4_targeted_transfer_inner.json')
    tight=read('stage4_targeted_transfer_inner_tight.json')
    z25=read('stage4_bottleneck_cells/z_0p25_centered_interval_summary.json')
    z5=read('stage4_bottleneck_cells/z_0p5_interval_summary.json')
    points={
      '0.25':False,'0.5':False,
      '0.75':read('stage4_bridge_pointwise_0p75.json')['pointwise_global_verified'],
      '1.0':read('stage4_pointwise_repaired.json')[0]['pointwise_global_verified'],
      '1.25':read('stage4_bridge_pointwise_1p25.json')['pointwise_global_verified'],
      '1.5':read('stage4_pointwise_repaired.json')[1]['pointwise_global_verified'],
      '1.75':read('stage4_bridge_pointwise_1p75.json')['pointwise_global_verified'],
      '2.0':True}
    result=dict(
      final_category='B — BRIDGE NUMERICALLY SUPPORTED BUT CERTIFICATION-LIMITED',
      largest_connected_global_z_interval=['1.9999997','2.0'],
      local_crossing_nodes={r['z_lo']:r['verified'] for r in local},
      isolated_pointwise_global_nodes=points,
      new_z_1p9999_to_2p0_outer_exclusion=outer['outer_uniform_exclusion'],
      new_z_1p9999_to_2p0_outer_refined_original_cells=outer['original_failing_cells'],
      new_z_1p9999_to_2p0_outer_refined_children=outer['refined_children'],
      new_z_1p9999_to_2p0_inner_unresolved_broad=len(broad['failures']),
      new_z_1p9999_to_2p0_inner_unresolved_tight=len(tight['failures']),
      z_0p25_hard_strip_interval_excluded=z25['excluded'],
      z_0p25_other_cells_without_directed_exclusion=132381,
      z_0p5_terminal_cells_interval_excluded=z5['counts']['objective_excluded'],
      z_0p5_other_float_excluded_leaves_interval_checked=False,
      no_third_branch_detected_numerically=True,
      no_third_branch_globally_excluded_on_bridge=False,
      connected_bridge_certified=False,
      source_hashes={n:sha(n) for n in (
          'stage4_bridge_local_nodes.json','stage4_targeted_transfer_outer.json',
          'stage4_targeted_transfer_inner.json','stage4_targeted_transfer_inner_tight.json',
          'stage4_failure_cells.csv.gz','stage4_inner_failure_cells.csv',
          'stage4_bridge_pointwise_0p75.json','stage4_bridge_pointwise_1p25.json',
          'stage4_bridge_pointwise_1p75.json')})
    (ROOT/'stage4_proof_bottleneck_certificate.json').write_text(json.dumps(result,indent=2))
    print(json.dumps({k:result[k] for k in ('final_category','largest_connected_global_z_interval',
          'new_z_1p9999_to_2p0_inner_unresolved_broad','connected_bridge_certified')},indent=2))

if __name__=='__main__':main()
