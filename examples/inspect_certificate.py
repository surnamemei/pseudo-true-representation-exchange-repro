"""Print the key fields of the main certificate records (read-only)."""
import csv
import json
from pathlib import Path

C = Path(__file__).resolve().parent.parent / "certificates"


def show(name, keys):
    d = json.loads((C / name).read_text())
    print(f"== {name}")
    for k in keys:
        print(f"   {k}: {str(d.get(k))[:110]}")


show("stage12_replay/N21_summary.json", ["backend", "precision_bits", "predicates", "passed", "failures"])
show("stage12_replay/N31_summary.json", ["backend", "precision_bits", "predicates", "passed", "failures"])
show("continuum_global_certificate.json", ["verified_global", "coalescent_margin_lower", "tail_margin_lower"])
show("independent_N0_replay.json", ["backend", "N0", "partition_cells", "passed", "failed", "partition_cover"])
show("tightened_N0_certificate.json", ["new_sufficient_odd_N_floor", "all_uniform_conditions_pass", "scope", "minimality_claim"])
show("finite_width_epsilon_certificate.json", ["final_closed_epsilon_interval", "full_global_interval_certified",
                                               "unique_transverse_crossing", "scope"])
show("lower_transition_certificate.json", ["coalescent_excluded_for_all_positive_lambda", "scope"])
print("== generalN_certified_instances.csv")
for r in csv.DictReader((C / "generalN_certified_instances.csv").open()):
    print(f"   N={r['N']}: {r['status']}; lambda_N={r['lambda_N'][:14]}; c_N={r['c_N'][:15]}")
