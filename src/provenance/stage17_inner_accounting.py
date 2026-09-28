"""Exact archive/replay accounting for the Stage 17 inner-cell discrepancy."""
from collections import Counter,defaultdict
import csv,json
from pathlib import Path

root=Path(__file__).resolve().parent
files={21:root/"stage2_inner_cells.csv",31:root/"N31_inner_interval_cells.csv"}
keyfields=("u1_lo","u1_hi","u2_lo","u2_hi")
rows=[]
all_good=True
for n,file in files.items():
    original=list(csv.DictReader(file.open(newline="")))
    replay=list(csv.DictReader((root/f"stage12_replay/N{n}_predicates.csv").open(newline="")))
    leaves=list(csv.DictReader((root/f"stage12_replay/N{n}_refinement_leaves.csv").open(newline="")))
    summary=json.loads((root/f"stage12_replay/N{n}_summary.json").read_text())
    inner=[r for r in replay if r["kind"]=="inner"]
    assert len(inner)==len(original)
    for i,r in enumerate(original):
        p=inner[i]
        assert int(p["index"])==i and p["case"]==r.get("label","crossing")
        assert p["branch"]==r["branch"] and p["original_predicate"]==r["status"]
        assert p["agreement"]=="PASS"
    for label,branch in sorted({(r.get("label","crossing"),r["branch"]) for r in original}):
        group=[r for r in original if r.get("label","crossing")==label and r["branch"]==branch]
        gp=[r for r in inner if r["case"]==label and r["branch"]==branch]
        gl=[r for r in leaves if r["case"]==label and r["branch"]==branch]
        tree=next((r for r in json.loads((root/"stage2_inner_global_summary.json").read_text())
                   if n==21 and r["label"]==label and r["branch"]==branch),None)
        if n==31:
            tree=next(r for r in json.loads((root/"N31_global_certificate.json").read_text())["inner"]
                      if r["branch"]==branch)
        assert tree is not None and tree["complete"] and tree["remaining"]==0
        coords=[tuple(r[k] for k in keyfields) for r in group]
        duplicates=len(coords)-len(set(coords))
        parent_keys={r["parent"] for r in gl}
        orig_keys={",".join(r[k] for k in keyfields) for r in group}
        assert parent_keys<=orig_keys
        child_counts=Counter(r["parent"] for r in gl)
        assert all(v==2 for v in child_counts.values())
        refined=len(parent_keys)
        terminals=len(group)
        visited=int(tree["visited"]);internal=visited-terminals
        assert terminals==int(tree["terminal"]) and visited==2*terminals-1
        assert len(gp)==terminals and duplicates==0
        coverage=json.loads((root/"stage12_replay/coverage_audit.json").read_text())
        region=next(r for r in coverage if r["N"]==n and r["case"]==label
                    and r["kind"]=="inner_"+branch)
        assert region["cells"]==terminals and region["complete"]
        assert region["uncovered_slabs"]==[] and region["expected_neighborhood_matches"]
        # Every replay refinement has two accepted child leaves and a coverage audit.
        cov=json.loads((root/"stage12_replay/refinement_coverage_audit.json").read_text())
        cove={(r["N"],r["case"],r["branch"],r["parent"]):r for r in cov}
        assert all(cove[(n,label,branch,p)]["complete"]
                   and cove[(n,label,branch,p)]["parent_domain_matches"]
                   and cove[(n,label,branch,p)]["cells"]==2 for p in parent_keys)
        statuses=Counter(r["status"] for r in group)
        alternative=sum(p["original_predicate"]!=p["arb_predicate"] for p in gp)
        root_accepted=sum("inside" in r["status"].lower() for r in group)
        rows.append(dict(N=n,case=label,branch=branch,
                         archived_visited=visited,archived_internal_parents=internal,
                         archived_terminal_obligations=terminals,
                         archived_duplicate_terminal_coordinates=duplicates,
                         accepted_root_box_terminals=root_accepted,
                         arb_replayed_archived_terminals=len(gp),
                         arb_refined_archived_parents=refined,
                         arb_refined_child_leaves=len(gl),
                         arb_effective_terminal_leaves=terminals-refined+len(gl),
                         arb_alternate_predicate_proofs=alternative,
                         objective_excluded=statuses["objective_excluded"],
                         gradient_excluded=statuses["gradient_excluded"],
                         krawczyk_excluded=statuses["krawczyk_excluded"],
                         all_pass=True))
    assert sum(x["arb_refined_archived_parents"] for x in rows if x["N"]==n)==summary["refined_archived_inner_cells"]
    assert 3*summary["refined_archived_inner_cells"]==summary["refinement_nodes"]
with (root/"inner_cell_accounting.csv").open("w",newline="",encoding="utf-8") as f:
    writer=csv.DictWriter(f,fieldnames=list(rows[0]))
    writer.writeheader();writer.writerows(rows)
print(json.dumps({
    "groups":len(rows),
    "all_terminal":sum(r["archived_terminal_obligations"] for r in rows),
    "crossing_visited":sum(r["archived_visited"] for r in rows if r["case"]=="crossing"),
    "crossing_terminal":sum(r["archived_terminal_obligations"] for r in rows if r["case"]=="crossing"),
    "auxiliary_terminal":sum(r["archived_terminal_obligations"] for r in rows if r["case"]!="crossing"),
    "refined_parents":sum(r["arb_refined_archived_parents"] for r in rows),
    "refined_children":sum(r["arb_refined_child_leaves"] for r in rows),
    "arb_effective_leaves":sum(r["arb_effective_terminal_leaves"] for r in rows),
    "root_box_terminal":sum(r["accepted_root_box_terminals"] for r in rows),
    "duplicates":sum(r["archived_duplicate_terminal_coordinates"] for r in rows)
},indent=2))
