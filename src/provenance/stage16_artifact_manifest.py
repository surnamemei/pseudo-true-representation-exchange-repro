"""Write a deterministic SHA-256 manifest of the active submission package."""
from pathlib import Path
import csv
import hashlib

root=Path(__file__).resolve().parent
groups={
 "anonymous manuscript":[
  "paper/main_tsp_reviewfriendly.tex","paper/main_tsp_reviewfriendly.pdf",
  "paper/sections/body_reviewfriendly.tex","paper/sections/continuum_stage14.tex",
  "paper/sections/evidence_table.tex","paper/sections/main_certificate_table_10page.tex",
  "paper/refs.bib","paper/figures/fig1_continuation.pdf","paper/figures/fig2_objective.pdf",
  "paper/figures/fig3_critical_law.pdf","paper/figures/fig7_phase_resolution.pdf",
  "continuum_vs_finiteN_figure.pdf"],
 "supplement":[
  "supplement/supplement.tex","supplement/supplement.pdf",
  "supplement/S1_complete_local_proof.tex","supplement/S1b_continuum_largeN.tex",
  "supplement/S2_interval_certificate_protocol.tex","supplement/S3_certificate_tables.tex",
  "supplement/S4_extended_generality.tex","supplement/S5_independent_replay.tex"],
 "finite length certificates":[
  "stage12_generalN/N11_certificate.json","stage12_generalN/N15_certificate.json",
  "stage12_generalN/N21_certificate.json","stage12_generalN/N31_certificate.json",
  "stage12_generalN/N41_certificate.json",
  "validated_globality_certificate.json","N31_global_certificate.json",
  "certificate_replay_report.md","certificate_replay_comparison.csv"],
 "continuum and threshold certificates":[
  "stage14_continuum_candidate.json","stage14_continuum_local_interval.json",
  "stage14_continuum_global_partition.csv","stage14_continuum_global_replay.json",
  "continuum_global_certificate.json","stage14_continuum_endpoint_interval.json",
  "stage15_uniform_interval.py","stage15_finalize_certificate.py",
  "stage15_interval_trial_N35377.json","explicit_N0_certificate.json",
  "stage16_independent_arb.py","independent_N0_replay.json",
  "independent_N0_replay.csv","independent_N0_replay.md"],
 "asymptotic derivation and validation":[
  "largeN_expansion_derivation.md","lambdaN_asymptotic_theorem.md",
  "stage14_asymptotics.py","stage14_asymptotic_coefficients.json",
  "stage15_asymptotic_data.py","stage16_symbolic_asymptotics.py",
  "stage16_symbolic_asymptotics.json","lambdaN_asymptotics.csv",
  "lambdaN_largeN_data.csv","lambdaN_residual_scaling.csv",
  "stage16_residual_scaling.csv","lambda_asymptotic_independent_check.md",
  "largeN_theory_vs_data.pdf"],
 "numerical context":[
  "phase_crossing_map.csv","resolution_sweep.csv","noise_branch_trials.csv",
  "branch_probability.csv","noise_model_fit.csv","transition_width.csv"],
 "submission documentation":[
  "REPRODUCIBILITY_README.md","CODE_README.md","certificate_archive_index.md",
  "final_hostile_review.md","final_submission_blockers.md",
  "cover_letter_factual_summary.md","author_metadata_checklist.md",
  "paper/CURRENT_VERSION.md","stage16_artifact_manifest.py"],
}
rows=[]
seen=set()
for role,names in groups.items():
    for name in names:
        if name in seen: continue
        seen.add(name)
        file=root/name
        if not file.is_file(): raise FileNotFoundError(name)
        data=file.read_bytes()
        rows.append({"role":role,"path":name.replace("\\","/"),
                     "bytes":len(data),"sha256":hashlib.sha256(data).hexdigest()})
rows.sort(key=lambda r:r["path"])
with (root/"artifact_manifest.csv").open("w",newline="",encoding="utf-8") as f:
    writer=csv.DictWriter(f,fieldnames=("role","path","bytes","sha256"))
    writer.writeheader();writer.writerows(rows)
print(f"{len(rows)} artifacts hashed; no duplicate paths")
