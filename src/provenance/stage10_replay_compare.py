"""Compare isolated Stage-10 replays with frozen certificate artifacts."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parent
RPL=ROOT/"stage10_replay"
OUT=ROOT/"certificate_replay_comparison.csv"

files={
"n21":[
"stage2_reference.json","stage2_local_interval_checks.json",
"stage2_large_krawczyk_checks.json",
"stage2_crossing_interval.json","crossing_certificate.csv",
"stage2_outer_global_summary.json","interval_boxes.csv",
"stage2_inner_global_summary.json","stage2_inner_cells.csv",
"stage2_cell_audit.json","stage2_crossing_global.json",
"validated_globality_certificate.json"],
"n31":["N31_crossing_certificate.csv","N31_interval_cells.csv",
       "N31_global_certificate.json","N31_inner_interval_cells.csv"],
}
rows=[]
def add(case,art,metric,frozen,replay,match,note=""):
    rows.append(dict(case=case,artifact=art,metric=metric,frozen=str(frozen),
                     replay=str(replay),match=str(match).lower(),notes=note))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest().upper()
def nrows(p):
    with p.open(encoding="utf-8",newline="") as f:return sum(1 for _ in f)-1
def cmp(case,art,metric,a,b):add(case,art,metric,json.dumps(a,sort_keys=True),json.dumps(b,sort_keys=True),a==b)

for case,names in files.items():
    for name in names:
        p=ROOT/name;q=RPL/case/name
        if not q.exists():
            add(case,name,"file",sha(p),"MISSING",False,"Replay incomplete or script did not emit artifact")
            continue
        add(case,name,"sha256",sha(p),sha(q),sha(p)==sha(q))
        if p.suffix==".csv":add(case,name,"data_rows",nrows(p),nrows(q),nrows(p)==nrows(q))

for case,name in (("n21","validated_globality_certificate.json"),
                  ("n31","N31_global_certificate.json")):
    q=RPL/case/name
    if not q.exists():continue
    a=json.loads((ROOT/name).read_text());b=json.loads(q.read_text())
    add(case,name,"parsed_JSON_all_fields_equal",True,a==b,a==b,
        "Exact parsed structures, including every root and inequality")
    canon=lambda d:hashlib.sha256(json.dumps(d,sort_keys=True,separators=(",",":"),ensure_ascii=True).encode()).hexdigest().upper()
    add(case,name,"canonical_JSON_sha256",canon(a),canon(b),canon(a)==canon(b),
        "Key ordering and whitespace ignored")
    if case=="n21":
        for i,label in enumerate(("below","crossing","above")):
            for key in ("epsilon","A_box","B_box","A_Krawczyk","B_Krawczyk","A_H_positive","B_H_positive"):
                cmp(case,name,f"{label}.{key}",a["cases"][i][key],b["cases"][i][key])
            for key in ("terminal_cells","excluded_cells","accepted_cells","failed_interval_cells"):
                cmp(case,name,f"{label}.outer.{key}",a["cases"][i]["outer"][key],b["cases"][i]["outer"][key])
        for key in ("epsilon_bracket","outer_excluded_cells","outer_failed","outer_min_margin","full_bracket_global"):
            cmp(case,name,f"crossing_global.{key}",a["crossing_global_bracket"][key],b["crossing_global_bracket"][key])
        cmp(case,name,"verified",a["verified"],b["verified"])
    else:
        for key in ("high_precision_numerical_root","epsilon_bracket","crossing_checks","large_root_boxes","slope_interval","crossing_bracket_certified","outer","inner","full_global_certificate"):
            cmp(case,name,key,a[key],b[key])

with OUT.open("w",encoding="utf-8",newline="") as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print("wrote",OUT,"rows",len(rows),"failed",sum(r["match"]!="true" for r in rows))
