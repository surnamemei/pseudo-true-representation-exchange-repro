"""Recreate the archived N=21 larger Krawczyk root-box log.

The original archive did not keep the one-shot driver that selected the two
larger radii. This calls the original interval predicate, with the radii
recorded in the frozen log, at each of the original three amplitudes.
"""
import json
from pathlib import Path
from three_to_two_tone_stage2_krawczyk import check, ref

ROOT=Path(__file__).resolve().parent
rows=[]
for rec in ref["records"]:
    for which,radius in (("A","0.0003"),("B","0.01")):
        rows.append(check(rec["label"],rec["epsilon"],which,radius=radius))
out=ROOT/"stage2_large_krawczyk_checks.json"
out.write_text(json.dumps(rows,indent=2),encoding="utf-8")
print("regenerated",out.name,"rows",len(rows),"all included",
      all(r["krawczyk_inclusion"] and r["positive_definite"] for r in rows))
