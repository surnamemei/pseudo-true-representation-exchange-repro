"""Replay the frozen interval scripts in an isolated directory.

The source certificate files in the project root are never opened for writing.
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
BASE = ROOT / "stage10_replay"

N21 = [
    "three_to_two_tone_stage2_reference.py",
    "three_to_two_tone_stage2_krawczyk.py",
    "stage10_replay_n21_large.py",
    "three_to_two_tone_stage2_crossing.py",
    "three_to_two_tone_stage2_outer_global.py",
    "three_to_two_tone_stage2_inner_global.py",
    "three_to_two_tone_stage2_cell_audit.py",
    "three_to_two_tone_stage2_crossing_global.py",
    "three_to_two_tone_stage2_certificate_build.py",
]
N31 = ["stage7_n31_certificate.py", "stage7_n31_inner_log.py"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case", choices=["n21", "n31"])
    args = ap.parse_args()
    dest = BASE / args.case
    dest.mkdir(parents=True, exist_ok=True)
    if not dest.resolve().is_relative_to(BASE.resolve()):
        raise RuntimeError("Replay path escaped the workspace replay directory")
    for path in ROOT.glob("*.py"):
        if path.name != Path(__file__).name:
            shutil.copy2(path, dest / path.name)
    if args.case == "n31":
        # The shared inner-domain function is imported from the Stage-2 module,
        # whose module-level constants read these N=21 logs. N=31 supplies its
        # own roots and interval predicates to the function at call time.
        for name in ("stage2_reference.json", "stage2_local_interval_checks.json",
                     "stage3_tangent_reference.json"):
            shutil.copy2(ROOT / name, dest / name)
    scripts = N21 if args.case == "n21" else N31
    emitted = (
        ["stage2_reference.json", "stage2_local_interval_checks.json",
         "stage2_large_krawczyk_checks.json", "stage2_crossing_interval.json",
         "crossing_certificate.csv", "stage2_outer_global_summary.json",
         "interval_boxes.csv", "stage2_inner_global_summary.json",
         "stage2_inner_cells.csv", "stage2_cell_audit.json",
         "stage2_crossing_global.json", "validated_globality_certificate.json"]
        if args.case == "n21" else
        ["N31_crossing_certificate.csv", "N31_interval_cells.csv",
         "N31_global_certificate.json", "N31_inner_interval_cells.csv"]
    )
    for name in emitted:
        (dest/name).unlink(missing_ok=True)
    log = dest / "replay_log.txt"
    with log.open("w", encoding="utf-8") as f:
        for script in scripts:
            print("RUN", script, flush=True)
            f.write("RUN " + script + "\n")
            f.flush()
            p = subprocess.run([sys.executable, "-u", script], cwd=dest,
                               stdout=f, stderr=subprocess.STDOUT)
            print("EXIT", script, p.returncode, flush=True)
            f.write("EXIT " + script + " " + str(p.returncode) + "\n")
            f.flush()
            if p.returncode:
                return p.returncode
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
