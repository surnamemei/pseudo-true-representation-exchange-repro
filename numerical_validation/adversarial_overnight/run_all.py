"""Master resumable runner for the adversarial overnight validation.

    python -u run_all.py                 # all stages in priority order (0,1,2,3,6,4,5,7,8)
    python -u run_all.py --stages 0,1    # subset
Each stage writes stage*/DONE.json on completion and is skipped on rerun; inside a
stage, per-configuration scan/analysis files are checkpoints, so an interrupted
run resumes where it stopped.  STATUS.json is updated after every major unit.
"""
from __future__ import annotations

import os

for _k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_k, "1")
os.environ.setdefault("MKL_CBWR", "COMPATIBLE")
os.environ.setdefault("MPLBACKEND", "Agg")

import argparse
import importlib
import sys
import time
import traceback
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from advrun import log, status_update, status_event, status_load, jdump, code_hashes, shutdown_pool, sha  # noqa: E402

PRIORITY = ["0", "1", "2", "3", "6", "4", "5", "7", "8"]
MODULES = {"0": "s0_baseline", "1": "s1_window", "2": "s2_phase", "3": "s3_bwidth", "4": "s4_imbalance",
           "5": "s5_matching", "6": "s6_noise", "7": "s7_numinv", "8": "s8_real"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stages", default=",".join(PRIORITY))
    args = ap.parse_args()
    st = status_load()
    if "prereg_sha256" not in st:
        status_update(prereg_sha256=sha(HERE / "PREREGISTRATION.md"), freeze_sha256=sha(HERE / "FREEZE.md"))
    status_update(runner_started=time.strftime("%Y-%m-%dT%H:%M:%S%z"), runner_pid=os.getpid(),
                  code_sha256=code_hashes())
    status_event(f"runner start stages={args.stages}")
    for sid in args.stages.split(","):
        sid = sid.strip()
        mod = MODULES.get(sid)
        if not mod or not (HERE / f"{mod}.py").exists():
            log(f"stage {sid}: module {mod} not present yet; skipped", sid)
            continue
        t0 = time.time()
        log(f"=== stage {sid} ({mod}) start", sid)
        status_update(sid, started=time.strftime("%Y-%m-%dT%H:%M:%S%z"), module=mod)
        try:
            m = importlib.import_module(mod)
            m.run()
            status_update(sid, finished=time.strftime("%Y-%m-%dT%H:%M:%S%z"), seconds=round(time.time() - t0, 1))
            log(f"=== stage {sid} done in {time.time() - t0:.1f} s", sid)
        except Exception as exc:
            log(f"stage {sid} FAILED: {exc!r}\n{traceback.format_exc()}", sid)
            status_update(sid, status="error", error=repr(exc))
            status_event(f"stage {sid} error {exc!r}")
    shutdown_pool()
    status_event("runner end")
    log("runner finished")


if __name__ == "__main__":
    main()
