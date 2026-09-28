"""Stage 0 certificate replays in a sandbox copy (original scripts unmodified).

Usage:  python tools/replay_sandbox.py prepare
        python tools/replay_sandbox.py run            (runs all chains in parallel, resumable)
        python tools/replay_sandbox.py compare
The sandbox lives in stage0_baseline/replay_sandbox.  Every chain runs with
cwd = sandbox, so relative-path writes (e.g. stage12_cover_audit.py) stay inside it.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys
import threading
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
ROOT = HERE.parent.parent
SB = HERE / "stage0_baseline" / "replay_sandbox"
LOGS = HERE / "stage0_baseline" / "replay_logs"
STATE = HERE / "stage0_baseline" / "replay_state.json"
PY = sys.executable
SIZE_LIMIT = 50 * 1024 * 1024

CHAINS = {
    "R1_arb_finite_z": dict(
        cmds=[["stage12_full_arb_replay.py", "21"], ["stage12_full_arb_replay.py", "31"], ["stage12_cover_audit.py"]],
        outputs=["stage12_replay/N21_predicates.csv", "stage12_replay/N21_refinement_leaves.csv",
                 "stage12_replay/N21_summary.json", "stage12_replay/N31_predicates.csv",
                 "stage12_replay/N31_refinement_leaves.csv", "stage12_replay/N31_summary.json",
                 "stage12_replay/coverage_audit.json", "stage12_replay/refinement_coverage_audit.json"]),
    "R2_generalN": dict(
        cmds=[["stage12_generalN.py", str(n)] for n in (11, 15, 21, 31, 41)],
        outputs=[f"stage12_generalN/N{n}_{k}" for n in (11, 15, 21, 31, 41)
                 for k in ("certificate.json", "partition.csv", "reference.json", "roots.csv")]),
    "R3_continuum": dict(
        cmds=[["stage14_continuum.py"], ["stage14_continuum_interval.py"], ["stage14_continuum_global_hybrid.py"],
              ["stage14_replay_global.py"], ["stage14_continuum_endpoint_interval.py"]],
        outputs=["stage14_continuum_candidate.json", "stage14_continuum_candidate_roots.csv",
                 "stage14_continuum_local_interval.json", "continuum_global_certificate.json",
                 "stage14_continuum_global_partition.csv", "stage14_continuum_global_roots.csv",
                 "stage14_continuum_global_replay.json", "stage14_continuum_endpoint_interval.json"]),
    "R4_independent_arb_N0": dict(
        cmds=[["stage16_independent_arb.py"]],
        outputs=["independent_N0_replay.json", "independent_N0_replay.csv"], flint_path=True),
    "R5_largeN_floor_10001": dict(
        cmds=[["stage17_uniform_trial.py", "10001"], ["stage17_refine_uniform.py", "10001"],
              ["stage17_largeN_margin_summary.py"], ["stage17_write_largeN.py"],
              ["stage17_annotate_condition_thresholds.py"]],
        outputs=["stage17_interval_trial_N10001.json", "stage17_refinement_N10001.json",
                 "stage17_largeN_margins.json", "tightened_N0_certificate.json", "largeN_condition_thresholds.csv"]),
    "R6_finite_width_eps_0p035": dict(
        cmds=[["stage17_epsilon_refine.py", "0.035", "14287"], ["stage17_epsilon_root_tiling.py", "0.035", "4096"],
              ["stage17_epsilon_inner_tiled.py", "_0p035", "3000000"], ["stage17_epsilon_joints.py", "_0p035"],
              ["stage17_epsilon_coverage_audit.py", "0.035"], ["stage17_write_epsilon_certificate.py", "0.035"],
              ["stage17_epsilon_coeff_crosscheck.py"]],
        outputs=["stage17_epsilon_outer_refinement_0p035_14287.json", "stage17_epsilon_root_tiling_0p035.json",
                 "stage17_epsilon_inner_tiled_0p035.json", "stage17_epsilon_joints_0p035.json",
                 "stage17_epsilon_coverage_audit_0p035.json", "finite_width_epsilon_certificate_0p035.json",
                 "finite_width_epsilon_certificate_0p035.md", "epsilon_interval_obligations_0p035.csv",
                 "stage17_epsilon_coeff_crosscheck.json"]),
    "R7_lower_transition": dict(
        cmds=[["stage13_build.py"]],
        outputs=["lower_transition_certificate.json", "tangent_cost_vs_lambda.csv", "smallz_global_regimes.csv"]),
}


def sha(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda: f.read(8 << 20), b""):
            h.update(b)
    return h.hexdigest()


def prepare():
    if SB.exists():
        shutil.rmtree(SB)
    SB.mkdir(parents=True)
    n = 0
    for p in ROOT.iterdir():
        if p.is_file() and p.stat().st_size <= SIZE_LIMIT and p.suffix.lower() != ".zip":
            shutil.copy2(p, SB / p.name)
            n += 1
    for d in ("stage12_deps",):
        shutil.copytree(ROOT / d, SB / d)
    for d in ("stage12_replay", "stage12_generalN", "paper/figures"):
        (SB / d).mkdir(parents=True, exist_ok=True)
    # inputs needed from subdirectories
    for rel in ("stage12_generalN",):
        for p in (ROOT / rel).glob("*"):
            if p.is_file():
                shutil.copy2(p, SB / rel / p.name)
    for rel in ("stage12_replay",):
        for p in (ROOT / rel).glob("*"):
            if p.is_file():
                shutil.copy2(p, SB / rel / p.name)
    # remove every chain output so a crash cannot masquerade as a reproduction
    removed = []
    for c in CHAINS.values():
        for o in c["outputs"]:
            q = SB / o
            if q.exists():
                q.unlink()
                removed.append(o)
    json.dump(dict(prepared=time.strftime("%Y-%m-%dT%H:%M:%S%z"), copied_files=n, removed_outputs=removed),
              open(SB / "_SANDBOX_INFO.json", "w"), indent=1)
    print(f"sandbox prepared: {n} files copied, {len(removed)} outputs removed")


def _state():
    return json.loads(STATE.read_text()) if STATE.exists() else {}


_lock = threading.Lock()


def _save(name, rec):
    with _lock:
        st = _state()
        st[name] = rec
        STATE.write_text(json.dumps(st, indent=1))


def run_chain(name):
    c = CHAINS[name]
    st = _state().get(name, {})
    if st.get("status") == "done":
        return
    LOGS.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    env.update(OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", PYTHONUNBUFFERED="1",
               MPLBACKEND="Agg")
    if c.get("flint_path"):
        env["PYTHONPATH"] = str(SB / "stage12_deps")
    rec = dict(status="running", started=time.strftime("%Y-%m-%dT%H:%M:%S%z"), steps=[])
    _save(name, rec)
    ok = True
    for cmd in c["cmds"]:
        t0 = time.time()
        logf = LOGS / f"{name}__{'_'.join(cmd).replace('.py', '').replace('.', 'p')}.log"
        with logf.open("w", encoding="utf-8") as lf:
            r = subprocess.run([PY, "-u"] + cmd, cwd=SB, env=env, stdout=lf, stderr=subprocess.STDOUT)
        rec["steps"].append(dict(cmd=cmd, returncode=r.returncode, seconds=round(time.time() - t0, 1), log=logf.name))
        _save(name, rec)
        if r.returncode != 0:
            ok = False
            break
    rec["status"] = "done" if ok else "failed"
    rec["finished"] = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    _save(name, rec)


DEPS = {"R4_independent_arb_N0": ["R3_continuum"], "R5_largeN_floor_10001": ["R3_continuum"]}


def _run_with_deps(name):
    for d in DEPS.get(name, []):
        while _state().get(d, {}).get("status") not in ("done", "failed"):
            time.sleep(5)
        if _state().get(d, {}).get("status") != "done":
            _save(name, dict(status="skipped_dependency_failed", dependency=d))
            return
    run_chain(name)


def run(names=None):
    names = names or list(CHAINS)
    th = [threading.Thread(target=_run_with_deps, args=(n,)) for n in names]
    for t in th:
        t.start()
    for t in th:
        t.join()


# ----------------------------------------------------------------------------- comparison
VOLATILE = ("second", "elapsed", "runtime", "time", "timestamp", "date", "prepared", "started", "finished",
            "path", "python", "platform", "version")


def _flatten(o, p=""):
    out = {}
    if isinstance(o, dict):
        for k, v in o.items():
            out.update(_flatten(v, f"{p}.{k}" if p else str(k)))
    elif isinstance(o, list):
        for i, v in enumerate(o):
            out.update(_flatten(v, f"{p}[{i}]"))
    else:
        out[p] = o
    return out


def compare():
    st = _state()
    report = {}
    for name, c in CHAINS.items():
        rows = []
        for o in c["outputs"]:
            a, b = ROOT / o, SB / o
            r = dict(output=o, frozen_exists=a.exists(), replay_exists=b.exists())
            if a.exists() and b.exists():
                r["identical_bytes"] = sha(a) == sha(b)
                if not r["identical_bytes"] and o.endswith(".json"):
                    fa, fb = _flatten(json.loads(a.read_text(encoding="utf-8"))), _flatten(json.loads(b.read_text(encoding="utf-8")))
                    diffs = [k for k in set(fa) | set(fb) if fa.get(k) != fb.get(k)]
                    vol = [k for k in diffs if any(t in k.lower() for t in VOLATILE)]
                    sub = [k for k in diffs if k not in vol]
                    r["json_diff_keys_total"] = len(diffs)
                    r["json_diff_volatile"] = len(vol)
                    r["json_diff_substantive"] = sorted(sub)[:40]
                    r["json_diff_substantive_count"] = len(sub)
                    flags = [k for k in set(fa) if any(t in k.lower() for t in ("verified", "pass", "all_", "complete"))
                             and isinstance(fa.get(k), bool)]
                    r["flag_mismatches"] = sorted(k for k in flags if fa.get(k) != fb.get(k))
                elif not r["identical_bytes"] and o.endswith(".csv"):
                    ra = list(csv.reader(a.open(encoding="utf-8", newline="")))
                    rb = list(csv.reader(b.open(encoding="utf-8", newline="")))
                    r["rows_frozen"], r["rows_replay"] = len(ra), len(rb)
                    r["rows_identical"] = sum(x == y for x, y in zip(ra, rb))
            rows.append(r)
        report[name] = dict(state=st.get(name), outputs=rows)
    (HERE / "stage0_baseline" / "replay_comparison.json").write_text(json.dumps(report, indent=1))
    print(json.dumps({k: dict(status=(v["state"] or {}).get("status"),
                              identical=sum(1 for r in v["outputs"] if r.get("identical_bytes")),
                              n=len(v["outputs"])) for k, v in report.items()}, indent=1))


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "run"
    if cmd == "prepare":
        prepare()
    elif cmd == "run":
        run(sys.argv[2:] or None)
    elif cmd == "compare":
        compare()
