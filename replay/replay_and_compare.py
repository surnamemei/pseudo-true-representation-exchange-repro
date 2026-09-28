#!/usr/bin/env python3
"""Replay the seven theorem-level certificate chains (R1-R7) in a scratch workspace and compare every
output with the archived frozen record.

Usage, from the directory that contains this tool's parent (the archive or repository root):

    python <tools>/replay_and_compare.py prepare [--work DIR]
    python <tools>/replay_and_compare.py run     [--work DIR] [R1 R2 ... T1] # default: all chains, in parallel
    python <tools>/replay_and_compare.py compare [--work DIR]
    python <tools>/replay_and_compare.py all     [--work DIR] [R1 R2 ...]   # prepare + run + compare
    python <tools>/replay_and_compare.py figures [--work DIR]              # redraw the manuscript figures

The workspace (default ./_replay_work) is assembled from workspace_map.csv, which sits next to this tool
and lists `workspace_path,source_path` pairs. The frozen scripts run unmodified, one Python process per
command, with the workspace as the working directory. Every declared chain output is deleted from the
workspace before the run, so a crash cannot masquerade as a reproduction. The archived files are never
written.

Comparison rules. An output passes when it is byte-identical to the archived file, or when the only
differences are (a) JSON key order, (b) timing, path or version fields, (c) SHA-256 strings of
regenerated intermediate files, or (d) numeric enclosure values whose relative difference is at most
NUMERIC_TOL. Floating-point preconditioners and candidate centres can differ between platforms and
BLAS builds; they move enclosure bounds in far digits but must not change any verdict. Booleans,
predicate names, labels and integer counts must match exactly. Outputs too large to archive are compared
by SHA-256 with frozen_output_hashes.json, both as frozen (CRLF line endings written on Windows) and after
CRLF->LF conversion.

The chain definitions R1-R7 are those of validation/adversarial_overnight/tools/replay_sandbox.py, the
Stage-0 tool used for the manuscript's replay record.

Chain T1 (added with the approved Stage-25 revision) recomputes the deterministic part of the sealed
first-order transition theory (paper/transition_theory/): the branch quantities, the envelope slope and its
checks, and every predicted selection probability and global error. The one-shot validation against the
Monte Carlo is not a chain and is never rerun. For T1, fields that are roundoff-level diagnostics of the
checks (finite-difference residuals, synthetic-sample moments) are reported but not held to NUMERIC_TOL.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import threading
import time
from decimal import Decimal, getcontext
from pathlib import Path

getcontext().prec = 200
TOOL_DIR = Path(__file__).resolve().parent
BASE = TOOL_DIR.parent
NUMERIC_TOL = Decimal("1e-10")

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
        outputs=["independent_N0_replay.json", "independent_N0_replay.csv"]),
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
    "T1_transition_theory": dict(
        cmds=[["paper/transition_theory/code/transition_theory.py"],
              ["paper/transition_theory/code/second_order_diagnostic.py"]],
        outputs=["results/transition_theory/deterministic/predictions.csv",
                 "results/transition_theory/deterministic/deterministic.json",
                 "results/transition_theory/deterministic/checks.json",
                 "results/transition_theory/deterministic/second_order_diagnostic.json"],
        diagnostic=("rel_diff", "abs_diff", "g_fd", "fd_check", "sample_", "var_ratio", "mean_over_sd", "delta0_resid",
                    "grad_norm", "environment", "seconds")),
}
# Declared outputs that the figure script reads: `figures` restores the archived copy when a prepared
# workspace lacks it (Fig. 7 draws the sealed predictions).
FIGURE_INPUTS = ("results/transition_theory/deterministic/predictions.csv",)
DEPS = {"R4_independent_arb_N0": ["R3_continuum"], "R5_largeN_floor_10001": ["R3_continuum"]}
SHORT = {k.split("_")[0]: k for k in CHAINS}
VOLATILE = ("second", "elapsed", "runtime", "time", "timestamp", "date", "prepared", "started", "finished",
            "path", "python", "platform", "version")
NUM = re.compile(r"^\s*[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?\s*$")
IVL = re.compile(r"^\s*\[\s*([^,\]]+)\s*,\s*([^\]]+)\]\s*$")


# ----------------------------------------------------------------------------- workspace
def workspace_map():
    with (TOOL_DIR / "workspace_map.csv").open(encoding="utf-8", newline="") as fh:
        return [(r["workspace_path"], r["source_path"]) for r in csv.DictReader(fh)]


def prepare(work: Path):
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    n = 0
    for wpath, spath in workspace_map():
        src = BASE / spath
        if not src.is_file():
            raise SystemExit(f"missing source file listed in workspace_map.csv: {spath}")
        dst = work / wpath
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        n += 1
    for d in ("stage12_replay", "stage12_generalN", "paper/figures", "supplement/figures"):
        (work / d).mkdir(parents=True, exist_ok=True)       # output directories the frozen scripts expect
    removed = []
    for c in CHAINS.values():
        for o in c["outputs"]:
            q = work / o
            if q.exists():
                q.unlink()
                removed.append(o)
    (work / "_workspace_info.json").write_text(json.dumps(dict(
        prepared=time.strftime("%Y-%m-%dT%H:%M:%S%z"), files_copied=n, outputs_removed=removed), indent=1))
    print(f"workspace {work}: {n} files copied, {len(removed)} declared outputs removed")


# ----------------------------------------------------------------------------- run
_lock = threading.Lock()


def _state(work):
    p = work / "_state.json"
    for _ in range(50):
        try:
            return json.loads(p.read_text()) if p.exists() else {}
        except json.JSONDecodeError:
            time.sleep(0.1)
    raise RuntimeError("state file unreadable")


def _save(work, name, rec):
    with _lock:
        st = _state(work)
        st[name] = rec
        tmp = work / "_state.tmp"
        tmp.write_text(json.dumps(st, indent=1))
        os.replace(tmp, work / "_state.json")


def run_chain(work, name):
    logs = work / "_logs"
    logs.mkdir(exist_ok=True)
    env = dict(os.environ)
    env.update(OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", PYTHONUNBUFFERED="1",
               MPLBACKEND="Agg", PYTHONDONTWRITEBYTECODE="1")
    rec = dict(status="running", started=time.strftime("%Y-%m-%dT%H:%M:%S%z"), steps=[])
    _save(work, name, rec)
    for o in CHAINS[name]["outputs"]:          # a crash must never leave an older copy to be compared
        (work / o).unlink(missing_ok=True)
    ok = True
    for cmd in CHAINS[name]["cmds"]:
        tag = f"{name}__{'_'.join(cmd).replace('.py', '').replace('.', 'p').replace('/', '_')}"
        t0 = time.time()
        with (logs / f"{tag}.log").open("w", encoding="utf-8") as lf:
            p = subprocess.Popen([sys.executable, "-u"] + cmd, cwd=work, env=env, stdout=lf, stderr=subprocess.STDOUT)
            if hasattr(os, "wait4"):                      # POSIX: exit status and peak memory of this command
                _, status, usage = os.wait4(p.pid, 0)
                p.returncode = os.waitstatus_to_exitcode(status)
                peak = round(usage.ru_maxrss / (1024 * 1024 if sys.platform == "darwin" else 1024), 1)  # bytes on macOS, KiB on Linux
            else:
                p.wait()
                peak = None
        step = dict(cmd=cmd, returncode=p.returncode, seconds=round(time.time() - t0, 1), log=f"{tag}.log")
        if peak is not None:
            step["peak_rss_mib"] = peak
        rec["steps"].append(step)
        _save(work, name, rec)
        if p.returncode != 0:
            ok = False
            break
    rec["status"] = "done" if ok else "failed"
    rec["finished"] = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    _save(work, name, rec)


def _run_with_deps(work, name, selected):
    for d in DEPS.get(name, []):
        if d not in selected and _state(work).get(d, {}).get("status") != "done":
            _save(work, name, dict(status="skipped", reason=f"requires {d} in the same workspace"))
            return
        while _state(work).get(d, {}).get("status") not in ("done", "failed", "skipped"):
            time.sleep(5)
        if _state(work).get(d, {}).get("status") != "done":
            _save(work, name, dict(status="skipped", reason=f"dependency {d} did not complete"))
            return
    run_chain(work, name)


def run(work, names):
    if not (work / "_workspace_info.json").exists():
        raise SystemExit("run 'prepare' first")
    th = [threading.Thread(target=_run_with_deps, args=(work, n, set(names))) for n in names]
    for t in th:
        t.start()
    for t in th:
        t.join()
    for n in names:
        s = _state(work).get(n, {})
        steps = ", ".join(f"{' '.join(x['cmd'])} -> exit {x['returncode']} ({x['seconds']} s)" for x in s.get("steps", []))
        print(f"{n}: {s.get('status')}  [{steps}]")


# ----------------------------------------------------------------------------- compare
def _nums(v):
    if isinstance(v, bool) or v is None or isinstance(v, int):
        return None
    if isinstance(v, float):
        return [Decimal(repr(v))]
    if isinstance(v, str):
        if NUM.match(v) and not re.fullmatch(r"\s*[-+]?\d+\s*", v):
            return [Decimal(v)]
        m = IVL.match(v)
        if m and NUM.match(m.group(1)) and NUM.match(m.group(2)):
            return [Decimal(m.group(1)), Decimal(m.group(2))]
    return None


def _classify(pairs, diagnostic=()):
    res = dict(numeric=0, hash=0, volatile=0, diagnostic=0, other=0, max_rel=Decimal(0), worst=None, other_examples=[])
    for key, a, b in pairs:
        if a == b:
            continue
        if any(t in key.lower() for t in VOLATILE):
            res["volatile"] += 1
            continue
        if any(t in key.lower() for t in diagnostic):
            res["diagnostic"] += 1
            continue
        if isinstance(a, str) and isinstance(b, str) and re.fullmatch(r"[0-9a-f]{64}", a) and re.fullmatch(r"[0-9a-f]{64}", b):
            res["hash"] += 1
            continue
        na, nb = _nums(a), _nums(b)
        if na is not None and nb is not None and len(na) == len(nb):
            res["numeric"] += 1
            for x, y in zip(na, nb):
                s = max(abs(x), abs(y))
                r = abs(x - y) / s if s else abs(x - y)
                if r > res["max_rel"]:
                    res["max_rel"], res["worst"] = r, [key, str(a)[:80], str(b)[:80]]
            continue
        res["other"] += 1
        if len(res["other_examples"]) < 5:
            res["other_examples"].append([key, str(a)[:100], str(b)[:100]])
    return res


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


def _sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(8 << 20), b""):
            h.update(b)
    return h.hexdigest()


def _sha_lf(p):
    h, carry = hashlib.sha256(), b""
    with open(p, "rb") as f:
        for block in iter(lambda: f.read(8 << 20), b""):
            block = carry + block
            carry = b"\r" if block.endswith(b"\r") else b""
            if carry:
                block = block[:-1]
            h.update(block.replace(b"\r\n", b"\n"))
    h.update(carry)
    return h.hexdigest()


def compare_file(frozen: Path | None, replay: Path, frozen_hash: dict | None, diagnostic=()):
    if not replay.exists():
        return dict(verdict="FAIL", reason="output missing after replay")
    if frozen is None or not frozen.exists():
        if frozen_hash:
            h = _sha(replay)
            if h == frozen_hash.get("sha256"):
                return dict(verdict="PASS", reason="not archived; SHA-256 equals the frozen hash")
            if _sha_lf(replay) == frozen_hash.get("sha256_lf"):
                return dict(verdict="PASS", reason="not archived; SHA-256 after CRLF->LF equals the frozen file's")
            return dict(verdict="UNVERIFIED", reason="not archived and its SHA-256 differs from the frozen hash; "
                        "a content comparison needs the frozen file (the downstream certificate records are compared)")
        return dict(verdict="UNVERIFIED", reason="no archived copy")
    if _sha(frozen) == _sha(replay):
        return dict(verdict="PASS", reason="byte-identical")
    big = frozen.stat().st_size > 150e6
    if not big and frozen.read_bytes().replace(b"\r\n", b"\n") == replay.read_bytes().replace(b"\r\n", b"\n"):
        return dict(verdict="PASS", reason="identical after normalizing line endings (CRLF from Windows text mode)")
    if frozen.suffix == ".json":
        if big:                                     # indent=2 layout: compare key/value lines, newlines normalized
            pairs = []
            with open(frozen) as fa, open(replay) as fb:
                for i, (la, lb) in enumerate(zip(fa, fb)):
                    if la != lb:
                        split = lambda s: s.split('":', 1)[1] if '":' in s else s
                        key = la.split('":', 1)[0].strip().strip('"') if '":' in la else f"line{i}"
                        pairs.append((key, split(la).strip().rstrip(",").strip('"'), split(lb).strip().rstrip(",").strip('"')))
                if next(fa, None) is not None or next(fb, None) is not None:
                    return dict(verdict="FAIL", reason="line count differs")
            if not pairs:
                return dict(verdict="PASS", reason="identical after normalizing line endings (CRLF from Windows text mode)")
            res = _classify(pairs, diagnostic)
        else:
            fa = _flatten(json.loads(frozen.read_text()))
            fb = _flatten(json.loads(replay.read_text()))
            if fa == fb:
                return dict(verdict="PASS", reason="same content; JSON key order differs")
            if set(fa) != set(fb):
                return dict(verdict="FAIL", reason="JSON keys differ",
                            keys=sorted(set(fa) ^ set(fb))[:10])
            res = _classify(((k, fa[k], fb[k]) for k in sorted(fa)), diagnostic)
    elif frozen.suffix == ".csv":
        ra = list(csv.reader(open(frozen, encoding="utf-8", newline="")))
        rb = list(csv.reader(open(replay, encoding="utf-8", newline="")))
        if len(ra) != len(rb) or any(len(x) != len(y) for x, y in zip(ra, rb)):
            return dict(verdict="FAIL", reason="CSV shape differs")
        hdr = ra[0] if ra else []
        res = _classify(((f"row{i}.{hdr[j] if j < len(hdr) else j}", u, v)
                         for i, (x, y) in enumerate(zip(ra, rb)) for j, (u, v) in enumerate(zip(x, y))), diagnostic)
    elif frozen.suffix == ".md":
        la, lb = frozen.read_text().splitlines(), replay.read_text().splitlines()
        if len(la) != len(lb):
            return dict(verdict="FAIL", reason="line count differs")
        pairs = []
        for i, (x, y) in enumerate(zip(la, lb)):
            tx, ty = re.split(r"(\s+|[(),;=`|])", x), re.split(r"(\s+|[(),;=`|])", y)
            if len(tx) != len(ty):
                pairs.append((f"line{i}", x, y))
            else:
                pairs += [(f"line{i}.tok{k}", u, v) for k, (u, v) in enumerate(zip(tx, ty))]
        res = _classify(pairs, diagnostic)
    else:
        return dict(verdict="UNVERIFIED", reason="binary file differs")
    ok = res["other"] == 0 and res["max_rel"] <= NUMERIC_TOL
    return dict(verdict="PASS" if ok else "FAIL",
                reason=("differences confined to " + ", ".join(k for k in ("numeric", "hash", "volatile", "diagnostic") if res[k])
                        if ok else "non-numeric difference or numeric difference above tolerance"),
                numeric_values_differing=res["numeric"], max_relative_numeric_difference=f"{res['max_rel']:.3e}",
                worst_numeric=res["worst"], hash_fields=res["hash"], volatile_fields=res["volatile"],
                diagnostic_fields=res["diagnostic"],
                other_differences=res["other"], other_examples=res["other_examples"])


def compare(work, names=None):
    src_of = {w: s for w, s in workspace_map()}
    hashes_p = TOOL_DIR / "frozen_output_hashes.json"
    hashes = json.loads(hashes_p.read_text()) if hashes_p.exists() else {}
    state = _state(work)
    judged = set(names) if names else ({n for n in CHAINS if n in state} or set(CHAINS))
    report, all_pass = {}, True
    for name, c in CHAINS.items():
        st = state.get(name, {})
        if name not in judged:
            report[name] = dict(status=st.get("status"), verdict="NOT REQUESTED", outputs={})
            continue
        rows = {}
        for o in c["outputs"]:
            frozen = BASE / src_of[o] if o in src_of else None
            rows[o] = compare_file(frozen, work / o, hashes.get(o), c.get("diagnostic", ()))
        exit_ok = st.get("status") == "done"
        chain_pass = exit_ok and all(r["verdict"] == "PASS" for r in rows.values())
        all_pass &= chain_pass
        report[name] = dict(status=st.get("status"), steps=st.get("steps", []),
                            verdict="PASS" if chain_pass else ("NOT RUN" if not st else "REVIEW"), outputs=rows)
    (work / "_comparison.json").write_text(json.dumps(report, indent=1))
    for name, r in report.items():
        if r["verdict"] == "NOT REQUESTED":
            print(f"{name}: not requested")
            continue
        n_ok = sum(1 for v in r["outputs"].values() if v["verdict"] == "PASS")
        print(f"{name}: {r['verdict']} (run status {r['status']}; outputs passing {n_ok}/{len(r['outputs'])})")
        for o, v in r["outputs"].items():
            if v["verdict"] != "PASS":
                print(f"    {o}: {v['verdict']} - {v['reason']}")
    print("OVERALL:", ("PASS" if all_pass else "REVIEW (see _comparison.json)") + f" for {len(judged)} chain(s)")
    return all_pass


def figures(work):
    """Redraw the manuscript figures from the archived results (no new search or solver run)."""
    script = work / "paper/figure_scripts/make_revision_figures.py"
    if not script.exists():
        raise SystemExit("run 'prepare' first (the figure script is copied into the workspace)")
    for d in ("paper/figures", "supplement/figures"):
        (work / d).mkdir(parents=True, exist_ok=True)
    src_of = {w: s for w, s in workspace_map()}
    for rel in FIGURE_INPUTS:
        if not (work / rel).exists() and rel in src_of:
            (work / rel).parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(BASE / src_of[rel], work / rel)
            print(f"restored the archived {rel} for the figure script")
    env = dict(os.environ, MPLBACKEND="Agg", PYTHONDONTWRITEBYTECODE="1")
    r = subprocess.run([sys.executable, str(script.relative_to(work))], cwd=work, env=env)
    outs = sorted(str(p.relative_to(work)) for d in ("paper/figures", "supplement/figures")
                  for p in (work / d).glob("fig*.pdf"))
    print("figure script exit", r.returncode, "; outputs:", ", ".join(outs))
    return r.returncode == 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("action", choices=["prepare", "run", "compare", "all", "figures"])
    ap.add_argument("chains", nargs="*", help="R1..R7, T1 (default: all)")
    ap.add_argument("--work", default="_replay_work")
    a = ap.parse_args()
    work = Path(a.work).resolve()
    names = [SHORT.get(c, c) for c in a.chains] or list(CHAINS)
    for n in names:
        if n not in CHAINS:
            raise SystemExit(f"unknown chain {n}; choose from {', '.join(SHORT)}")
    if a.action in ("prepare", "all"):
        prepare(work)
    if a.action in ("run", "all"):
        run(work, names)
    if a.action == "figures":
        if not (work / "_workspace_info.json").exists():
            prepare(work)
        sys.exit(0 if figures(work) else 1)
    if a.action in ("compare", "all"):
        sys.exit(0 if compare(work, names if a.chains else None) else 1)


if __name__ == "__main__":
    main()
