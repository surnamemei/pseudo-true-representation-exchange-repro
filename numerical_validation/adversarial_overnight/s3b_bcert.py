"""Stage 3B driver: b-parameterised continuum certificate attempt (primary mpmath.iv + Arb replay).

Candidates (pre-registered order): [9.75,10.25], [9.90,10.10], [9.95,10.05], [9.98,10.02], [9.99,10.01].
Blocks of width 0.01 are certified centre-out; a candidate is certified iff all its blocks pass
(equivalent to descending-order testing with stop at the widest success). Total 3B cap: 150 min.
"""
from __future__ import annotations

import gzip
import json
import time
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
from decimal import Decimal
from pathlib import Path

from advrun import HERE, log, jdump, jload, status_update, MAX_WORKERS

SD = HERE / "stage3_bwidth" / "bcert"
CANDIDATES = [("9.75", "10.25"), ("9.90", "10.10"), ("9.95", "10.05"), ("9.98", "10.02"), ("9.99", "10.01")]
CAP_MIN = 150.0
PRIMARY_SHARE = 0.72


def blocks_centre_out():
    out = []
    for k in range(25):
        up = Decimal("10.00") + Decimal("0.01") * k
        dn = Decimal("10.00") - Decimal("0.01") * (k + 1)
        out += [str(up), str(dn)]
    return out


def _primary(b0):
    import bcert_block as bb
    t0 = time.time()
    try:
        out, rep = bb.certify_block(b0)
    except Exception as exc:  # pragma: no cover
        import traceback
        return b0, dict(ok=False, reason="exception", error=repr(exc), trace=traceback.format_exc(),
                        seconds=time.time() - t0), None
    return b0, out, rep


def _replay(b0):
    import bcert_arb as ba
    t0 = time.time()
    rep = json.loads(gzip.open(SD / "blocks" / f"{b0}_replay.json.gz", "rt").read())
    try:
        r = ba.replay_block(rep)
    except Exception as exc:
        import traceback
        r = dict(ok=False, error=repr(exc), trace=traceback.format_exc())
    r["seconds"] = time.time() - t0
    return b0, r


def run_pool(fn, items, deadline, desc):
    done = {}
    ex = ProcessPoolExecutor(max_workers=MAX_WORKERS)
    futs = {ex.submit(fn, it): it for it in items}
    pending = set(futs)
    while pending and time.time() < deadline:
        fin, pending = wait(pending, timeout=max(1.0, min(60.0, deadline - time.time())), return_when=FIRST_COMPLETED)
        for f in fin:
            res = f.result()
            done[res[0]] = res
            log(f"{desc}: block {res[0]} ok={res[1].get('ok')} ({res[1].get('seconds', 0):.0f}s) "
                f"[{len(done)}/{len(items)}]", 3)
            yield res
    for f in pending:
        f.cancel()
    ex.shutdown(wait=False, cancel_futures=True)


def block_set(lo, hi):
    a, b = Decimal(lo), Decimal(hi)
    out = []
    x = a
    while x < b:
        out.append(str(x))
        x += Decimal("0.01")
    return out


def run():
    SD.mkdir(parents=True, exist_ok=True)
    (SD / "blocks").mkdir(exist_ok=True)
    fin = HERE / "stage3_bwidth" / "DONE_3B.json"
    if fin.exists():
        return jload(fin)
    t_start = time.time()
    state_f = SD / "state.json"
    state = jload(state_f) if state_f.exists() else dict(started=time.time())
    t_start = state["started"]
    jdump(state_f, state)
    deadline_primary = t_start + CAP_MIN * 60 * PRIMARY_SHARE
    deadline_all = t_start + CAP_MIN * 60
    status_update(3, status="running_3B")
    order = blocks_centre_out()
    todo = [b for b in order if not (SD / "blocks" / f"{b}.json").exists()]
    for b0, out, rep in run_pool(_primary, todo, deadline_primary, "3B primary"):
        jdump(SD / "blocks" / f"{b0}.json", out)
        if rep is not None:
            with gzip.open(SD / "blocks" / f"{b0}_replay.json.gz", "wt") as f:
                f.write(json.dumps(rep))
    res = {b: jload(SD / "blocks" / f"{b}.json") for b in order if (SD / "blocks" / f"{b}.json").exists()}
    cand = []
    for lo, hi in CANDIDATES:
        bs = block_set(lo, hi)
        status = "certified_primary" if all(b in res and res[b].get("ok") for b in bs) else (
            "failed" if any(b in res and not res[b].get("ok") for b in bs) else "not_completed_budget")
        cand.append(dict(interval=[lo, hi], blocks=len(bs), status_primary=status))
    widest = next((c for c in cand if c["status_primary"] == "certified_primary"), None)
    # Arb replay of the blocks in the widest primary-certified candidate
    rep_res = {}
    if widest:
        bs = block_set(*widest["interval"])
        todo = [b for b in bs if not (SD / "blocks" / f"{b}_arb.json").exists()]
        for b0, r in run_pool(_replay, todo, deadline_all, "3B Arb replay"):
            jdump(SD / "blocks" / f"{b0}_arb.json", r)
        rep_res = {b: jload(SD / "blocks" / f"{b}_arb.json") for b in bs if (SD / "blocks" / f"{b}_arb.json").exists()}
        widest["arb_replay"] = ("all_pass" if len(rep_res) == len(bs) and all(r.get("ok") for r in rep_res.values())
                                else "incomplete_or_failed")
    final = None
    for c in cand:
        if c["status_primary"] == "certified_primary":
            if c is widest and widest.get("arb_replay") == "all_pass":
                final = c
            break
    blocks_ok = [res[b] for b in (block_set(*final["interval"]) if final else []) if b in res]
    summary = dict(candidates=cand, widest_primary=widest, final_certified=final,
                   verdict=("CERTIFIED_B_INTERVAL" if final else "NO_EXPLICIT_B_INTERVAL_CERTIFIED"),
                   blocks_attempted=len(res), blocks_ok=sum(1 for r in res.values() if r.get("ok")),
                   runtime_minutes=(time.time() - t_start) / 60, cap_minutes=CAP_MIN)
    if final:
        mm = [r["min_margins"] for r in blocks_ok]
        summary["certificate_record"] = dict(
            b_interval=final["interval"],
            lambda_range=[min(m["lambda_range"][0] for m in mm), max(m["lambda_range"][1] for m in mm)],
            min_Rvv_A=min(m["Rvv_A"] for m in mm), min_Rvv_B=min(m["Rvv_B"] for m in mm),
            min_beta_A=min(m["beta_A"] for m in mm), min_abs_beta_B=min(m["beta_B"] for m in mm),
            min_crossing_slope=min(m["slope"] for m in mm), min_gram=min(m["gram"] for m in mm),
            max_contraction=max(m["contraction_max"] for m in mm),
            min_coalescent_gap=min(m["coalescent"] for m in mm), min_tail_gap=min(m["tail"] for m in mm),
            terminal_cells=sum(r["leaves"] for r in blocks_ok), refinement_visits=sum(r["visits"] for r in blocks_ok),
            subtiles=sum(r["n_sub"] for r in blocks_ok),
            competitor_gap_note="other regular stationary points excluded cell-wise by cost/gradient/concavity/monotone predicates; minimum cost-cell margin is not separately tracked",
            replay="python-flint/Arb 256-bit independent replay: " + widest.get("arb_replay", "not run"),
            arb_replay_min_margins={b: rep_res[b].get("min_margins") for b in list(rep_res)[:3]},
            runtime_minutes=summary["runtime_minutes"])
    jdump(fin, summary)
    status_update(3, status="3B_done", verdict_3B=summary["verdict"])
    log(f"3B verdict {summary['verdict']}; candidates {[(c['interval'], c['status_primary']) for c in cand]}", 3)
    return summary
