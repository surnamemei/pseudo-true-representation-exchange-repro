"""Runner utilities: status file, logging, resumable parallel configuration scans."""
from __future__ import annotations

import hashlib
import json
import math
import os
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent.parent
STATUS = HERE / "STATUS.json"
MAX_WORKERS = int(os.environ.get("ADV_WORKERS", "24"))


def now():
    return time.strftime("%Y-%m-%dT%H:%M:%S%z")


def log(msg, stage=None):
    line = f"[{time.strftime('%H:%M:%S')}]" + (f"[S{stage}]" if stage is not None else "") + " " + str(msg)
    print(line, flush=True)
    if stage is not None:
        d = HERE / "logs"
        d.mkdir(exist_ok=True)
        with (d / f"stage{stage}.log").open("a", encoding="utf-8") as f:
            f.write(line + "\n")


def jdump(path: Path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(obj, indent=1, default=_default), encoding="utf-8")
    os.replace(tmp, path)


def _default(o):
    import numpy as np
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (complex,)):
        return [o.real, o.imag]
    return str(o)


def jload(path: Path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def sha(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


# ----------------------------------------------------------------------------- status
def status_load():
    if STATUS.exists():
        try:
            return jload(STATUS)
        except Exception:
            pass
    return {"run_started": now(), "stages": {}, "events": []}


def status_update(stage=None, **kw):
    st = status_load()
    st["last_update"] = now()
    if stage is not None:
        s = st["stages"].setdefault(str(stage), {})
        s.update(kw)
    else:
        st.update(kw)
    jdump(STATUS, st)
    return st


def status_event(msg):
    st = status_load()
    st.setdefault("events", []).append(f"{now()} {msg}")
    st["last_update"] = now()
    jdump(STATUS, st)


def code_hashes():
    return {p.name: sha(p) for p in sorted(HERE.glob("*.py"))}


# ----------------------------------------------------------------------------- pool
_POOL = None


def pool():
    global _POOL
    if _POOL is None:
        _POOL = ProcessPoolExecutor(max_workers=MAX_WORKERS)
    return _POOL


def shutdown_pool():
    global _POOL
    if _POOL is not None:
        _POOL.shutdown(wait=True, cancel_futures=True)
        _POOL = None


def pmap(fn, tasks, desc="", stage=None, chunks=1):
    """Parallel map with progress; results returned in task order."""
    from tqdm import tqdm
    ex = pool()
    futs = {ex.submit(fn, t): i for i, t in enumerate(tasks)}
    out = [None] * len(tasks)
    with tqdm(total=len(tasks), desc=desc, file=sys.stdout, ascii=True, mininterval=10, ncols=100) as bar:
        for f in as_completed(futs):
            i = futs[f]
            try:
                out[i] = f.result()
            except Exception as exc:
                out[i] = {"error": repr(exc), "trace": traceback.format_exc()}
                log(f"task {i} failed: {exc!r}", stage)
            bar.update(1)
    return out


# ----------------------------------------------------------------------------- config scans
def run_configs(stage, sdir: Path, configs, validate_level="R", desc="scan"):
    """Resumable: scan (GS(P) at 101 points), analyse, apply boundary rule once."""
    from advcross import scan_task, scan_grid
    sdir = Path(sdir)
    (sdir / "scan").mkdir(parents=True, exist_ok=True)
    (sdir / "analysis").mkdir(parents=True, exist_ok=True)
    for i, c in enumerate(configs):
        c.setdefault("cid_int", int(hashlib.sha256(c["cid"].encode()).hexdigest()[:7], 16))
    # ---- base scans
    todo = [c for c in configs if not (sdir / "scan" / f"{c['cid']}.json").exists()]
    if todo:
        tasks = [(c, k, p, "P", False) for c in todo for k, p in enumerate(scan_grid(c))]
        log(f"{desc}: {len(todo)} configs x 101 points = {len(tasks)} GS(P) tasks", stage)
        # process in batches of configs so partial progress is checkpointed
        B = max(1, min(len(todo), 2000 // 101 + 1))
        for j in range(0, len(todo), B):
            batch = todo[j:j + B]
            bt = [(c, k, p, "P", False) for c in batch for k, p in enumerate(scan_grid(c))]
            res = pmap(_scan_wrap, bt, desc=f"{desc} {j // B + 1}/{math.ceil(len(todo) / B)}", stage=stage)
            for c in batch:
                pts = [r for (t, r) in zip(bt, res) if t[0]["cid"] == c["cid"]]
                if all(r and "error" not in r for r in pts):
                    jdump(sdir / "scan" / f"{c['cid']}.json", dict(config=c, points=pts))
                else:
                    log(f"scan incomplete for {c['cid']}; not checkpointed", stage)
            status_update(stage, scans_done=sum((sdir / "scan" / f"{c['cid']}.json").exists() for c in configs),
                          scans_total=len(configs))
    # ---- analysis
    todo = [c for c in configs if not (sdir / "analysis" / f"{c['cid']}.json").exists()
            and (sdir / "scan" / f"{c['cid']}.json").exists()]
    if todo:
        res = pmap(_analyse_wrap, [(c, str(sdir), validate_level) for c in todo], desc=f"{desc} analyse", stage=stage)
        for c, r in zip(todo, res):
            if r and "error" not in r:
                jdump(sdir / "analysis" / f"{c['cid']}.json", r)
            else:
                log(f"analysis failed for {c['cid']}: {r}", stage)
    # ---- boundary rule
    need_ext = []
    for c in configs:
        p = sdir / "analysis" / f"{c['cid']}.json"
        if p.exists():
            a = jload(p)
            if a["verdict"] == "absent" and not a.get("extended"):
                need_ext.append(c)
    if need_ext:
        log(f"boundary rule: extending {len(need_ext)} configs to [hi, 3hi]", stage)
        bt = [(c, k, p, "P", True) for c in need_ext for k, p in enumerate(scan_grid(c, extended=True))]
        res = pmap(_scan_wrap, bt, desc=f"{desc} extension", stage=stage)
        for c in need_ext:
            pts = [r for (t, r) in zip(bt, res) if t[0]["cid"] == c["cid"]]
            if all(r and "error" not in r for r in pts):
                jdump(sdir / "scan" / f"{c['cid']}__ext.json", dict(config=c, points=pts))
        need_ext = [c for c in need_ext if (sdir / "scan" / f"{c['cid']}__ext.json").exists()]
        res = pmap(_analyse_wrap, [(c, str(sdir), validate_level) for c in need_ext],
                   desc=f"{desc} analyse ext", stage=stage)
        for c, r in zip(need_ext, res):
            if r and "error" not in r:
                jdump(sdir / "analysis" / f"{c['cid']}.json", r)
    out = {}
    for c in configs:
        p = sdir / "analysis" / f"{c['cid']}.json"
        out[c["cid"]] = jload(p) if p.exists() else None
    return out


def _scan_wrap(args):
    from advcross import scan_task
    return scan_task(args)


def _analyse_wrap(args):
    from advcross import analyse_config, classify
    c, sdir, level = args
    sdir = Path(sdir)
    base = jload(sdir / "scan" / f"{c['cid']}.json")["points"]
    extended = (sdir / "scan" / f"{c['cid']}__ext.json").exists()
    pts = list(base)
    if extended:
        ext = jload(sdir / "scan" / f"{c['cid']}__ext.json")["points"]
        for d in ext:
            d = dict(d)
            d["k"] = d["k"] + len(base)
            pts.append(d)
    t0 = time.time()
    lists, results = analyse_config(c, pts, validate_level=level)
    verdict = classify(results)
    curves = [dict(k=pt["k"], eps=pt["eps"], minima=[q.as_dict() for q in L[:3]])
              for pt, L in zip(sorted(pts, key=lambda d: d["k"]), lists)]
    return dict(config=c, verdict=verdict, extended=extended, switches=results, curves=curves,
                seconds=time.time() - t0)
