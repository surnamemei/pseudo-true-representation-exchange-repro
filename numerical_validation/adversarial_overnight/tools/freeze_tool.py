"""Freeze tool: hash the pre-existing project state before any new experiment.

Writes freeze_manifest.json next to this tool's parent directory. Read-only with
respect to the project; the only file written is the manifest.
"""
from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE.parent
ROOT = OUT.parent.parent  # F:/Research

# Directories that belong to this project and whose files are hashed recursively.
RECURSIVE = ["paper", "supplement", "stage12_generalN", "stage12_replay", "manuscript"]
# Unrelated projects or bulky historical snapshots that are listed but not hashed file-by-file.
EXCLUDED_DIRS = ["GaussSt_TopoID", "RydbergComms", "graph_lasso_study", "matlab_prefs", "matlab_tmp",
                 "container_tools", "__pycache__", "validation", "TSP_SUBMISSION_FREEZE_v1",
                 "stage10_replay", "stage11_before_revision", "stage12_before_revision",
                 "stage12_deps", "stage12_pdf_preview", "stage13_before_revision",
                 "stage14_before_revision", "stage15_before_revision", "stage19_preview",
                 "stage19_replay_smoke", "stage23_before_revision", "stage4_bottleneck_cells",
                 "three_to_two_tone_figures", "three_to_two_tone_stage1_figures",
                 "three_to_two_tone_stage2_figures", "three_to_two_tone_stage3_figures",
                 "three_to_two_tone_stage4_figures", "three_to_two_tone_stage5_figures", "tmp"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def entry(path: Path) -> dict:
    st = path.stat()
    return {"path": path.relative_to(ROOT).as_posix(), "bytes": st.st_size,
            "mtime": time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime(st.st_mtime)),
            "sha256": sha256(path)}


def versions() -> dict:
    out = {"python": sys.version, "platform": platform.platform(), "machine": platform.machine(),
           "processor": platform.processor(), "cpu_count": os.cpu_count()}
    for mod in ["numpy", "scipy", "mpmath", "matplotlib", "sympy", "numba", "tqdm"]:
        try:
            m = __import__(mod)
            out[mod] = getattr(m, "__version__", "?")
        except Exception as exc:  # pragma: no cover
            out[mod] = f"unavailable: {exc}"
    try:
        sys.path.insert(0, str(ROOT / "stage12_deps"))
        import flint  # type: ignore
        out["python_flint"] = flint.__version__
        out["python_flint_location"] = str(Path(flint.__file__).parent)
    except Exception as exc:
        out["python_flint"] = f"unavailable: {exc}"
    try:
        import numpy as np
        cfg = np.show_config(mode="dicts")
        out["numpy_blas"] = cfg.get("Build Dependencies", {}).get("blas", {}).get("name")
        out["numpy_lapack"] = cfg.get("Build Dependencies", {}).get("lapack", {}).get("name")
    except Exception as exc:
        out["numpy_blas"] = f"unknown: {exc}"
    try:
        smi = subprocess.run(["nvidia-smi", "--query-gpu=name,driver_version,memory.used,memory.total,utilization.gpu",
                              "--format=csv,noheader"], capture_output=True, text=True, timeout=20)
        out["gpu"] = smi.stdout.strip()
    except Exception as exc:
        out["gpu"] = f"unavailable: {exc}"
    return out


def main() -> None:
    t0 = time.time()
    files = []
    for p in sorted(ROOT.iterdir()):
        if p.is_file():
            files.append(entry(p))
    for d in RECURSIVE:
        base = ROOT / d
        if base.exists():
            for p in sorted(base.rglob("*")):
                if p.is_file() and "__pycache__" not in p.parts:
                    files.append(entry(p))
    dirs = sorted(p.name for p in ROOT.iterdir() if p.is_dir())
    manifest = {
        "freeze_time_local": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "project_root": str(ROOT),
        "git": "NO_GIT_REPOSITORY (F:/Research is not a git work tree; SHA-256 manifest replaces HEAD)",
        "hashed_file_count": len(files),
        "hashed_total_bytes": sum(f["bytes"] for f in files),
        "recursive_dirs_hashed": RECURSIVE,
        "top_level_dirs_present": dirs,
        "dirs_not_hashed_file_by_file": [d for d in dirs if d in EXCLUDED_DIRS],
        "tsp_freeze_zip_sha256": sha256(ROOT / "TSP_SUBMISSION_FREEZE_v1.zip") if (ROOT / "TSP_SUBMISSION_FREEZE_v1.zip").exists() else None,
        "software": versions(),
        "files": files,
        "hash_seconds": round(time.time() - t0, 2),
    }
    (OUT / "freeze_manifest.json").write_text(json.dumps(manifest, indent=1), encoding="utf-8")
    print(f"hashed {len(files)} files, {manifest['hashed_total_bytes']/1e9:.3f} GB in {manifest['hash_seconds']} s")


if __name__ == "__main__":
    main()
