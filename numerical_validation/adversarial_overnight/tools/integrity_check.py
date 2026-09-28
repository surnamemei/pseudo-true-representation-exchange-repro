"""Re-hash every frozen file listed in freeze_manifest.json and report any change."""
import hashlib
import json
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
ROOT = HERE.parent.parent
man = json.loads((HERE / "freeze_manifest.json").read_text(encoding="utf-8"))
changed, missing = [], []
for f in man["files"]:
    p = ROOT / f["path"]
    if not p.exists():
        missing.append(f["path"])
        continue
    h = hashlib.sha256()
    with p.open("rb") as fh:
        for b in iter(lambda: fh.read(8 << 20), b""):
            h.update(b)
    if h.hexdigest() != f["sha256"]:
        changed.append(f["path"])
new_top = sorted(p.name for p in ROOT.iterdir() if p.is_file())
known = {Path(f["path"]).name for f in man["files"] if "/" not in f["path"]}
added = [n for n in new_top if n not in known]
out = dict(checked=len(man["files"]), changed=changed, missing=missing, new_top_level_files=added,
           all_frozen_files_unchanged=not changed and not missing, time=time.strftime("%Y-%m-%dT%H:%M:%S%z"))
(HERE / "integrity_check.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
print(json.dumps(out, indent=1))
