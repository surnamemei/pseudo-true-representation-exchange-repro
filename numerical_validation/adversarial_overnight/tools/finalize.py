"""Finalize: propagate the Stage-3B outcome into the reports, regenerate tables, integrity check, final summary."""
import json
import subprocess
import sys
import time
from pathlib import Path

H = Path(__file__).resolve().parent.parent
py = sys.executable


def J(p):
    p = H / p
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else None


s3b = J("stage3_bwidth/DONE_3B.json")
if s3b is None:
    print("3B not finished")
    sys.exit(1)
cands = {tuple(c["interval"]): c for c in s3b["candidates"]}
done_blocks = s3b.get("blocks_attempted", 0)
ok_blocks = s3b.get("blocks_ok", 0)
if s3b["verdict"] == "CERTIFIED_B_INTERVAL":
    cr = s3b["certificate_record"]
    iv_ = cr["b_interval"]
    res3b = (f"**explicit certificate for b ∈ [{iv_[0]}, {iv_[1]}]** (continuum tangent problem, uniform in b; "
             f"mpmath.iv primary + independent Arb replay: all pass; λ∞ ∈ [{cr['lambda_range'][0]:.10f}, {cr['lambda_range'][1]:.10f}], "
             f"min R_vv,A {cr['min_Rvv_A']:.3e}, min slope {cr['min_crossing_slope']:.4f}, min coalescent gap {cr['min_coalescent_gap']:.3e}, "
             f"min tail gap {cr['min_tail_gap']:.3e}, {cr['terminal_cells']} terminal cells)")
    item3b = (f"Replace the existential b-neighbourhood statement by the certified continuum interval b ∈ [{iv_[0]}, {iv_[1]}] "
              f"(Stage 3B), keeping the existential statement for the finite-N families.")
    claim3b = f"explicit continuum b-interval [{iv_[0]}, {iv_[1]}] CERTIFIED (primary + Arb replay)"
else:
    reason = ("the pre-registered 150-min cap expired before any 0.01-wide block finished (mpmath.iv evaluation of the "
              "near-well flanks with b-ranges needs ~1e5–1e6 cell visits per block)") if done_blocks == 0 else (
              f"{ok_blocks}/{done_blocks} completed blocks passed; no candidate interval had all of its blocks certified within the cap")
    res3b = f"**NO_EXPLICIT_B_INTERVAL_CERTIFIED** — {reason}; the existential implicit-function statement is preserved"
    item3b = ("Keep the existential b-neighbourhood statement (Stage 3B produced no explicit certified interval). A future attempt "
              "should use Arb as the fast evaluator (≈40× faster here) or centred forms in b.")
    claim3b = "NO_EXPLICIT_B_INTERVAL_CERTIFIED (budget)" if done_blocks == 0 else "NO_EXPLICIT_B_INTERVAL_CERTIFIED"

mr = (H / "MORNING_REPORT.md").read_text(encoding="utf-8")
mr = mr.replace("Stage 3B explicit certificate: **[3B RESULT]**.", f"Stage 3B explicit certificate: {res3b}.")
mr = mr.replace("6. **[3B ITEM]**", f"6. {item3b}")
(H / "MORNING_REPORT.md").write_text(mr, encoding="utf-8")

ca = (H / "CLAIM_AUDIT.md").read_text(encoding="utf-8")
ca = ca.replace("| weak location b | [9.5, 10.5] continuum + finite z = 2 | persists; λ∞ varies < 1.6% | GLOBAL_NUMERICAL (diagnostic); 3B certificate: see Stage 3 |",
                f"| weak location b | [9.5, 10.5] continuum + finite z = 2 | persists; λ∞ varies < 1.6% | GLOBAL_NUMERICAL (diagnostic); 3B: {claim3b} |")
ca = ca.replace("4. b-robustness: numerical over [9.5, 10.5]; Stage 3B result: see `stage3_bwidth/STAGE3_REPORT.md` (explicit certified interval if listed there).",
                f"4. b-robustness: numerical over [9.5, 10.5]; Stage 3B: {claim3b}.")
ca = ca.replace("* Explicit parameter radii exist only where Stage 3B certifies them (b); κ and amplitude-ratio radii remain existential.",
                "* Explicit parameter radii: " + ("b only (Stage 3B certified interval); κ and amplitude-ratio radii remain existential."
                                                  if s3b["verdict"] == "CERTIFIED_B_INTERVAL" else
                                                  "none certified (Stage 3B did not produce one); b, κ and amplitude-ratio radii remain existential."))
(H / "CLAIM_AUDIT.md").write_text(ca, encoding="utf-8")

hr = (H / "HOSTILE_REVIEW.md").read_text(encoding="utf-8")
hr = hr.replace("explicit certificate: see Stage 3B in `stage3_bwidth/STAGE3_REPORT.md`.", f"explicit certificate: {claim3b}.")
(H / "HOSTILE_REVIEW.md").write_text(hr, encoding="utf-8")

for script in ("tools/report_stages346.py", "tools/report_stage2.py", "tools/make_reports.py", "tools/integrity_check.py"):
    subprocess.run([py, str(H / script)], cwd=H, check=True)

from_status = J("STATUS.json")
from_status["final"] = dict(completed=time.strftime("%Y-%m-%dT%H:%M:%S%z"), overall_status="NEEDS ONE FOCUSED VALIDATION PASS",
                            stage3B=s3b["verdict"])
(H / "STATUS.json").write_text(json.dumps(from_status, indent=1), encoding="utf-8")
print("finalized; 3B:", claim3b)
