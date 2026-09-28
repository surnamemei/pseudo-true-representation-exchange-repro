"""Write stage2_phase/STAGE2_REPORT.md from stage outputs (no new computation)."""
import csv
import json
from pathlib import Path

H = Path(__file__).resolve().parent.parent
d = json.loads((H / "stage2_phase/DONE.json").read_text())
ph = json.loads((H / "stage2_phase/POSTHOC_fixed_delta_all_switches.json").read_text())
rowsA = list(csv.DictReader((H / "stage2_phase/stage2A_finite.csv").open()))
t = d["tangent"]


def fr(k, z):
    return next((r for r in rowsA if abs(float(r["kappa"]) - k) < 1e-12 and abs(float(r["z"]) - z) < 1e-12), None)


def tan(k):
    for key, v in t.items():
        if abs(float(key) - k) < 1e-12:
            return v
    return None


L = ["# Stage 2 — record-origin / strong-pair phase invariance of the z² law", "",
     "Pre-registered classification (z² law only): **INCONCLUSIVE** — reason below. Existence of the exchange in general is assessed separately.", "",
     "## 2A. Spacing-scaled offset δ = κz (weak tone fixed at φ = π)", "",
     "Tangent level: independent realified least-squares implementation (not the frozen closed forms), N = 21. "
     "Finite level: full-torus crossing analysis at z ∈ {0.25, 0.5, 1, 2}, ε = λz², λ ∈ [0, 1] (+ one boundary extension when no qualifying crossing).", "",
     "| κ | tangent λ₂₁(κ) | v_A | v_B | coalescent margin | third margin | finite ε_c/z² at z = 0.25 / 0.5 / 1 / 2 | rel. dev. from λ₂₁(κ) at z = 0.25 / 0.5 / 1 | qualifying at z = 0.25 / 0.5 / 1 / 2 |",
     "|---|---|---|---|---|---|---|---|---|"]
for k in [0.0, 0.05, 0.10, 0.20, 0.30, 0.3882278594, 0.40, 0.50]:
    tt = tan(k)
    rr = [fr(k, z) for z in (0.25, 0.5, 1.0, 2.0)]
    ez = " / ".join(f"{float(r['eps_c_over_z2']):.6f}" if r and r.get("eps_c_over_z2") else "—" for r in rr)
    dv = " / ".join(f"{float(r['rel_dev_from_tangent']) * 100:+.3f}%" if r and r.get("rel_dev_from_tangent") else "—" for r in rr[:3])
    q = " / ".join(("yes" if r and r["verdict"] == "present" else "no") for r in rr)
    if tt:
        L.append(f"| {k:+.4g} | {tt['lam']:.8f} | {tt['vA']:.4f} | {tt['vB']:.4f} | {tt['coalescent_margin']:.3e} | {tt['third_margin']:.4f} | {ez} | {dv} | {q} |")
    else:
        L.append(f"| {k:+.4g} | no regular A/B root (A at the central chart) | — | — | — | — | {ez} | {dv} | {q} |")
q = d["quadratic"]
L += ["", "Negative κ gives identical values (exact conjugation/time-reversal symmetry): tangent relative differences ≤ 2e-16, finite-record relative differences ≤ 2.6e-11 (M2A.3 pass).",
      f"Quadratic shift (M2A.4): λ₂₁(0.05) − λ₂₁(0) = {q['0.05']['observed_shift']:.6e} vs frozen prediction 1.92838866512018·κ² = {q['0.05']['predicted_shift']:.6e} "
      f"(rel. error {q['0.05']['rel_err'] * 100:.2f}%; {q['0.1']['rel_err'] * 100:.2f}% at κ = 0.1; {q['0.2']['rel_err'] * 100:.1f}% at κ = 0.2 where O(κ⁴) matters).",
      "", "**2A (pre-registered M2A.1–M2A.4): pass for every |κ| ≤ 0.20.** A qualifying separated global exchange exists at every tested z, and ε_c/z² converges to the independently computed tangent λ₂₁(κ) (deviation 0.10–0.19% at z = 0.25, shrinking ∝ z²). "
      "For |κ| = 0.3 the small-z exchanges (z = 0.25, 0.5) fail the pre-registered regularity rule (A's within-fit separation 0.45–0.63 < 1.0 normalized units: A is near-coalescent) although they remain global switches; "
      "for |κ| ≥ 0.388 the A participant sits exactly on the confluent line (h = 0, unbounded two-tone amplitudes), matching the frozen analytical collision κ ≈ 0.3882 at N = 21.", "",
      "## 2B. Fixed nonzero offset δ (weak tone at φ = π)", "",
      "| δ/π | tangent ε/z | z = 0.05 | z = 0.1 | z = 0.2 | z = 0.4 | z = 0.8 (outside fit window) |", "|---|---|---|---|---|---|---|"]
for dpi in ("0.05", "0.1", "0.2", "0.3", "0.4"):
    o = ph[dpi]
    cells = []
    for r in o["rows"]:
        if "eps_over_z" in r:
            cells.append(f"ε/z = {r['eps_over_z']:.5f}; A sep {r['P_sep']:.2f} ({'regular' if r['P_cls'] == 'regular' else 'coalescent-class'})")
        else:
            cells.append("—")
    L.append(f"| {dpi} | {o['tangent_eps_over_z']:.5f} | " + " | ".join(cells) + " |")
L += ["", "Pre-registered 2B metric: exponent fit on the four in-window z using the *qualifying* crossing nearest the tangent pair. Qualifying crossings exist only at z ≤ 0.1 for δ/π ≤ 0.2, at z = 0.05 for δ/π = 0.3, and at no in-window z for δ/π = 0.4, because the lower-ε participant (A) has within-fit separation < 1.0 or lies on the confluent line. "
      "**The pre-registered exponent fit is therefore not evaluable → Stage-2 classification INCONCLUSIVE.**", "",
      "### POST_HOC_SENSITIVITY (does not replace the pre-registered result)", "",
      "Using every global switch into the B-like fit (satellite near 10.3), regardless of A's regularity, over z ∈ {0.05, 0.1, 0.2, 0.4}:", "",
      "| δ/π | α (log–log) | rel. RMS residual ε_c = C₁z | rel. RMS residual ε_c = C₂z² | α without smallest z | α without largest z |", "|---|---|---|---|---|---|"]
for dpi in ("0.05", "0.1", "0.2", "0.3", "0.4"):
    F = ph[dpi]["fits"]
    L.append(f"| {dpi} | {F['full']['alpha']:.4f} | {F['full']['relres_rms_z1']:.4f} | {F['full']['relres_rms_z2']:.4f} | {F['drop_smallest']['alpha']:.4f} | {F['drop_largest']['alpha']:.4f} |")
L += ["", "The global switch location scales as ε_c ∝ z (α = 0.986–0.994) and ε_c/z matches the tangent prediction 0.6274474|sin δ| to better than 1% at z = 0.05. "
      "The z² law does **not** hold at fixed δ; the fixed-δ exchange is a first-order (ε ∝ z) phenomenon whose lower participant becomes near-coalescent already at z ≈ 0.1–0.2.", "",
      "## 2C. Reinterpretation as window placement", "",
      "A shift of the record centre by τ samples gives κ = τ/N (fraction of record = κ). A *true* shift also rotates the weak tone by bκ rad.", "",
      "| κ | fraction of record | τ (N=21) | τ (N=31) | τ (N=101) | τ (N=1001) | weak-tone rotation bκ (rad) |", "|---|---|---|---|---|---|---|"]
for r in csv.DictReader((H / "stage2_phase/kappa_to_record_shift.csv").open()):
    L.append(f"| {float(r['kappa']):.2f} | {float(r['fraction_of_record']):.2f} | {float(r['tau_samples_N21']):.2f} | {float(r['tau_samples_N31']):.2f} | "
             f"{float(r['tau_samples_N101']):.1f} | {float(r['tau_samples_N1001']):.0f} | {float(r['weak_phase_rotation_rad']):.2f} |")
tt = d["true_record_shift_tangent"]
good = [(float(k), v) for k, v in tt.items() if v and float(k) > 0 and v["third_margin"] > 0 and v["coalescent_margin"] > 0 and abs(v["vA"] - v["vB"]) > 0.1]
L += ["", "Secondary (pre-registered) tangent check of a *true* record shift (δ = κz and φ = π + bκ): a valid global A/B tangent crossing exists for κ = "
      + ", ".join(f"{k:.2f} (λ₂₁ = {v['lam']:.5f})" for k, v in sorted(good, key=lambda kv: kv[0]))
      + "; for |κ| ≥ 0.4 the seeded solver returned no valid global A/B pair (degenerate or non-global roots) — recorded as FAILED for this secondary diagnostic.", "",
      "## Interpretation for the z² claim", "",
      "* The quadratic law is **locally stable in the manuscript's κ-scaled sense**: it survives record-centre displacements up to 0.2N samples (4.2 samples at N = 21) with λ(κ) = λ(0) + 1.93κ² + O(κ⁴), numerically at finite z to 0.1% agreement.",
      "* It is **alignment-specific**: it needs |δ| ≲ 0.2–0.3·z. For a fixed relative phase — the generic case when the beat period far exceeds the record (small z) — the crossing scales as ε ∝ z. The fraction of relative phases compatible with the z² regime shrinks ∝ z as z → 0.",
      "* The pre-registered label stays **INCONCLUSIVE** only because the fixed-δ fit required *regular* (separation ≥ 1.0) participants, which fixed-δ exchanges do not provide beyond z ≈ 0.1. The post-hoc evidence points to `QUADRATIC_LAW_LOCALLY_STABLE_BUT_ALIGNMENT_SPECIFIC`; it does not override the pre-registered label.", "",
      "Evidence levels: tangent rows GLOBAL_NUMERICAL (1-D full-period scan); finite rows GLOBAL_NUMERICAL; post-hoc fits POST_HOC_SENSITIVITY. Figures: `stage2_phase/phase_scaling.png|pdf` (pre-registered; its fixed-δ legend shows α = n/a because the pre-registered fit was not evaluable) and `stage2_phase/POSTHOC_fixed_delta_scaling.png` (POST_HOC_SENSITIVITY; filled markers = regular A participant, open = near-coalescent)."]
(H / "stage2_phase/STAGE2_REPORT.md").write_text("\n".join(L), encoding="utf-8")
print("ok", len(L))
