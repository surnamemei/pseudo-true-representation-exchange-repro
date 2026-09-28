"""One-shot validation of the sealed transition theory against the frozen Stage-6 Monte Carlo.

Written and sealed (SEAL_2) before any Monte Carlo outcome was read. The real run:
  python3 validate.py
  - verifies the sealed predictions (SEAL_2.sha256) and the frozen outcome files (SHA-256 in transition_spec.json);
  - refuses to run if results/transition_theory/validation/validation.json already exists (one shot);
  - applies the criteria of 00_THEORY_PLAN.md §8 exactly; no quantity is re-fitted.
Test mode, on synthetic stand-in files only:
  python3 validate.py --mock SUMMARY.csv DONE.json OUTDIR
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
TT = ROOT / "paper/transition_theory"
DET = ROOT / "results/transition_theory/deterministic"
SPEC = json.loads((TT / "transition_spec.json").read_text())
ETAS = SPEC["setting"]["eta_grid"]
SNRS = SPEC["setting"]["snr_db"]
PRIMARY = SPEC["setting"]["primary_snr_db"]
S6 = "validation/adversarial_overnight/stage6_noise/"


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def fnum(v):
    return float(v) if v not in (None, "") else float("nan")


def load_inputs(summary_path, done_path):
    summ = list(csv.DictReader(open(summary_path, newline="")))
    done = json.loads(Path(done_path).read_text())
    pred = list(csv.DictReader(open(DET / "predictions.csv", newline="")))
    det = json.loads((DET / "deterministic.json").read_text())
    return summ, done, pred, det


def pick(rows, snr, eta):
    hit = [r for r in rows if int(round(float(r["snr"]))) == snr and abs(float(r["eta"]) - eta) < 1e-9]
    if len(hit) != 1:
        raise SystemExit(f"expected one row for snr={snr}, eta={eta}; found {len(hit)}")
    return hit[0]


def probit_entry(done, snr):
    pw = done["probit_widths"]
    return pw[str(snr)] if str(snr) in pw else pw[snr]


def evaluate(summ, done, pred, det):
    table, per_snr = [], {}
    for snr in SNRS:
        rows = []
        for eta in ETAS:
            s, p = pick(summ, snr, eta), pick(pred, snr, eta)
            n = int(round(float(s["n"])))
            PA_emp = float(s["P_A"])
            nA = int(round(PA_emp * n))
            PA, PAl = float(p["P_A"]), float(p["P_A_local"])
            mse_emp, mse_pred = float(s["MSE_to_winner"]), float(p["MSE_pred"])
            r = dict(snr=snr, eta=eta, n=n, n_A=nA, P_A_emp=PA_emp, P_A_emp_lo=fnum(s.get("P_A_lo")),
                     P_A_emp_hi=fnum(s.get("P_A_hi")), P_B_emp=float(s["P_B"]), P_OUT_emp=float(s["P_OUT"]),
                     P_A_pred=PA, P_A_pred_local=PAl, err=PA - PA_emp, err_local=PAl - PA_emp,
                     MSE_emp=mse_emp, MSE_pred=mse_pred, within_pred=float(p["within"]), between_pred=float(p["between"]),
                     ratio_emp_over_pred=mse_emp / mse_pred, log10_ratio=math.log10(mse_emp / mse_pred),
                     MSE_pred_local=float(p["MSE_pred_local"]), ratio_emp_over_pred_local=mse_emp / float(p["MSE_pred_local"]),
                     reference=p["reference"], trC_A_pred=float(p["trC_A"]), trC_B_pred=float(p["trC_B"]),
                     trC_A_frozen=fnum(s.get("local_trace_A")), trC_B_frozen=fnum(s.get("local_trace_B")))
            if abs(eta) < 1e-12:
                r["MSE_emp_Aref"] = float(s["MSE_to_A"])
                r["MSE_pred_Aref"] = float(p["MSE_pred_Aref"])
                r["ratio_Aref"] = r["MSE_emp_Aref"] / r["MSE_pred_Aref"]
            rows.append(r)
        table += rows
        tot = sum(r["n"] for r in rows)
        brier = sum(r["n_A"] * (1 - r["P_A_pred"]) ** 2 + (r["n"] - r["n_A"]) * r["P_A_pred"] ** 2 for r in rows) / tot
        brier_l = sum(r["n_A"] * (1 - r["P_A_pred_local"]) ** 2 + (r["n"] - r["n_A"]) * r["P_A_pred_local"] ** 2
                      for r in rows) / tot
        brier_floor = sum(r["n_A"] * (1 - r["P_A_emp"]) ** 2 + (r["n"] - r["n_A"]) * r["P_A_emp"] ** 2 for r in rows) / tot
        pw = probit_entry(done, snr)
        dp = det["per_snr"][str(snr)]
        W_pred, W_exact, W_emp = dp["W_closed_form"], dp["W_exact_curve"], float(pw["width_10_90"])
        eta50_emp, eta50_pred = float(pw["eta50"]), det["eta50_pred_exact"]
        central = [r for r in rows if abs(r["eta"]) <= 0.1 + 1e-12]
        ok_c5 = sum(0.5 <= r["ratio_emp_over_pred"] <= 2.0 for r in central)
        bad_f3 = sum(not (0.1 <= r["ratio_emp_over_pred"] <= 10.0) for r in central)
        p_m20 = pick(summ, snr, -0.20)
        p_p20 = pick(summ, snr, 0.20)
        orient = float(p_m20["P_A"]) > 0.5 > float(p_p20["P_A"])
        slope = float(pw["slope"])
        mid_dev = abs(eta50_emp - eta50_pred)
        per_snr[snr] = dict(
            role="primary" if snr in PRIMARY else "stress test (not used in the classification)",
            max_abs_err=max(abs(r["err"]) for r in rows), mean_abs_err=sum(abs(r["err"]) for r in rows) / len(rows),
            max_abs_err_local=max(abs(r["err_local"]) for r in rows),
            mean_abs_err_local=sum(abs(r["err_local"]) for r in rows) / len(rows),
            brier_theory=brier, brier_theory_local=brier_l, brier_empirical_frequency_floor=brier_floor,
            W_pred_closed_form=W_pred, W_pred_exact_curve=W_exact, W_emp_probit=W_emp, W_ratio_pred_over_emp=W_pred / W_emp,
            W_ratio_exact_over_emp=W_exact / W_emp, eta50_emp=eta50_emp, eta50_pred=eta50_pred, eta50_abs_dev=mid_dev,
            eta50_dev_over_W_pred=mid_dev / W_pred, slope_emp=slope, intercept_emp=float(pw["intercept"]),
            P_OUT_max=max(r["P_OUT_emp"] for r in rows), P_OUT_mean=sum(r["P_OUT_emp"] for r in rows) / len(rows),
            MSE_ratio_min=min(r["ratio_emp_over_pred"] for r in rows), MSE_ratio_max=max(r["ratio_emp_over_pred"] for r in rows),
            MSE_central_within_0p5_2=ok_c5, MSE_central_outside_0p1_10=bad_f3,
            MSE_all_within_0p1_10=sum(0.1 <= r["ratio_emp_over_pred"] <= 10.0 for r in rows),
            local_trace_max_rel_diff=max(max(abs(r["trC_A_pred"] - r["trC_A_frozen"]) / r["trC_A_frozen"],
                                             abs(r["trC_B_pred"] - r["trC_B_frozen"]) / r["trC_B_frozen"]) for r in rows),
            C1=bool(orient and slope < 0 and mid_dev <= W_pred / 2),
            C1_parts=dict(orientation=bool(orient), slope_negative=bool(slope < 0), midpoint_within_half_W=bool(mid_dev <= W_pred / 2)),
            C2=bool(abs(W_pred / W_emp - 1) <= 0.25),
            C3=bool(sum(abs(r["err"]) for r in rows) / len(rows) <= 0.08),
            C4=bool(all(0.1 <= r["ratio_emp_over_pred"] <= 10.0 for r in rows)),
            C5=bool(ok_c5 >= 5),
            F2=bool((not orient) or slope >= 0 or mid_dev > W_pred),
            F3=bool(bad_f3 >= 5))
    w = {s: per_snr[s]["W_emp_probit"] for s in SNRS}
    wp = {s: det["per_snr"][str(s)]["W_closed_form"] for s in SNRS}
    scaling = dict(W30_over_W40_emp=w[30] / w[40], W30_over_W40_pred=wp[30] / wp[40],
                   W30_over_W40_pred_exact=det["ratios"]["W30_over_W40_exact_curve"],
                   W20_over_W30_emp=w[20] / w[30], W20_over_W30_pred=wp[20] / wp[30])
    scaling["F1"] = bool(abs(scaling["W30_over_W40_emp"] / scaling["W30_over_W40_pred"] - 1) > 0.25)
    strong = all(per_snr[s][c] for s in PRIMARY for c in ("C1", "C2", "C3", "C4", "C5"))
    fails = scaling["F1"] or any(per_snr[s]["F2"] or per_snr[s]["F3"] for s in PRIMARY)
    outcome = "FAILS_PREDICTIVELY" if fails else ("STRONG_PREDICTION" if strong else "USEFUL_FIRST_ORDER")
    local_vs_exact = {s: dict(mean_abs_err_exact=per_snr[s]["mean_abs_err"], mean_abs_err_local=per_snr[s]["mean_abs_err_local"])
                      for s in SNRS}
    return table, per_snr, scaling, outcome, local_vs_exact


def figures(table, per_snr, det, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.special import ndtr
    col = {20: "#166da2", 30: "#d77519", 40: "#32834b"}
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
    fine = [e / 1000 for e in range(-220, 221)]
    for snr in SNRS:
        rr = [r for r in table if r["snr"] == snr]
        e = [r["eta"] for r in rr]
        ax[0].errorbar(e, [r["P_A_emp"] for r in rr], yerr=[[max(0, r["P_A_emp"] - r["P_A_emp_lo"]) for r in rr],
                       [max(0, r["P_A_emp_hi"] - r["P_A_emp"]) for r in rr]], fmt="o", ms=3, color=col[snr],
                       label=f"Monte Carlo {snr} dB")
        ax[0].plot(e, [r["P_A_pred"] for r in rr], "-", color=col[snr], lw=1.2, label=f"theory (exact gap) {snr} dB")
        K = det["per_snr"][str(snr)]["K"]
        ax[0].plot(fine, [float(ndtr(-K * x)) for x in fine], ":", color=col[snr], lw=0.9)
        ax[1].semilogy(e, [r["MSE_emp"] for r in rr], "o", ms=3, color=col[snr], label=f"Monte Carlo {snr} dB")
        ax[1].semilogy(e, [r["MSE_pred"] for r in rr], "-", color=col[snr], lw=1.2, label=f"two-branch theory {snr} dB")
        ax[1].semilogy(e, [r["within_pred"] for r in rr], "--", color=col[snr], lw=0.7)
    ax[0].set(xlabel=r"$\eta$", ylabel="P(A)", title="Branch-A selection (dotted: local law $\\Phi(-K\\eta)$)")
    ax[1].set(xlabel=r"$\eta$", ylabel="global MSE (frequency pair, u units)", title="Global MSE (dashed: within-mode term)")
    ax[0].legend(fontsize=6)
    ax[1].legend(fontsize=6)
    fig.tight_layout()
    fig.savefig(out / "transition_theory_validation.png", dpi=140)
    fig.savefig(out / "transition_theory_validation.pdf")
    plt.close(fig)


def main(argv):
    mock = "--mock" in argv
    if mock:
        i = argv.index("--mock")
        summary_path, done_path, out = Path(argv[i + 1]), Path(argv[i + 2]), Path(argv[i + 3])
    else:
        out = ROOT / "results/transition_theory/validation"
        if (out / "validation.json").exists():
            raise SystemExit("validation already performed (one-shot); refusing to run again")
        seal = (TT / "SEAL_2.sha256").read_text().split("\n")
        bad = [ln for ln in seal if ln.strip() and sha(ROOT / ln.split(None, 1)[1].strip()) != ln.split()[0]]
        if bad:
            raise SystemExit(f"sealed files changed since SEAL_2: {bad}")
        summary_path, done_path = ROOT / S6 / "summary.csv", ROOT / S6 / "DONE.json"
        for rel, h in SPEC["validation_inputs_read_only_after_seal_2"].items():
            if sha(ROOT / rel) != h:
                raise SystemExit(f"frozen outcome file {rel} differs from its sealed fingerprint")
    out.mkdir(parents=True, exist_ok=True)
    summ, done, pred, det = load_inputs(summary_path, done_path)
    table, per_snr, scaling, outcome, local_vs_exact = evaluate(summ, done, pred, det)
    with (out / "comparison.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(dict.fromkeys(k for r in table for k in r)))
        w.writeheader()
        w.writerows(table)
    res = dict(mode="mock (synthetic stand-in data)" if mock else "one-shot validation on the frozen Stage-6 Monte Carlo",
               run_at=time.strftime("%Y-%m-%dT%H:%M:%S%z"), outcome=outcome, per_snr={str(k): v for k, v in per_snr.items()},
               width_scaling=scaling, local_vs_exact={str(k): v for k, v in local_vs_exact.items()},
               inputs=dict(summary=str(summary_path.relative_to(ROOT)) if not mock else str(summary_path),
                           done=str(done_path.relative_to(ROOT)) if not mock else str(done_path),
                           summary_sha256=sha(summary_path), done_sha256=sha(done_path)),
               criteria=SPEC["criteria"])
    (out / "validation.json").write_text(json.dumps(res, indent=1) + "\n")
    figures(table, per_snr, det, out)
    print(json.dumps(dict(outcome=outcome, width_scaling=scaling,
                          per_snr={s: {k: v for k, v in d.items() if k in ("C1", "C2", "C3", "C4", "C5", "F2", "F3",
                                                                          "W_pred_closed_form", "W_emp_probit", "W_ratio_pred_over_emp",
                                                                          "mean_abs_err", "max_abs_err", "brier_theory",
                                                                          "MSE_central_within_0p5_2", "MSE_ratio_min", "MSE_ratio_max")}
                                   for s, d in res["per_snr"].items()}), indent=1))


if __name__ == "__main__":
    main(sys.argv[1:])
