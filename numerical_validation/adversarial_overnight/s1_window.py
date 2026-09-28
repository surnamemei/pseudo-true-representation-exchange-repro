"""Stage 1: window / metric invariance."""
from __future__ import annotations

import csv
import math

import numpy as np
from scipy.optimize import brentq

from advrun import HERE, log, jdump, jload, status_update, run_configs
from advcore import Objective, make_signal, weights, refine_from, pair_distance, global_search

SD = HERE / "stage1_window"
PRIMARY = ["rect", "hann", "hamming", "dpss"]
SECONDARY = ["hann2"]
THETAS = [round(0.1 * k, 1) for k in range(11)]
VERDICT = {"present": "PERSISTS", "absent": "DOES_NOT_PERSIST", "inconclusive": "INCONCLUSIVE"}


def cfg(window, tag=None):
    return dict(cid=f"S1_{tag or window}".replace(":", "_").replace(".", "p"), stage=1, N=21, z=2.0, b=10.0,
                phi=math.pi, delta=0.0, amp_minus=1.0, amp_plus=1.0, window=window, scale="abs", lo=0.0, hi=1.0)


def first_present(a):
    ex = [r for r in a["switches"] if r.get("type") == "exchange"]
    pres = [r for r in ex if r["status"] == "present"]
    pick = pres[0] if pres else (ex[0] if ex else None)
    return pick, ex


def homotopy_continuation(rect_cross):
    """Local continuation of the rectangular A/B crossing along w_theta (theta step 0.01)."""
    N = 21
    A = tuple(rect_cross["P"][k] for k in ("u1", "u2"))
    B = tuple(rect_cross["Q"][k] for k in ("u1", "u2"))
    ec = rect_cross["eps_c"]
    rows = []
    for j in range(0, 101):
        th = j / 100
        w = weights(f"homotopy:{th}", N)
        cache = {}

        def branches(e):
            o = Objective(make_signal(N, 2.0, 10.0, e), w)
            a, b = refine_from(o, [A, B], tol_g=1e-12)
            cache[e] = (a, b)
            return a, b

        def gap(e):
            a, b = branches(e)
            return a.J - b.J
        lo, hi = ec - 0.02, ec + 0.02
        ok = True
        for _ in range(30):
            glo, ghi = gap(lo), gap(hi)
            if glo < 0 < ghi:
                break
            lo, hi = lo - 0.02, hi + 0.02
            if lo < 0:
                lo = 1e-6
        else:
            ok = False
        if not ok:
            rows.append(dict(theta=th, status="bracket_failed"))
            break
        e = brentq(gap, lo, hi, xtol=1e-14, rtol=1e-15)
        a, b = branches(e)
        jumpA = float(pair_distance(a.pair(), A, N))
        jumpB = float(pair_distance(b.pair(), B, N))
        rows.append(dict(theta=th, eps_c=e, A_u1=a.u1, A_u2=a.u2, B_u1=b.u1, B_u2=b.u2, J=0.5 * (a.J + b.J),
                         A_eig_min=a.eig_u[0], B_eig_min=b.eig_u[0], A_sep=a.sep, B_sep=b.sep,
                         jump_A=jumpA, jump_B=jumpB, status="ok" if max(jumpA, jumpB) < 0.5 else "jump"))
        A, B, ec = a.pair(), b.pair(), e
    return rows


def run(ctx=None):
    SD.mkdir(parents=True, exist_ok=True)
    done = SD / "DONE.json"
    if done.exists():
        log("Stage 1 already complete", 1)
        return jload(done)
    status_update(1, status="running")
    cfgs = [cfg(w) for w in PRIMARY + SECONDARY] + [cfg(f"homotopy:{t}", f"homotopy_{t}") for t in THETAS]
    res = run_configs(1, SD, cfgs, desc="S1 windows")
    rows = []
    for c in cfgs:
        a = res[c["cid"]]
        pick, ex = first_present(a)
        r = dict(config=c["cid"], window=c["window"], role=("primary" if c["window"] in PRIMARY else "secondary"),
                 verdict=VERDICT[a["verdict"]], extended=a["extended"], n_switches=len(a["switches"]),
                 n_exchanges=len(ex))
        if pick:
            r.update(eps_c=pick["eps_c"], J_star=pick["J_star"], A_u1=pick["P"]["u1"], A_u2=pick["P"]["u2"],
                     B_u1=pick["Q"]["u1"], B_u2=pick["Q"]["u2"], A_sep=pick["P"]["sep"], B_sep=pick["Q"]["sep"],
                     A_eig_u_min=pick["P"]["eig_u_min"], B_eig_u_min=pick["Q"]["eig_u_min"],
                     A_amp_max=pick["P"]["amp_max"], B_amp_max=pick["Q"]["amp_max"], slope=pick["slope"],
                     third_margin=pick["third_margin"], third_margin_rel=pick["third_margin_rel"],
                     coalescent_margin=pick["coalescent_margin"], branch_distance=pick["branch_distance"],
                     status=pick["status"], C1=pick["C1"], C2=pick["C2"], C3=pick["C3"], C4=pick["C4"], C5=pick["C5"],
                     approaches_coalescence=pick["approaches_coalescence"],
                     third_u1=(pick["third"] or {}).get("u1"), third_u2=(pick["third"] or {}).get("u2"))
        rows.append(r)
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with (SD / "window_summary.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    # minima arrays
    arrays = {}
    for c in cfgs:
        a = res[c["cid"]]
        eps = np.array([cv["eps"] for cv in a["curves"]])
        J = np.full((len(eps), 3), np.nan)
        U = np.full((len(eps), 3, 2), np.nan)
        cls = np.zeros((len(eps), 3), int)
        for i, cv in enumerate(a["curves"]):
            for j, q in enumerate(cv["minima"][:3]):
                J[i, j] = q["J"]
                U[i, j] = (q["u1"], q["u2"])
                cls[i, j] = 1 if q["cls"] == "regular" else 0
        key = c["cid"]
        arrays[f"{key}__eps"] = eps
        arrays[f"{key}__J"] = J
        arrays[f"{key}__U"] = U
        arrays[f"{key}__regular"] = cls
    np.savez_compressed(SD / "window_minima.npz", **arrays)
    # homotopy continuation
    rect = next(r for r in res[cfg("rect")["cid"]]["switches"] if r.get("type") == "exchange")
    hom = homotopy_continuation(rect)
    jdump(SD / "homotopy_continuation.json", hom)
    # matching of continued crossing with CA at theta grid
    match = []
    for t in THETAS:
        a = res[cfg(f"homotopy:{t}", f"homotopy_{t}")["cid"]]
        pick, _ = first_present(a)
        hc = next((h for h in hom if abs(h.get("theta", -1) - t) < 1e-9 and h.get("status") in ("ok", "jump")), None)
        if pick and hc:
            dA = float(pair_distance((pick["P"]["u1"], pick["P"]["u2"]), (hc["A_u1"], hc["A_u2"]), 21))
            dB = float(pair_distance((pick["Q"]["u1"], pick["Q"]["u2"]), (hc["B_u1"], hc["B_u2"]), 21))
            match.append(dict(theta=t, ca_eps_c=pick["eps_c"], cont_eps_c=hc["eps_c"], dA=dA, dB=dB,
                              ca_status=pick["status"], same_branches=bool(dA < 1e-4 and dB < 1e-4)))
        else:
            match.append(dict(theta=t, ca_eps_c=pick["eps_c"] if pick else None,
                              cont_eps_c=hc["eps_c"] if hc else None, same_branches=False))
    jdump(SD / "homotopy_match.json", match)
    plots(res, cfgs, hom)
    summary = dict(rows=rows, homotopy_match=match,
                   homotopy_complete=bool(hom and hom[-1].get("theta") == 1.0 and all(h.get("status") == "ok" for h in hom)),
                   primary_verdicts={r["window"]: r["verdict"] for r in rows if r["role"] == "primary"})
    summary["n_nonrect_persist"] = sum(v == "PERSISTS" for k, v in summary["primary_verdicts"].items() if k != "rect")
    jdump(done, summary)
    status_update(1, status="done", verdicts=summary["primary_verdicts"])
    log(f"Stage 1 verdicts {summary['primary_verdicts']}", 1)
    return summary


def plots(res, cfgs, hom):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    show = [c for c in cfgs if c["window"] in PRIMARY + SECONDARY]
    fig, axes = plt.subplots(2, len(show), figsize=(3.1 * len(show), 5.6), sharex=True)
    for j, c in enumerate(show):
        a = res[c["cid"]]
        eps = [cv["eps"] for cv in a["curves"]]
        for r in range(3):
            axes[0, j].plot(eps, [cv["minima"][r]["J"] if len(cv["minima"]) > r else np.nan for cv in a["curves"]],
                            lw=1, label=f"rank {r + 1}")
            axes[1, j].plot(eps, [cv["minima"][0]["u1"] for cv in a["curves"]], "b.", ms=2)
            axes[1, j].plot(eps, [cv["minima"][0]["u2"] for cv in a["curves"]], "r.", ms=2)
        for s in a["switches"]:
            if s.get("type") == "exchange":
                for ax in axes[:, j]:
                    ax.axvline(s["eps_c"], color="k" if s["status"] == "present" else "orange", ls="--", lw=0.8)
        axes[0, j].set_title(f"{c['window']} ({VERDICT[a['verdict']]})", fontsize=9)
        axes[1, j].set_xlabel(r"$\epsilon$")
    axes[0, 0].set_ylabel("J of three lowest minima")
    axes[1, 0].set_ylabel("best-fit frequencies u")
    axes[0, 0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(SD / "window_branches.png", dpi=130)
    fig.savefig(SD / "window_branches.pdf")
    plt.close(fig)
    ok = [h for h in hom if "eps_c" in h]
    if ok:
        fig, ax = plt.subplots(1, 2, figsize=(8, 3))
        ax[0].plot([h["theta"] for h in ok], [h["eps_c"] for h in ok], "-")
        ax[0].set(xlabel=r"homotopy $\theta$ (0=rect, 1=Hann)", ylabel=r"continued $\epsilon_c$")
        for k, lab in (("A_u1", "A u1"), ("A_u2", "A u2"), ("B_u1", "B u1"), ("B_u2", "B u2")):
            ax[1].plot([h["theta"] for h in ok], [h[k] for h in ok], label=lab)
        ax[1].legend(fontsize=7)
        ax[1].set(xlabel=r"$\theta$", ylabel="u at crossing")
        fig.tight_layout()
        fig.savefig(SD / "window_homotopy.png", dpi=130)
        plt.close(fig)
