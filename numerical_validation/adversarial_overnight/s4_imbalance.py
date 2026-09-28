"""Stage 4: strong-pair amplitude imbalance A- = 1-eta (tone at -z), A+ = 1+eta (tone at +z)."""
from __future__ import annotations

import csv
import math

from advrun import HERE, log, jdump, jload, status_update, run_configs

SD = HERE / "stage4_imbalance"
ETAS = [0.0, 0.02, -0.02, 0.05, -0.05, 0.10, -0.10, 0.20, -0.20]


def cfg(eta):
    return dict(cid=f"S4_eta{eta:+.2f}".replace(".", "p").replace("+", "P").replace("-", "M"), stage=4, N=21,
                z=2.0, b=10.0, phi=math.pi, delta=0.0, amp_minus=1.0 - eta, amp_plus=1.0 + eta, window="rect",
                scale="abs", lo=0.0, hi=1.0, eta=eta)


def run(ctx=None):
    SD.mkdir(parents=True, exist_ok=True)
    done = SD / "DONE.json"
    if done.exists():
        log("Stage 4 already complete", 4)
        return jload(done)
    status_update(4, status="running")
    cfgs = [cfg(e) for e in ETAS]
    res = run_configs(4, SD, cfgs, desc="S4 imbalance")
    rows = []
    for c in cfgs:
        a = res[c["cid"]]
        ex = [s for s in a["switches"] if s.get("type") == "exchange"]
        pres = [s for s in ex if s["status"] == "present"]
        p = pres[0] if pres else (ex[0] if ex else None)
        r = dict(eta=c["eta"], verdict=a["verdict"], extended=a["extended"], n_exchanges=len(ex))
        if p:
            r.update(eps_c=p["eps_c"], J_star=p["J_star"], A_u1=p["P"]["u1"], A_u2=p["P"]["u2"], B_u1=p["Q"]["u1"],
                     B_u2=p["Q"]["u2"], A_eig=p["P"]["eig_u_min"], B_eig=p["Q"]["eig_u_min"], slope=p["slope"],
                     third_margin=p["third_margin"], third_margin_rel=p["third_margin_rel"],
                     coalescent_margin=p["coalescent_margin"], status=p["status"], status_orig=p.get("status_orig"),
                     approaches_coalescence=p["approaches_coalescence"],
                     min_sep=min(p["P"]["sep"], p["Q"]["sep"]))
        rows.append(r)
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with (SD / "imbalance_summary.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    v = {r["eta"]: r["verdict"] for r in rows}
    if any(x == "inconclusive" for x in v.values()):
        decisive_inconclusive = True
    else:
        decisive_inconclusive = False
    if all(x == "present" for x in v.values()):
        cls = "PERSISTS_OVER_TESTED_IMBALANCE"
    elif any(v[e] == "absent" for e in (0.02, -0.02)):
        cls = "DOES_NOT_PERSIST"
    elif all(v[e] == "present" for e in (0.0, 0.02, -0.02, 0.05, -0.05)) and any(
            v[e] == "absent" for e in (0.10, -0.10, 0.20, -0.20)):
        cls = "NARROW_IMBALANCE_ROBUSTNESS"
    else:
        cls = "INCONCLUSIVE"
    out = dict(classification=cls, rows=rows, any_inconclusive=decisive_inconclusive,
               evidence_level="GLOBAL_NUMERICAL")
    plots(rows)
    jdump(done, out)
    status_update(4, status="done", classification=cls)
    log(f"Stage 4 classification {cls}", 4)
    return out


def plots(rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rr = sorted([r for r in rows if "eps_c" in r], key=lambda r: r["eta"])
    fig, ax = plt.subplots(1, 2, figsize=(8, 3))
    ax[0].plot([r["eta"] for r in rr], [r["eps_c"] for r in rr], "o-")
    ax[0].set(xlabel=r"imbalance $\eta$ ($A_\pm=1\pm\eta$)", ylabel=r"$\epsilon_c$")
    ax[1].plot([r["eta"] for r in rr], [r["third_margin_rel"] for r in rr], "o-")
    ax[1].set(xlabel=r"$\eta$", ylabel="third margin / J*")
    fig.tight_layout()
    fig.savefig(SD / "imbalance.png", dpi=130)
    plt.close(fig)
