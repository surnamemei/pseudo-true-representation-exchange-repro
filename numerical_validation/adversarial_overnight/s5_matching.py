"""Stage 5: matching convention across N (fixed normalized geometry vs fixed physical frequencies)."""
from __future__ import annotations

import csv
import math

from advrun import HERE, log, jdump, jload, status_update, run_configs

SD = HERE / "stage5_matching"
NS = [11, 21, 31, 41, 61]
DELTA_PHYS = 2.0 / 21
OMEGA3_PHYS = 10.0 / 21


def cfg(n, conv):
    if conv == "A":
        z, b = 2.0, 10.0
    else:
        z, b = n * DELTA_PHYS, n * OMEGA3_PHYS
    return dict(cid=f"S5_{conv}_N{n}", stage=5, N=n, z=z, b=b, phi=math.pi, delta=0.0, amp_minus=1.0,
                amp_plus=1.0, window="rect", scale="abs", lo=0.0, hi=1.0, convention=conv)


def run(ctx=None):
    SD.mkdir(parents=True, exist_ok=True)
    done = SD / "DONE.json"
    if done.exists():
        log("Stage 5 already complete", 5)
        return jload(done)
    status_update(5, status="running")
    cfgs = [cfg(n, c) for c in ("A", "B") for n in NS]
    res = run_configs(5, SD, cfgs, desc="S5 matching")
    rows = []
    for c in cfgs:
        a = res[c["cid"]]
        ex = [s for s in a["switches"] if s.get("type") == "exchange"]
        pres = [s for s in ex if s["status"] == "present"]
        p = pres[0] if pres else (ex[0] if ex else None)
        r = dict(convention=c["convention"], N=c["N"], z=c["z"], b=c["b"], Delta_phys=c["z"] / c["N"],
                 omega3_phys=c["b"] / c["N"], pair_separation_bins=c["z"] / math.pi, verdict=a["verdict"],
                 extended=a["extended"], n_exchanges=len(ex))
        if p:
            r.update(eps_c=p["eps_c"], J_star=p["J_star"], third_margin=p["third_margin"],
                     third_margin_rel=p["third_margin_rel"], A=(p["P"]["u1"], p["P"]["u2"]),
                     B=(p["Q"]["u1"], p["Q"]["u2"]), status=p["status"],
                     approaches_coalescence=p["approaches_coalescence"])
        rows.append(r)
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with (SD / "matching_table.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    plots(rows)
    out = dict(rows=rows, evidence_level="GLOBAL_NUMERICAL")
    jdump(done, out)
    status_update(5, status="done")
    log("Stage 5 done", 5)
    return out


def plots(rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(5.5, 3.4))
    for conv, mk, lab in (("A", "o-", "A: fixed z=NΔ=2, b=Nω₃=10"),
                          ("B", "s--", "B: fixed Δ=2/21, ω₃=10/21 rad/sample")):
        rr = [r for r in rows if r["convention"] == conv and "eps_c" in r]
        ax.plot([r["N"] for r in rr], [r["eps_c"] for r in rr], mk, label=lab)
    ax.set(xlabel="record length N", ylabel=r"numerical $\epsilon_c$")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(SD / "matching.png", dpi=130)
    fig.savefig(SD / "matching.pdf")
    plt.close(fig)
