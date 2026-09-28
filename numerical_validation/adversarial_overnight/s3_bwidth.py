"""Stage 3A: parameter-width scan in b (continuum tangent + finite N21 z2).  3B lives in s3b_bcert.py."""
from __future__ import annotations

import csv
import math
import sys

import numpy as np

from advrun import HERE, PROJECT, log, jdump, jload, status_update, run_configs, pmap
from advtangent import TangentProblem, analyse_crossing, continuum_tail_margin

SD = HERE / "stage3_bwidth"
BGRID = [round(9.5 + 0.1 * k, 1) for k in range(11)]
ANCHOR = (-4.060030216180106892789, 12.42390445039620145, 0.0665069703923955354559)


def cfgF(b):
    return dict(cid=f"S3F_b{b}".replace(".", "p"), stage=3, N=21, z=2.0, b=b, phi=math.pi, delta=0.0,
                amp_minus=1.0, amp_plus=1.0, window="rect", scale="abs", lo=0.0, hi=1.0)


def _cont_task(args):
    b, seed = args
    tp = TangentProblem(21, b, math.pi, 0.0, continuum=True)
    r = analyse_crossing(tp, *seed, scan_pts=20001)
    tail, base = continuum_tail_margin(r["lam"], b, r["R_star"])
    r["tail_margin"] = tail
    r["b"] = b
    return r


def passes(r):
    return bool(r["Rvv_A"] > 0 and r["Rvv_B"] > 0 and r["beta_A_abs"] >= 1e-3 and r["beta_B_abs"] >= 1e-3
                and r["slope"] >= 1e-3 and r["third_margin"] > 0 and r["coalescent_margin"] > 0
                and r["tail_margin"] > 0 and abs(r["vA"]) >= 0.1 and abs(r["vB"]) >= 0.1
                and abs(r["R_diff"]) < 1e-12)


def continuum_scan():
    # sequential continuation seeds (cheap float64 solves), then parallel full analyses
    seeds = {10.0: ANCHOR}
    tp = None
    for direction in (+1, -1):
        prev = ANCHOR
        for k in range(1, 6):
            b = round(10.0 + direction * 0.1 * k, 1)
            tp = TangentProblem(21, b, math.pi, 0.0, continuum=True)
            x = tp.crossing(*prev)
            seeds[b] = tuple(x)
            prev = tuple(x)
    res = pmap(_cont_task, [(b, seeds[b]) for b in BGRID], desc="S3A continuum", stage=3)
    return res


def crosscheck_b10(r10):
    """Frozen closed form (stage22_b_map.loss) vs the independent GL implementation."""
    sys.path.insert(0, str(PROJECT))
    import stage22_b_map as s22
    from mpmath import mp
    out = {}
    for v in (r10["vA"], r10["vB"], 1.0, 7.3, -20.0):
        Rs = float(s22.loss(mp.mpf(v), mp.mpf(r10["lam"]), mp.mpf(10)))
        Rg = float(TangentProblem(21, 10.0, math.pi, 0.0, continuum=True).R([v], r10["lam"])[0][0])
        out[str(v)] = dict(closed_form=Rs, gauss_legendre=Rg, abs_diff=abs(Rs - Rg))
    out["max_abs_diff"] = max(d["abs_diff"] for d in out.values() if isinstance(d, dict))
    out["lambda_diff_vs_frozen"] = abs(r10["lam"] - float("0.0665069703923955354559"))
    out["pass"] = bool(out["max_abs_diff"] < 1e-10 and out["lambda_diff_vs_frozen"] < 1e-10)
    return out


def run(ctx=None):
    SD.mkdir(parents=True, exist_ok=True)
    done = SD / "DONE_3A.json"
    if done.exists():
        log("Stage 3A already complete", 3)
        out = jload(done)
    else:
        status_update(3, status="running_3A")
        rows = continuum_scan()
        rows.sort(key=lambda r: r["b"])
        for r in rows:
            r["pass"] = passes(r)
        r10 = next(r for r in rows if r["b"] == 10.0)
        xc = crosscheck_b10(r10)
        radius = 0.0
        for rr in (0.1, 0.2, 0.3, 0.4, 0.5):
            inside = [r for r in rows if abs(r["b"] - 10.0) <= rr + 1e-9]
            if all(r["pass"] for r in inside):
                radius = rr
            else:
                break
        keys = ["b", "lam", "vA", "vB", "R_star", "Rvv_A", "Rvv_B", "beta_A_re", "beta_B_re", "slope",
                "coalescent_margin", "third_margin", "third_v", "tail_margin", "n_regular_minima", "best_scan_is_AB",
                "pass"]
        with (SD / "b_continuum_scan.csv").open("w", newline="") as f:
            w = csv.writer(f)
            w.writerow(keys)
            for r in rows:
                w.writerow([r[k] for k in keys])
        # finite N21 z2 secondary
        fin = run_configs(3, SD, [cfgF(b) for b in BGRID], desc="S3 finite b")
        frows = []
        for b in BGRID:
            a = fin[cfgF(b)["cid"]]
            ex = [s for s in a["switches"] if s.get("type") == "exchange"]
            pres = [s for s in ex if s["status"] == "present"]
            p = pres[0] if pres else (ex[0] if ex else None)
            fr = dict(b=b, verdict=a["verdict"], n_exchanges=len(ex))
            if p:
                fr.update(eps_c=p["eps_c"], J_star=p["J_star"], third_margin=p["third_margin"],
                          third_margin_rel=p["third_margin_rel"], A=(p["P"]["u1"], p["P"]["u2"]),
                          B=(p["Q"]["u1"], p["Q"]["u2"]), min_eig=min(p["P"]["eig_u_min"], p["Q"]["eig_u_min"]),
                          status=p["status"], approaches_coalescence=p["approaches_coalescence"])
            frows.append(fr)
        out = dict(continuum=rows, crosscheck_b10=xc, diagnostic_symmetric_radius=radius,
                   diagnostic_interval=[10 - radius, 10 + radius], finite_N21_z2=frows,
                   evidence_level="GLOBAL_NUMERICAL (tangent/continuum problem; diagnostic, not certified)")
        plots(rows, frows)
        jdump(done, out)
        status_update(3, status="3A_done", diagnostic_radius=radius)
        log(f"Stage 3A diagnostic radius {radius}; crosscheck pass {xc['pass']}", 3)
    # 3B certificate attempt (separate module, may be absent)
    try:
        import s3b_bcert
        out3b = s3b_bcert.run()
        out["certificate_3B"] = out3b
    except ImportError:
        log("s3b_bcert not present; 3B skipped for now", 3)
    return out


def plots(rows, frows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(11, 3.2))
    b = [r["b"] for r in rows]
    ax[0].plot(b, [r["lam"] for r in rows], "o-")
    ax[0].set(xlabel="b", ylabel=r"$\lambda_\infty(b)$")
    for k, lab in (("coalescent_margin", "coalescent gap"), ("third_margin", "regular competitor gap"),
                   ("tail_margin", "tail margin")):
        ax[1].plot(b, [r[k] for r in rows], "o-", label=lab)
    ax[1].axhline(0, color="k", lw=0.6)
    ax[1].legend(fontsize=7)
    ax[1].set(xlabel="b", ylabel="margin (continuum)")
    fb = [r["b"] for r in frows if "eps_c" in r]
    ax[2].plot(fb, [r["eps_c"] for r in frows if "eps_c" in r], "s-")
    ax[2].set(xlabel="b", ylabel=r"finite N21, z=2 $\epsilon_c$")
    fig.tight_layout()
    fig.savefig(SD / "b_width_scan.png", dpi=130)
    plt.close(fig)
