"""Stage 2: record-origin / strong-pair phase invariance of the z^2 law."""
from __future__ import annotations

import csv
import math

import numpy as np

from advrun import HERE, log, jdump, jload, status_update, run_configs, pmap
from advtangent import TangentProblem, analyse_crossing

SD = HERE / "stage2_phase"
KAPPAS = [0.0, 0.05, -0.05, 0.10, -0.10, 0.20, -0.20, 0.30, -0.30, 0.40, -0.40, 0.50, -0.50]
KDIAG = [0.3882278594, -0.3882278594]
ZS_A = [0.25, 0.5, 1.0, 2.0]
DELTAS = [0.05, 0.10, 0.20, 0.30, 0.40]  # delta / pi
ZS_B = [0.05, 0.1, 0.2, 0.4, 0.8]
ZFIT = [0.05, 0.1, 0.2, 0.4]
QUAD_COEF = 1.92838866512018
FROZEN_KMAP = HERE.parent.parent / "kappa_crossing_map.csv"
FIXED_V = (1.61825508737687492, 10.3047708177652624)
B = 10.0


def cfgA(kappa, z):
    return dict(cid=f"S2A_k{kappa:+.10f}_z{z}".replace(".", "p").replace("+", "P").replace("-", "M"), stage=2, N=21,
                z=z, b=B, phi=math.pi, delta=kappa * z, amp_minus=1.0, amp_plus=1.0, window="rect", scale="z2",
                lo=0.0, hi=1.0, kappa=kappa)


def cfgB(dpi, z):
    return dict(cid=f"S2B_d{dpi}_z{z}".replace(".", "p"), stage=2, N=21, z=z, b=B, phi=math.pi,
                delta=dpi * math.pi, amp_minus=1.0, amp_plus=1.0, window="rect", scale="z", lo=0.0, hi=1.0,
                delta_over_pi=dpi)


# ----------------------------------------------------------------------------- tangent level
def _tangent_task(args):
    kappa, phi, seeds = args
    tp = TangentProblem(21, B, phi, kappa)
    out = None
    errs = []
    for sd in seeds:
        try:
            r = analyse_crossing(tp, *sd)
            if abs(r["R_diff"]) < 1e-10 and r["Rvv_A"] > 0 and r["Rvv_B"] > 0 and abs(r["vA"]) > 1e-3:
                r["seed"] = list(sd)
                out = r
                break
            errs.append(f"seed {sd}: invalid (vA={r['vA']:.4g}, Rvv=({r['Rvv_A']:.3g},{r['Rvv_B']:.3g}))")
        except Exception as exc:
            errs.append(f"seed {sd}: {exc!r}")
    return dict(kappa=kappa, phi=phi, result=out, errors=errs)


def tangent_level():
    frozen = {}
    if FROZEN_KMAP.exists():
        for r in csv.DictReader(FROZEN_KMAP.open()):
            frozen[float(r["kappa"])] = (float(r["v_A"]), float(r["v_B"]), float(r["lambda_c"]))
    base = (-4.153440486648769, 12.468094844300916, 0.06598041112699136)
    tasks = []
    for k in KAPPAS + KDIAG:
        seeds = [base]
        if round(k, 2) in frozen:
            seeds.insert(0, frozen[round(k, 2)])
        for kk in sorted(frozen, key=lambda q: abs(abs(q) - abs(k))):
            if frozen[kk] not in seeds:
                seeds.append(frozen[kk])
        tasks.append((k, math.pi, seeds[:6]))
        # secondary: true record shift rotates the weak tone by b*kappa
        tasks.append((k, math.pi + B * k, seeds[:6]))
    res = pmap(_tangent_task, tasks, desc="S2 tangent", stage=2)
    return res


# ----------------------------------------------------------------------------- analysis helpers
def pick_crossing(a, vA, vB):
    ex = [r for r in a["switches"] if r.get("type") == "exchange"]
    if not ex:
        return None, None

    def sat(p):
        return p["u1"] if abs(p["u1"]) > abs(p["u2"]) else p["u2"]

    def dist(r):
        sp, sq = sat(r["P"]), sat(r["Q"])
        return min(max(abs(sp - vA), abs(sq - vB)), max(abs(sp - vB), abs(sq - vA)))
    pres = [r for r in ex if r["status"] == "present"]
    pool_ = pres if pres else ex
    best = min(pool_, key=dist)
    return best, dist(best)


def fits(z, e):
    z, e = np.asarray(z, float), np.asarray(e, float)
    out = {}
    A = np.vstack([np.log(z), np.ones_like(z)]).T
    coef, *_ = np.linalg.lstsq(A, np.log(e), rcond=None)
    out["alpha"] = float(coef[0])
    out["log_resid_rms"] = float(np.sqrt(np.mean((A @ coef - np.log(e)) ** 2)))
    for k in (1, 2):
        C = float(np.sum(z ** k / e) / np.sum(z ** (2 * k) / e ** 2))
        rel = (C * z ** k - e) / e
        out[f"C{k}"] = C
        out[f"relres_rms_z{k}"] = float(np.sqrt(np.mean(rel ** 2)))
        out[f"relres_max_z{k}"] = float(np.max(np.abs(rel)))
    return out


def run(ctx=None):
    SD.mkdir(parents=True, exist_ok=True)
    done = SD / "DONE.json"
    if done.exists():
        log("Stage 2 already complete", 2)
        return jload(done)
    status_update(2, status="running_tangent")
    tf = SD / "tangent_results.json"
    if not tf.exists():
        jdump(tf, tangent_level())
    tang = jload(tf)
    tan_main = {t["kappa"]: t["result"] for t in tang if abs(t["phi"] - math.pi) < 1e-12}
    tan_true = {t["kappa"]: t["result"] for t in tang if abs(t["phi"] - math.pi) >= 1e-12}
    status_update(2, status="running_finite")
    cA = [cfgA(k, z) for k in KAPPAS + KDIAG for z in ZS_A]
    cB = [cfgB(d, z) for d in DELTAS for z in ZS_B]
    res = run_configs(2, SD, cA + cB, desc="S2 finite")
    # ---------------- 2A table
    rowsA = []
    for c in cA:
        k, z = c["kappa"], c["z"]
        t = tan_main.get(k)
        a = res[c["cid"]]
        vA, vB = (t["vA"], t["vB"]) if t else (-4.15, 12.47)
        pick, d = pick_crossing(a, vA, vB)
        r = dict(kappa=k, z=z, delta=k * z, verdict=a["verdict"], tangent_lambda=(t or {}).get("lam"),
                 tangent_vA=(t or {}).get("vA"), tangent_vB=(t or {}).get("vB"))
        if pick:
            r.update(eps_c=pick["eps_c"], eps_c_over_z2=pick["eps_c"] / z ** 2, status=pick["status"],
                     sat_match_distance=d, P_u=(pick["P"]["u1"], pick["P"]["u2"]),
                     Q_u=(pick["Q"]["u1"], pick["Q"]["u2"]), third_margin=pick["third_margin"],
                     third_margin_rel=pick["third_margin_rel"], coalescent_margin=pick["coalescent_margin"],
                     min_eig=min(pick["P"]["eig_u_min"], pick["Q"]["eig_u_min"]),
                     approaches_coalescence=pick["approaches_coalescence"], slope=pick["slope"],
                     J_star=pick["J_star"])
            if t:
                r["rel_dev_from_tangent"] = (pick["eps_c"] / z ** 2 - t["lam"]) / t["lam"]
        rowsA.append(r)
    # ---------------- metrics
    def row(k, z):
        return next((r for r in rowsA if abs(r["kappa"] - k) < 1e-12 and abs(r["z"] - z) < 1e-12), None)
    metrics = {}
    for k in KAPPAS + KDIAG:
        rr = [row(k, z) for z in (0.25, 0.5, 1.0)]
        present = [x is not None and x["verdict"] == "present" and x.get("status") == "present" for x in rr]
        devs = [abs(x.get("rel_dev_from_tangent", np.nan)) if x else np.nan for x in rr]
        metrics[k] = dict(M2A1_present_all_small_z=all(present), present=present, abs_rel_dev=devs,
                          M2A2_dev_decreasing=bool(np.all(np.isfinite(devs)) and devs[0] < devs[1] < devs[2]),
                          M2A2_dev_below_5pct_z025=bool(np.isfinite(devs[0]) and devs[0] < 0.05))
    even = {}
    for k in [x for x in KAPPAS if x > 0] + [KDIAG[0]]:
        tp, tm = tan_main.get(k), tan_main.get(-k)
        e = dict(tangent_rel_diff=(abs(tp["lam"] - tm["lam"]) / tp["lam"]) if tp and tm else None)
        for z in ZS_A:
            a, b = row(k, z), row(-k, z)
            if a and b and a.get("eps_c") and b.get("eps_c"):
                e[f"finite_z{z}_rel_diff"] = abs(a["eps_c"] - b["eps_c"]) / a["eps_c"]
        even[k] = e
    lam0 = tan_main[0.0]["lam"]
    quad = {}
    for k in (0.05, 0.10, 0.20):
        t = tan_main.get(k)
        if t:
            obs, pred = t["lam"] - lam0, QUAD_COEF * k * k
            quad[k] = dict(observed_shift=obs, predicted_shift=pred, rel_err=(obs - pred) / obs)
    # ---------------- 2B
    rowsB = []
    fitsB = {}
    for dpi in DELTAS:
        pts = []
        for z in ZS_B:
            c = cfgB(dpi, z)
            a = res[c["cid"]]
            pick, d = pick_crossing(a, *FIXED_V)
            r = dict(delta_over_pi=dpi, z=z, verdict=a["verdict"], tangent_mu_pred=0.6274474 * math.sin(dpi * math.pi))
            if pick:
                r.update(eps_c=pick["eps_c"], eps_c_over_z=pick["eps_c"] / z, eps_c_over_z2=pick["eps_c"] / z ** 2,
                         status=pick["status"], sat_match_distance=d, P_u=(pick["P"]["u1"], pick["P"]["u2"]),
                         Q_u=(pick["Q"]["u1"], pick["Q"]["u2"]), third_margin=pick["third_margin"],
                         third_margin_rel=pick["third_margin_rel"], coalescent_margin=pick["coalescent_margin"],
                         approaches_coalescence=pick["approaches_coalescence"], J_star=pick["J_star"],
                         n_exchanges=sum(1 for s in a["switches"] if s.get("type") == "exchange"))
                if z in ZFIT and pick["status"] == "present":
                    pts.append((z, pick["eps_c"]))
            rowsB.append(r)
        f = dict(n_points=len(pts))
        if len(pts) == len(ZFIT):
            zz, ee = zip(*pts)
            f["full"] = fits(zz, ee)
            f["drop_smallest"] = fits(zz[1:], ee[1:])
            f["drop_largest"] = fits(zz[:-1], ee[:-1])
            f["alpha_in_0p8_1p2"] = bool(0.8 <= f["full"]["alpha"] <= 1.2)
            f["z_fit_better_than_z2"] = bool(f["full"]["relres_rms_z1"] < f["full"]["relres_rms_z2"])
        fitsB[dpi] = f
    # ---------------- classification
    small = [k for k in KAPPAS if abs(k) <= 0.20 + 1e-12]
    cond_A = all(metrics[k]["M2A1_present_all_small_z"] and metrics[k]["M2A2_dev_decreasing"]
                 and metrics[k]["M2A2_dev_below_5pct_z025"] for k in small)
    cond_even = all((1.0 if even[k].get("tangent_rel_diff") is None else even[k]["tangent_rel_diff"]) < 1e-9 and
                    all(v < 1e-6 for kk, v in even[k].items() if kk.startswith("finite") and v is not None)
                    for k in (0.05, 0.10, 0.20))
    cond_quad = quad.get(0.05, {}).get("rel_err") is not None and abs(quad[0.05]["rel_err"]) < 0.10
    cond_B = all(fitsB[d].get("alpha_in_0p8_1p2") and fitsB[d].get("z_fit_better_than_z2") for d in DELTAS)
    fragile = not (metrics[0.05]["M2A1_present_all_small_z"] and metrics[-0.05]["M2A1_present_all_small_z"]
                   and metrics[0.05]["M2A2_dev_below_5pct_z025"])
    if cond_A and cond_even and cond_quad and cond_B:
        cls = "QUADRATIC_LAW_LOCALLY_STABLE_BUT_ALIGNMENT_SPECIFIC"
    elif fragile:
        cls = "QUADRATIC_LAW_EXTREMELY_FRAGILE"
    else:
        cls = "INCONCLUSIVE"
    # ---------------- 2C table
    table = []
    for k in [x for x in KAPPAS if x >= 0]:
        r = dict(kappa=k, fraction_of_record=k, weak_phase_rotation_rad=B * k,
                 weak_phase_rotation_over_pi=B * k / math.pi)
        for n in (21, 31, 101, 1001):
            r[f"tau_samples_N{n}"] = k * n
        table.append(r)
    with (SD / "kappa_to_record_shift.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(table[0]))
        w.writeheader()
        w.writerows(table)
    true_shift = {k: tan_true.get(k) for k in KAPPAS + KDIAG}
    for name, rows in (("stage2A_finite.csv", rowsA), ("stage2B_fixed_delta.csv", rowsB)):
        keys = list(dict.fromkeys(k for r in rows for k in r))
        with (SD / name).open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            w.writerows(rows)
    tkeys = ["kappa", "lam", "vA", "vB", "R_star", "Rvv_A", "Rvv_B", "beta_A_abs", "beta_B_abs", "slope",
             "coalescent_margin", "third_margin", "n_regular_minima", "best_scan_is_AB"]
    with (SD / "stage2A_tangent.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(tkeys + ["weak_phase"])
        for t in tang:
            r = t["result"]
            w.writerow([t["kappa"]] + ([r.get(k) for k in tkeys[1:]] if r else ["FAILED"] * (len(tkeys) - 1))
                       + ["pi" if abs(t["phi"] - math.pi) < 1e-12 else f"pi+b*kappa"])
    summary = dict(classification=cls, cond_kappa_small=cond_A, cond_even=cond_even, cond_quad=cond_quad,
                   cond_fixed_delta=cond_B, fragile_trigger=fragile, metrics=metrics, evenness=even,
                   quadratic=quad, fixed_delta_fits=fitsB, tangent=tan_main, true_record_shift_tangent=true_shift)
    plots(rowsA, rowsB, tan_main, fitsB)
    jdump(done, summary)
    status_update(2, status="done", classification=cls)
    log(f"Stage 2 classification {cls}", 2)
    return summary


def plots(rowsA, rowsB, tan, fitsB):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    ks = sorted(k for k in tan if tan[k])
    ax[0].plot(ks, [tan[k]["lam"] for k in ks], "k-o", ms=3, label=r"tangent $\lambda_{21}(\kappa)$")
    for z in ZS_A:
        rr = sorted([r for r in rowsA if r["z"] == z and r.get("eps_c")], key=lambda r: r["kappa"])
        ax[0].plot([r["kappa"] for r in rr], [r["eps_c_over_z2"] for r in rr], ".", ms=6, label=fr"finite $\epsilon_c/z^2$, z={z}")
    kk = np.linspace(-0.2, 0.2, 50)
    ax[0].plot(kk, tan[0.0]["lam"] + QUAD_COEF * kk ** 2, "g--", lw=0.8, label="frozen quadratic prediction")
    ax[0].set(xlabel=r"$\kappa=\delta/z$", ylabel=r"$\lambda$")
    ax[0].legend(fontsize=7)
    for d in DELTAS:
        rr = sorted([r for r in rowsB if r["delta_over_pi"] == d and r.get("eps_c")], key=lambda r: r["z"])
        ax[1].loglog([r["z"] for r in rr], [r["eps_c"] for r in rr], "o-", ms=3,
                     label=fr"$\delta/\pi$={d} ($\alpha$={fitsB[d].get('full', {}).get('alpha', float('nan')):.3f})")
    zz = np.array([0.05, 0.8])
    ax[1].loglog(zz, 0.3 * zz, "k:", lw=0.8, label=r"$\propto z$")
    ax[1].loglog(zz, 0.3 * zz ** 2, "k--", lw=0.8, label=r"$\propto z^2$")
    ax[1].set(xlabel="z", ylabel=r"$\epsilon_c$ (fixed $\delta$)")
    ax[1].legend(fontsize=6)
    fig.tight_layout()
    fig.savefig(SD / "phase_scaling.png", dpi=130)
    fig.savefig(SD / "phase_scaling.pdf")
    plt.close(fig)
