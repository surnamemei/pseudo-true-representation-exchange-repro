"""Stage 6: noisy estimator, solver invariance (6A), Monte Carlo (6B), local sandwich theory (6C)."""
from __future__ import annotations

import csv
import math
import time

import numpy as np

from advrun import HERE, log, jdump, jload, status_update, pmap
from advcore import Objective, make_signal, global_search, refine_from, pair_distance, rng_for, times

SD = HERE / "stage6_noise"
N, Z, BW = 21, 2.0, 10.0
EC = 0.2481906301722774
ETAS = [-0.20, -0.15, -0.10, -0.075, -0.05, -0.025, 0.0, 0.025, 0.05, 0.075, 0.10, 0.15, 0.20]
SNRS = [20, 30, 40]
SUB_ETAS = [-0.10, -0.025, 0.0, 0.025, 0.10]
BLOCK = 50
A_FROZEN = (-4.772447353918614, 0.605170320276862)
B_FROZEN = (-0.057653561926680, 12.46341935680166)


def noiseless_roots():
    """Continue A and B from the frozen crossing roots to every eta (float64)."""
    out = {}
    for eta in ETAS:
        e = EC * (1 + eta)
        a, b = A_FROZEN, B_FROZEN
        steps = 40
        for j in range(1, steps + 1):
            ee = EC + (e - EC) * j / steps
            o = Objective(make_signal(N, Z, BW, ee))
            qa, qb = refine_from(o, [a, b], tol_g=1e-12)
            a, b = qa.pair(), qb.pair()
        o = Objective(make_signal(N, Z, BW, e))
        qa, qb = refine_from(o, [a, b], tol_g=1e-12)
        out[eta] = dict(eps=e, A=qa.pair(), B=qb.pair(), J_A=qa.J, J_B=qb.J, A_eig=qa.eig_u, B_eig=qb.eig_u)
    return out


def label(pair, roots):
    dA = float(pair_distance(pair, roots["A"], N))
    dB = float(pair_distance(pair, roots["B"], N))
    lab = "A" if dA < 1.0 else ("B" if dB < 1.0 else "OUT")
    return lab, dA, dB


def block_task(args):
    snr, ei, blk, roots, level, with_ref = args
    eta = ETAS[ei]
    x = make_signal(N, Z, BW, EC * (1 + eta))
    sigma2 = float(np.vdot(x, x).real / N) * 10 ** (-snr / 10)
    rng = np.random.default_rng(np.random.SeedSequence([20260928, 6, snr, ei, blk]))
    rows = []
    for t in range(BLOCK):
        w = math.sqrt(sigma2 / 2) * (rng.normal(size=N) + 1j * rng.normal(size=N))
        o = Objective(x + w)
        mins = global_search(o, level, blind=True, rng_key=(6, snr, ei, blk, t))
        best = mins[0]
        lab, dA, dB = label(best.pair(), roots)
        r = dict(snr=snr, eta=eta, eta_index=ei, block=blk, trial=blk * BLOCK + t, J=best.J, u1=best.u1, u2=best.u2,
                 label=lab, dA=dA, dB=dB, second_J=(mins[1].J if len(mins) > 1 else None),
                 second_label=(label(mins[1].pair(), roots)[0] if len(mins) > 1 else None), n_minima=len(mins),
                 best_cls=best.cls, best_eig=best.eig_u[0])
        if with_ref:
            ref = global_search(o, "R", z=Z, b=BW, extra=[roots["A"], roots["B"]], rng_key=(66, snr, ei, blk, t))[0]
            rl = label(ref.pair(), roots)[0]
            leg = refine_from(o, [roots["A"], roots["B"]], tol_g=1e-11)
            lg = leg[0] if leg[0].J <= leg[1].J else leg[1]
            ll = label(lg.pair(), roots)[0]
            r.update(ref_J=ref.J, ref_u1=ref.u1, ref_u2=ref.u2, ref_label=rl, legacy_J=lg.J, legacy_label=ll,
                     prod_vs_ref_disagree=bool(abs(best.J - ref.J) > 1e-8 * max(1.0, abs(ref.J)) or lab != rl),
                     legacy_vs_ref_disagree=bool(abs(lg.J - ref.J) > 1e-8 * max(1.0, abs(ref.J)) or ll != rl))
        rows.append(r)
    return rows


# ----------------------------------------------------------------------------- local sandwich covariance
def local_cov(x, u, sigma2):
    s = times(N)
    v = np.exp(1j * np.outer(s, u))
    a = np.linalg.lstsq(v, x, rcond=None)[0]
    r = x - v @ a
    jac = np.column_stack([v[:, 0], 1j * v[:, 0], v[:, 1], 1j * v[:, 1], 1j * s * a[0] * v[:, 0], 1j * s * a[1] * v[:, 1]])
    second = np.zeros((N, 6, 6), complex)
    for k in range(2):
        second[:, 2 * k, 4 + k] = second[:, 4 + k, 2 * k] = 1j * s * v[:, k]
        second[:, 2 * k + 1, 4 + k] = second[:, 4 + k, 2 * k + 1] = -s * v[:, k]
        second[:, 4 + k, 4 + k] = -s * s * a[k] * v[:, k]
    jj = (jac.conj().T @ jac).real
    H = jj - np.einsum("n,nij->ij", np.conj(r), second).real
    Hi = np.linalg.inv(H)
    cov = sigma2 / 2 * Hi @ jj @ Hi
    # finite-difference validation of H against the loss gradient
    theta = np.r_[a[0].real, a[0].imag, a[1].real, a[1].imag, u]

    def grad(t):
        aa = np.array([t[0] + 1j * t[1], t[2] + 1j * t[3]])
        vv = np.exp(1j * np.outer(s, t[4:]))
        rr = x - vv @ aa
        j = np.column_stack([vv[:, 0], 1j * vv[:, 0], vv[:, 1], 1j * vv[:, 1], 1j * s * aa[0] * vv[:, 0],
                             1j * s * aa[1] * vv[:, 1]])
        return -2 * (j.conj().T @ rr).real
    st = 2e-6
    num = np.column_stack([(grad(theta + np.eye(6)[k] * st) - grad(theta - np.eye(6)[k] * st)) / (2 * st) for k in range(6)])
    fd = float(np.max(np.abs(num - 2 * H)) / np.max(np.abs(2 * H)))
    return cov[4:, 4:], float(np.linalg.eigvalsh(H)[0]), fd


def wilson(k, n):
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    zq = 1.959963984540054
    den = 1 + zq * zq / n
    c = (p + zq * zq / (2 * n)) / den
    h = zq * math.sqrt(p * (1 - p) / n + zq * zq / (4 * n * n)) / den
    return c - h, c + h


def run(ctx=None):
    SD.mkdir(parents=True, exist_ok=True)
    (SD / "blocks").mkdir(exist_ok=True)
    done = SD / "DONE.json"
    if done.exists():
        log("Stage 6 already complete", 6)
        return jload(done)
    status_update(6, status="running")
    rf = SD / "noiseless_roots.json"
    if not rf.exists():
        r = noiseless_roots()
        jdump(rf, {str(k): v for k, v in r.items()})
    roots = {float(k): v for k, v in jload(rf).items()}
    # ---------------- 6A solver invariance
    level_file = SD / "production_level.json"
    level = jload(level_file)["level"] if level_file.exists() else "P"
    sa = SD / "solver_invariance_6A.json"
    if not sa.exists():
        for attempt in range(2):
            tasks = [(snr, ETAS.index(e), 0, roots[e], level, True) for snr in SNRS for e in SUB_ETAS]
            t0 = time.time()
            res = pmap(block_task, tasks, desc=f"S6A level {level}", stage=6)
            rows = [x for blk in res for x in blk[:40]]
            n = len(rows)
            dis = sum(r["prod_vs_ref_disagree"] for r in rows)
            leg = sum(r["legacy_vs_ref_disagree"] for r in rows)
            rec = dict(level=level, records=n, prod_disagreements=dis, prod_rate=dis / n, legacy_disagreements=leg,
                       legacy_rate=leg / n, seconds=time.time() - t0,
                       ref_label_counts={k: sum(r["ref_label"] == k for r in rows) for k in ("A", "B", "OUT")},
                       details=[r for r in rows if r["prod_vs_ref_disagree"] or r["legacy_vs_ref_disagree"]][:60])
            jdump(SD / f"solver_invariance_6A_level{level}.json", rec)
            log(f"6A level {level}: prod disagreements {dis}/{n}, legacy {leg}/{n}", 6)
            if dis / n <= 0.01 or attempt == 1:
                rec["accepted_level"] = level
                jdump(sa, rec)
                jdump(level_file, dict(level=level))
                break
            level = "P2"
            with (HERE / "AMENDMENTS.md").open("a", encoding="utf-8") as f:
                f.write(f"\n## A-6A (automatic, {time.strftime('%Y-%m-%d %H:%M:%S%z')})\n6A production disagreement "
                        f"{dis}/{n} exceeded 1%; production upgraded once to level P2 (grid 2048x513, K_grid 128, "
                        f"n_rand 64) per the pre-registered rule, and 6A rerun.\n")
    sa_rec = jload(sa)
    level = sa_rec["accepted_level"]
    # ---------------- 6B main Monte Carlo
    n_trials = 1000
    nblocks = n_trials // BLOCK
    pending = [(snr, ei, b) for snr in SNRS for ei in range(len(ETAS)) for b in range(nblocks)
               if not (SD / "blocks" / f"s{snr}_e{ei}_b{b}.json").exists()]
    if pending:
        # pre-registered reduction rule: project after the first two blocks of every setting
        first = [p for p in pending if p[2] < 2]
        t0 = time.time()
        if first:
            res = pmap(block_task, [(s, e, b, roots[ETAS[e]], level, False) for s, e, b in first], desc="S6B first", stage=6)
            for (s, e, b), rr in zip(first, res):
                if isinstance(rr, list):
                    jdump(SD / "blocks" / f"s{s}_e{e}_b{b}.json", rr)
        el = time.time() - t0
        rest = [p for p in pending if p[2] >= 2]
        proj_min = (el / max(1, len(first))) * len(rest) / 60
        red = dict(first_blocks=len(first), seconds=el, remaining_blocks=len(rest), projected_minutes=proj_min,
                   reduced=bool(proj_min > 150))
        jdump(SD / "reduction_rule.json", red)
        if red["reduced"]:
            n_trials = 500
            nblocks = n_trials // BLOCK
            rest = [p for p in rest if p[2] < nblocks]
        log(f"6B projection {proj_min:.1f} min -> trials {n_trials}", 6)
        CH = 120
        for j in range(0, len(rest), CH):
            chunk = rest[j:j + CH]
            res = pmap(block_task, [(s, e, b, roots[ETAS[e]], level, False) for s, e, b in chunk],
                       desc=f"S6B {j // CH + 1}/{math.ceil(len(rest) / CH)}", stage=6)
            for (s, e, b), rr in zip(chunk, res):
                if isinstance(rr, list):
                    jdump(SD / "blocks" / f"s{s}_e{e}_b{b}.json", rr)
            status_update(6, blocks_done=len(list((SD / "blocks").glob("*.json"))))
    red = jload(SD / "reduction_rule.json") if (SD / "reduction_rule.json").exists() else dict(reduced=False)
    n_trials = 500 if red.get("reduced") else 1000
    out = analyse(roots, n_trials, sa_rec)
    jdump(done, out)
    status_update(6, status="done", hypotheses=out["hypotheses"])
    log(f"Stage 6 hypotheses {out['hypotheses']}", 6)
    return out


def analyse(roots, n_trials, sa_rec):
    nblocks = n_trials // BLOCK
    rows = []
    for snr in SNRS:
        for ei in range(len(ETAS)):
            for b in range(nblocks):
                rows += jload(SD / "blocks" / f"s{snr}_e{ei}_b{b}.json")
    keys = list(rows[0])
    with (SD / "trials.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    summ = []
    locs = {}
    fdmax = 0.0
    for snr in SNRS:
        for eta in ETAS:
            rr = [r for r in rows if r["snr"] == snr and abs(r["eta"] - eta) < 1e-12]
            x = make_signal(N, Z, BW, EC * (1 + eta))
            sigma2 = float(np.vdot(x, x).real / N) * 10 ** (-snr / 10)
            R = roots[eta]
            covA, hA, fdA = local_cov(x, np.array(R["A"]), sigma2)
            covB, hB, fdB = local_cov(x, np.array(R["B"]), sigma2)
            fdmax = max(fdmax, fdA, fdB)
            locs[(snr, eta)] = (covA, covB)
            U = np.array([[r["u1"], r["u2"]] for r in rr])
            L = np.array([r["label"] for r in rr])
            n = len(rr)
            s = dict(snr=snr, eta=eta, eps=EC * (1 + eta), n=n, sigma2=sigma2)
            for lab in ("A", "B", "OUT"):
                k = int(np.sum(L == lab))
                lo, hi = wilson(k, n)
                s[f"P_{lab}"] = k / n
                s[f"P_{lab}_lo"], s[f"P_{lab}_hi"] = lo, hi
            winner = np.array(R["A"] if eta < 0 else R["B"])
            s["MSE_to_winner"] = float(np.mean(np.sum((U - winner) ** 2, axis=1)))
            s["MSE_to_A"] = float(np.mean(np.sum((U - np.array(R["A"])) ** 2, axis=1)))
            s["MSE_to_B"] = float(np.mean(np.sum((U - np.array(R["B"])) ** 2, axis=1)))
            s["MSE_to_true_strong_pair"] = float(np.mean(np.sum((U - np.array([-Z, Z])) ** 2, axis=1)))
            for lab, cov, ref in (("A", covA, R["A"]), ("B", covB, R["B"])):
                q = U[L == lab]
                s[f"local_trace_{lab}"] = float(np.trace(cov))
                s[f"local_cov_{lab}"] = cov.tolist()
                if len(q) >= 2:
                    c = np.cov(q, rowvar=False)
                    s[f"cond_n_{lab}"] = len(q)
                    s[f"cond_mean_{lab}"] = q.mean(axis=0).tolist()
                    s[f"cond_cov_{lab}"] = c.tolist()
                    s[f"cond_trace_{lab}"] = float(np.trace(c))
                    s[f"cond_MSE_own_{lab}"] = float(np.mean(np.sum((q - np.array(ref)) ** 2, axis=1)))
                    s[f"ratio_trace_{lab}"] = float(np.trace(c) / np.trace(cov))
                else:
                    s[f"cond_n_{lab}"] = len(q)
            summ.append(s)
    with (SD / "summary.csv").open("w", newline="") as f:
        ks = list(dict.fromkeys(k for s in summ for k in s if not k.startswith(("local_cov", "cond_cov", "cond_mean"))))
        w = csv.DictWriter(f, fieldnames=ks, extrasaction="ignore")
        w.writeheader()
        w.writerows(summ)
    # hypotheses
    h1 = []
    for s in summ:
        if abs(s["eta"]) >= 0.10 - 1e-12 and s["snr"] in (30, 40):
            for lab in ("A", "B"):
                if s.get(f"cond_n_{lab}", 0) >= 100:
                    rt = s[f"ratio_trace_{lab}"]
                    h1.append(dict(snr=s["snr"], eta=s["eta"], branch=lab, n=s[f"cond_n_{lab}"], ratio=rt,
                                   ok=bool(0.8 <= rt <= 1.25)))
    h2 = []
    h3 = []
    for s in summ:
        if abs(s["eta"]) < 1e-12:
            ratio = s["MSE_to_A"] / s["local_trace_A"]
            h2.append(dict(snr=s["snr"], MSE_to_A=s["MSE_to_A"], local_trace_A=s["local_trace_A"], ratio=ratio,
                           ok=bool(ratio >= 10)))
            h3.append(dict(snr=s["snr"], P_A=s["P_A"], P_B=s["P_B"], ok=bool(min(s["P_A"], s["P_B"]) >= 0.2)))
    H1 = bool(h1) and all(x["ok"] for x in h1)
    H2 = all(x["ok"] for x in h2)
    H3 = all(x["ok"] for x in h3)
    # probit width per SNR (same MLE as archive; informative)
    widths = {}
    from scipy.special import ndtr, ndtri
    from scipy.optimize import minimize
    for snr in SNRS:
        ss = [s for s in summ if s["snr"] == snr]
        eta = np.array([s["eta"] for s in ss])
        k = np.array([s["P_A"] * s["n"] for s in ss])
        nn = np.array([s["n"] for s in ss], float)

        def nll(q):
            p = np.clip(ndtr(q[0] + q[1] * eta), 1e-12, 1 - 1e-12)
            return -np.sum(k * np.log(p) + (nn - k) * np.log1p(-p))
        f = minimize(nll, [0.0, -10.0], method="Nelder-Mead", options=dict(xatol=1e-10, fatol=1e-12, maxiter=4000))
        widths[snr] = dict(width_10_90=float((ndtri(0.9) - ndtri(0.1)) / abs(f.x[1])), intercept=float(f.x[0]),
                           slope=float(f.x[1]), eta50=float(-f.x[0] / f.x[1]))
    out = dict(n_trials_per_setting=n_trials, solver_invariance=dict((k, v) for k, v in sa_rec.items() if k != "details"),
               hypotheses=dict(H1=H1, H2=H2, H3=H3, TSP_consequence=bool(H1 and H2)), H1_detail=h1, H2_detail=h2,
               H3_detail=h3, probit_widths=widths, local_hessian_fd_max_rel_err=fdmax,
               evidence_level="GLOBAL_NUMERICAL (Monte Carlo with numerical global search)",
               local_theory_label="local Hessian (sandwich) covariance approximation - not a CRB/MCRB")
    plots(summ, rows, roots)
    return out


def plots(summ, rows, roots):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    col = {20: "#166da2", 30: "#d77519", 40: "#32834b"}
    fig, ax = plt.subplots(figsize=(5.5, 3.5))
    for snr in SNRS:
        ss = [s for s in summ if s["snr"] == snr]
        e = [s["eta"] for s in ss]
        ax.errorbar(e, [s["P_A"] for s in ss], yerr=[[max(0.0, s["P_A"] - s["P_A_lo"]) for s in ss], [max(0.0, s["P_A_hi"] - s["P_A"]) for s in ss]],
                    fmt="o-", ms=3, color=col[snr], label=f"P(A) {snr} dB")
        ax.plot(e, [s["P_OUT"] for s in ss], "x:", color=col[snr], label=f"P(out) {snr} dB")
    ax.set(xlabel=r"$\eta=(\epsilon-\epsilon_c)/\epsilon_c$", ylabel="selection probability")
    ax.legend(fontsize=6)
    fig.tight_layout()
    fig.savefig(SD / "selection_probability.png", dpi=130)
    fig.savefig(SD / "selection_probability.pdf")
    plt.close(fig)
    r0 = roots[0.0]
    d = np.array(r0["B"]) - np.array(r0["A"])
    d /= np.linalg.norm(d)
    fig, ax = plt.subplots(1, 3, figsize=(11, 3.2))
    for j, snr in enumerate(SNRS):
        rr = [r for r in rows if r["snr"] == snr and abs(r["eta"]) < 1e-12]
        U = np.array([[r["u1"], r["u2"]] for r in rr])
        c = (U - np.array(r0["A"])) @ d
        ax[j].hist(c, bins=60, color=col[snr])
        ax[j].set(title=f"{snr} dB, $\\eta=0$", xlabel="coordinate along A→B (normalized u)", ylabel="count")
    fig.tight_layout()
    fig.savefig(SD / "bimodality.png", dpi=130)
    plt.close(fig)
    fig, ax = plt.subplots(1, 3, figsize=(11, 3.2), sharey=True)
    for j, snr in enumerate(SNRS):
        ss = [s for s in summ if s["snr"] == snr]
        e = [s["eta"] for s in ss]
        ax[j].semilogy(e, [s["MSE_to_winner"] for s in ss], "k.-", label="global MSE to noiseless winner")
        ax[j].semilogy(e, [s.get("cond_MSE_own_A", np.nan) for s in ss], "b.-", label="A-conditioned MSE (own root)")
        ax[j].semilogy(e, [s.get("cond_MSE_own_B", np.nan) for s in ss], "r.-", label="B-conditioned MSE (own root)")
        ax[j].semilogy(e, [s["local_trace_A"] for s in ss], "b--", lw=0.8, label="local sandwich trace A")
        ax[j].semilogy(e, [s["local_trace_B"] for s in ss], "r--", lw=0.8, label="local sandwich trace B")
        ax[j].set(title=f"{snr} dB", xlabel=r"$\eta$")
    ax[0].set_ylabel("frequency-pair squared error")
    ax[0].legend(fontsize=6)
    fig.tight_layout()
    fig.savefig(SD / "global_vs_conditional_error.png", dpi=130)
    fig.savefig(SD / "global_vs_conditional_error.pdf")
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(5.5, 3.5))
    for snr in SNRS:
        ss = [s for s in summ if s["snr"] == snr]
        for lab, mk in (("A", "o"), ("B", "s")):
            pts = [(s["eta"], s[f"ratio_trace_{lab}"]) for s in ss if s.get(f"cond_n_{lab}", 0) >= 100]
            if pts:
                ax.plot(*zip(*pts), mk + "-", ms=3, color=col[snr], label=f"{lab}, {snr} dB")
    ax.axhspan(0.8, 1.25, color="0.9")
    ax.set_yscale("log")
    ax.set(xlabel=r"$\eta$", ylabel="MC conditional trace / local sandwich trace")
    ax.legend(fontsize=6, ncol=2)
    fig.tight_layout()
    fig.savefig(SD / "local_theory_comparison.png", dpi=130)
    fig.savefig(SD / "local_theory_comparison.pdf")
    plt.close(fig)
