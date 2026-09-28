"""Regenerate the revision figures from already-completed results (no new scientific computation).

Inputs: archived CSV/JSON results in the project root and in validation/adversarial_overnight/.
The only evaluation performed is the closed-form N=21 tangent loss for the explanatory Fig. 2,
which reproduces the archived figure (stage17_tangent_landscape.py) with the new branch notation.
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar

ROOT = Path(__file__).resolve().parents[2]
VAL = ROOT / "validation" / "adversarial_overnight"
FIG = ROOT / "paper" / "figures"
SUP = ROOT / "supplement" / "figures"
FIG.mkdir(exist_ok=True)
SUP.mkdir(exist_ok=True)

C1, C2, C3 = "#2a78d6", "#eb6834", "#1baf7a"      # validated categorical slots 1-3
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#d9d8d4"
plt.rcParams.update({"font.family": "serif", "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8,
                     "legend.fontsize": 6.8, "xtick.labelsize": 7, "ytick.labelsize": 7, "pdf.fonttype": 42,
                     "axes.edgecolor": MUTED, "axes.linewidth": 0.6, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.labelcolor": INK, "lines.linewidth": 1.4, "mathtext.fontset": "cm"})
BA, BB = r"$\mathcal{A}$", r"$\mathcal{B}$"
SINGLE, DOUBLE = 3.5, 7.16


def J(p):
    return json.loads(Path(p).read_text(encoding="utf-8"))


def style(ax):
    ax.grid(alpha=0.35, color=GRID, lw=0.5)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


def tag(ax, s):
    ax.text(-0.02, 1.02, s, transform=ax.transAxes, va="bottom", ha="right", fontsize=8, fontweight="bold", color=INK)


# ----------------------------------------------------------------------------- Fig. 1 continuation
def fig1():
    rows = list(csv.DictReader((ROOT / "branch_continuation.csv").open()))
    ec = float(J(ROOT / "stage2_reference.json")["epsilon_cross"])
    fig, ax = plt.subplots(figsize=(DOUBLE, 2.45))
    for br, col, lab in (("A", C1, BA), ("B", C2, BB)):
        rs = sorted([r for r in rows if r["branch"] == br], key=lambda r: float(r["epsilon"]))
        e = [float(r["epsilon"]) for r in rs]
        ax.plot(e, [float(r["u1"]) for r in rs], "-", color=col, lw=1.5, label=f"{lab}: lower fitted frequency")
        ax.plot(e, [float(r["u2"]) for r in rs], "--", color=col, lw=1.5, label=f"{lab}: upper fitted frequency")
    for val in (-2, 2, 10):
        ax.axhline(val, color="0.6", lw=0.7, ls=":")
        ax.text(0.302, val, f"{val:g}", va="center", color=MUTED, fontsize=7)
    ax.axvline(ec, color=INK, ls="-.", lw=0.9)
    ax.text(ec - 0.003, 14.6, r"certified $\epsilon_c$", ha="right", fontsize=7, color=INK)
    ax.set(xlabel=r"Omitted-tone amplitude $\epsilon$", ylabel=r"Fitted frequency $u=N\nu$", xlim=(0, 0.30))
    style(ax)
    ax.legend(ncol=2, loc="center", bbox_to_anchor=(0.47, 0.52), frameon=False)
    fig.tight_layout(pad=0.4)
    fig.savefig(FIG / "fig1_continuation.pdf")
    plt.close(fig)


# ----------------------------------------------------------------------------- Fig. 2 tangent landscape
def fig2():
    N = 21
    s = np.arange(-10, 11, dtype=float) / N
    S2, S4 = np.sum(s * s), np.sum(s ** 4)

    def D(v):
        return np.cos(np.multiply.outer(np.asarray(v), s)).sum(axis=-1)

    def D1(v):
        return -(np.sin(np.multiply.outer(np.asarray(v), s)) * s).sum(axis=-1)

    def D2(v):
        return -(np.cos(np.multiply.outer(np.asarray(v), s)) * s * s).sum(axis=-1)

    def R(v, lam):
        v = np.asarray(v)
        q0 = -S2 - lam * D(10)
        q1 = 2 * lam * D1(10)
        C = S4 - 2 * lam * D2(10) + N * lam * lam - q0 * q0 / N - q1 * q1 / (4 * S2)
        G = N - D(v) ** 2 / N - D1(v) ** 2 / S2
        H = D2(v) - lam * D(10 - v) - D(v) * q0 / N + D1(v) * q1 / (2 * S2)
        with np.errstate(divide="ignore", invalid="ignore"):
            out = C - H * H / G
        return np.where(np.abs(v) < 0.15, np.nan, out)

    lc = float(J(ROOT / "stage3_tangent_reference.json")["lambda_N"])
    lamvals = [lc - 0.002, lc, lc + 0.002]
    grid = np.linspace(-7, 16, 3000)
    fig, ax = plt.subplots(1, 3, figsize=(DOUBLE, 2.2), sharex=True, sharey=True)
    for k, (a, lam) in enumerate(zip(ax, lamvals)):
        mA = minimize_scalar(lambda v: float(R(v, lam)), bounds=(-5, -3), method="bounded")
        mB = minimize_scalar(lambda v: float(R(v, lam)), bounds=(11, 14), method="bounded")
        base = min(mA.fun, mB.fun)
        a.plot(grid, R(grid, lam) - base, color=C1, lw=1.4)
        a.plot([mA.x], [mA.fun - base], "o", ms=5, color=INK)
        a.plot([mB.x], [mB.fun - base], "s", ms=5, color=INK)
        a.text(mA.x - 0.9, mA.fun - base + 0.006, BA, fontsize=9)
        a.text(mB.x - 0.4, mB.fun - base + 0.006, BB, fontsize=9)
        a.set_title([r"$\lambda=\lambda_{21}-0.002$", r"$\lambda=\lambda_{21}$", r"$\lambda=\lambda_{21}+0.002$"][k])
        a.set_xlabel(r"Satellite coordinate $v$")
        a.set_xlim(-7, 16)
        a.set_ylim(-0.001, 0.07)
        style(a)
    ax[0].set_ylabel(r"$R_{21}(v,\lambda)-\min_v R_{21}$")
    fig.tight_layout(pad=0.4)
    fig.savefig(FIG / "fig2_tangent_landscape.pdf")
    plt.close(fig)


# ----------------------------------------------------------------------------- Fig. 3 fixed-geometry large N
def fig3():
    rows = list(csv.DictReader((ROOT / "lambdaN_asymptotics.csv").open()))
    co = J(ROOT / "stage14_asymptotic_coefficients.json")
    l0, a2, a4 = float(co["lambda_infinity"]), float(co["a_1_over_N2"]), float(co["b_1_over_N4"])
    x = np.array([1 / int(r["N"]) ** 2 for r in rows])
    y = np.array([float(r["lambda_N"]) for r in rows])
    cert = np.array([r["origin"].startswith("existing certified") for r in rows])
    grid = np.linspace(0, 1 / 11 ** 2 * 1.02, 200)
    fig, ax = plt.subplots(figsize=(SINGLE, 2.35))
    ax.plot(grid, l0 + a2 * grid + a4 * grid ** 2, color=MUTED, lw=1.1, label=r"$\lambda_\infty+a_2N^{-2}+a_4N^{-4}$")
    ax.plot(x[cert], y[cert], "o", ms=5, color=C1, label="interval-certified $N=11,15,21,31,41$")
    ax.plot(x[~cert], y[~cert], "^", ms=5, color=C2, label="numerical check $N=101,\\ldots,1001$")
    ax.plot([0], [l0], "*", ms=9, color=INK, label=r"certified continuum $\lambda_\infty$")
    ax.set(xlabel=r"$1/N^2$ (with $z=N\Delta$ and $b=N\omega_3$ held fixed)", ylabel=r"Critical coefficient $\lambda_N$")
    style(ax)
    ax.legend(frameon=False, loc="lower left")
    fig.tight_layout(pad=0.4)
    fig.savefig(FIG / "fig3_largeN_fixed_geometry.pdf")
    plt.close(fig)


# ----------------------------------------------------------------------------- Fig. 4 z law
def fig4():
    rows = list(csv.DictReader((ROOT / "z_continuation_vs_asymptotics.csv").open()))
    z = np.array([float(r["z"]) for r in rows])
    fig, ax = plt.subplots(figsize=(SINGLE, 2.4))
    ax.plot(z, [float(r["numerical_epsilon"]) for r in rows], color=C1, lw=2.0,
            label=r"numerical $\mathcal{A}/\mathcal{B}$ equal-cost continuation")
    ax.plot(z, [float(r["quadratic_epsilon"]) for r in rows], color=C2, ls="--", lw=1.3, label=r"$\lambda_{21}z^2$")
    ax.plot(z, [float(r["quartic_epsilon"]) for r in rows], color=C3, ls=":", lw=1.6, label=r"$\lambda_{21}z^2+c_{21}z^4$")
    ax.plot([2], [0.2481906301722774], "s", ms=5, color=INK, label=r"independently certified crossing, $z=2$")
    ax.set(xlabel=r"Normalized spacing $z=N\Delta$", ylabel=r"Crossing amplitude $\epsilon_c$", xlim=(0, 2.1), ylim=(0, 0.30))
    style(ax)
    ax.legend(loc="upper left", frameon=False)
    fig.tight_layout(pad=0.4)
    fig.savefig(FIG / "fig4_z_law.pdf")
    plt.close(fig)


# ----------------------------------------------------------------------------- Fig. 6 phase regimes (printed as Fig. 6)
def fig5():
    d = J(VAL / "stage2_phase/DONE.json")
    tan = {float(k): v for k, v in d["tangent"].items() if v}
    fin = list(csv.DictReader((VAL / "stage2_phase/stage2A_finite.csv").open()))
    fig, ax = plt.subplots(2, 1, figsize=(SINGLE, 4.3))
    ks = sorted(k for k in tan if k >= 0)
    ax[0].plot(ks, [tan[k]["lam"] for k in ks], "-", color=INK, lw=1.2, label=r"tangent $\lambda_{21}(\chi)$")
    kk = np.linspace(0, 0.2, 40)
    ax[0].plot(kk, tan[0.0]["lam"] + 1.92838866512018 * kk ** 2, ls="--", color=MUTED, lw=1.0,
               label=r"$\lambda_{21}(0)+1.928\,\chi^2$")
    for z, col, mk in ((0.25, C1, "o"), (0.5, C2, "s"), (1.0, C3, "^")):
        rr = [r for r in fin if abs(float(r["z"]) - z) < 1e-12 and float(r["kappa"]) >= 0 and r.get("eps_c_over_z2")]
        rr.sort(key=lambda r: float(r["kappa"]))
        q = [r for r in rr if r["verdict"] == "present"]
        n = [r for r in rr if r["verdict"] != "present"]
        ax[0].plot([float(r["kappa"]) for r in q], [float(r["eps_c_over_z2"]) for r in q], mk, ms=4.5, color=col,
                   label=fr"finite $\epsilon_c/z^2$, $z={z:g}$ (separated)")
        if n:
            ax[0].plot([float(r["kappa"]) for r in n], [float(r["eps_c_over_z2"]) for r in n], mk, ms=4.5, mfc="white",
                       color=col)
    ax[0].set(xlabel=r"Offset scale $\chi=\delta/z=\tau/N$ ($\tau$: centre shift in samples)",
              ylabel=r"Crossing coefficient")
    ax[0].text(0.31, 0.10, "open: lower fit\nnear-coalescent", fontsize=6.5, color=MUTED)
    style(ax[0])
    ax[0].legend(frameon=False, loc="upper left", fontsize=6.2)
    tag(ax[0], "(a)")
    # (b) z = 2 topology
    pts = []
    for f in sorted((VAL / "stage2_phase/analysis").glob("S2A_kP*_z2p0.json")):
        a = J(f)
        c = a["config"]
        ex = [s for s in a["switches"] if s.get("type") == "exchange"]
        if not ex:
            continue
        s = ex[0]
        pts.append((c["delta"], s["P"]["sep"], s["status"] == "present"))
    pts.sort()
    dl = np.array([p[0] for p in pts])
    sep = np.array([p[1] for p in pts])
    ok = np.array([p[2] for p in pts])
    ax[1].plot(dl, sep, "-", color=MUTED, lw=0.8)
    ax[1].plot(dl[ok], sep[ok], "o", ms=5, color=C1, label="separated exchange")
    ax[1].plot(dl[~ok], sep[~ok], "o", ms=5, mfc="white", color=C2, label="lower fit confluent (not separated)")
    ax[1].axhline(1.0, color=MUTED, ls=":", lw=0.8)
    ax[1].text(0.02, 1.15, "separation threshold", fontsize=6.5, color=MUTED)
    ax[1].set(xlabel=r"Strong-pair offset $\delta$ at $z=2$ (rad)",
              ylabel=r"Separation $|u_2-u_1|$")
    style(ax[1])
    ax[1].legend(frameon=False, loc="upper right")
    tag(ax[1], "(b)")
    fig.tight_layout(pad=0.4, h_pad=0.8)
    fig.savefig(FIG / "fig6_phase_regimes.pdf")
    plt.close(fig)


# ----------------------------------------------------------------------------- Fig. 5 robustness (printed as Fig. 5)
def fig6():
    s1 = J(VAL / "stage1_window/DONE.json")
    s3 = J(VAL / "stage3_bwidth/DONE_3A.json")
    s4 = J(VAL / "stage4_imbalance/DONE.json")
    fig, ax = plt.subplots(1, 3, figsize=(DOUBLE, 2.15), gridspec_kw=dict(width_ratios=[1.05, 1, 1]))
    wins = [("rect", "Rect."), ("hann", "Hann"), ("hamming", "Hamming"), ("dpss", "DPSS")]
    rows = {r["window"]: r for r in s1["rows"] if r["role"] == "primary"}
    x = np.arange(len(wins))
    ax[0].bar(x, [rows[w]["eps_c"] for w, _ in wins], width=0.55, color=C1, edgecolor="white", linewidth=1.5)
    for i, (w, _) in enumerate(wins):
        ax[0].text(i, rows[w]["eps_c"] + 0.006, f"{rows[w]['eps_c']:.3f}", ha="center", fontsize=6.5, color=INK)
    ax[0].set_xticks(x, [l for _, l in wins])
    ax[0].set(ylabel=r"Crossing amplitude $\epsilon_c$", ylim=(0, 0.29), title="Weighting window")
    style(ax[0])
    tag(ax[0], "(a)")
    fr = [r for r in s3["finite_N21_z2"] if "eps_c" in r]
    ax[1].plot([r["b"] for r in fr], [r["eps_c"] for r in fr], "o-", ms=4, color=C1)
    ax[1].plot([10.0], [next(r["eps_c"] for r in fr if r["b"] == 10.0)], "s", ms=6, color=INK, label="certified $b=10$")
    ax[1].set(xlabel=r"Weak-tone location $b$", title="Weak-tone location", ylim=(0.2, 0.27))
    ax[1].legend(frameon=False, loc="lower right")
    style(ax[1])
    tag(ax[1], "(b)")
    rr = sorted([r for r in s4["rows"] if "eps_c" in r], key=lambda r: r["eta"])
    ax[2].plot([r["eta"] for r in rr], [r["eps_c"] for r in rr], "o-", ms=4, color=C1)
    ax[2].plot([0.0], [next(r["eps_c"] for r in rr if r["eta"] == 0.0)], "s", ms=6, color=INK, label="equal amplitudes")
    ax[2].set(xlabel=r"Imbalance $\iota$ (amplitudes $1\mp\iota$)", title="Strong-pair imbalance", ylim=(0.2, 0.27))
    ax[2].legend(frameon=False, loc="lower right")
    style(ax[2])
    tag(ax[2], "(c)")
    fig.tight_layout(pad=0.4, w_pad=1.0)
    fig.savefig(FIG / "fig5_robustness.pdf")
    plt.close(fig)


# ----------------------------------------------------------------------------- Fig. 7 noise
def fig7():
    summ = list(csv.DictReader((VAL / "stage6_noise/summary.csv").open()))
    col = {20: C1, 30: C2, 40: C3}
    mk = {20: "o", 30: "s", 40: "^"}
    fig, ax = plt.subplots(1, 3, figsize=(DOUBLE, 2.55))
    for snr in (20, 30, 40):
        ss = sorted([r for r in summ if int(r["snr"]) == snr], key=lambda r: float(r["eta"]))
        e = np.array([float(r["eta"]) for r in ss])
        p = np.array([float(r["P_A"]) for r in ss])
        lo = np.clip(p - np.array([float(r["P_A_lo"]) for r in ss]), 0, None)
        hi = np.clip(np.array([float(r["P_A_hi"]) for r in ss]) - p, 0, None)
        ax[0].errorbar(e, p, yerr=[lo, hi], fmt=mk[snr] + "-", ms=3.5, lw=1.1, color=col[snr], capsize=1.5,
                       label=f"{snr} dB")
    ax[0].set(xlabel=r"$\eta=(\epsilon-\epsilon_c)/\epsilon_c$", ylabel=r"$P(\mathrm{select}\ \mathcal{A})$",
              title="Branch selection")
    ax[0].legend(frameon=False, loc="upper right")
    style(ax[0])
    tag(ax[0], "(a)")
    for a, snr, t in ((ax[1], 30, "(b)"), (ax[2], 40, "(c)")):
        ss = sorted([r for r in summ if int(r["snr"]) == snr], key=lambda r: float(r["eta"]))
        e = [float(r["eta"]) for r in ss]

        def val(r, k):
            return float(r[k]) if r.get(k) not in (None, "") else np.nan
        a.semilogy(e, [val(r, "MSE_to_winner") for r in ss], "-", color=INK, lw=1.6, label="global MSE (all records)")
        a.semilogy(e, [val(r, "cond_MSE_own_A") if int(float(r.get("cond_n_A") or 0)) >= 20 else np.nan for r in ss],
                   "o", ms=3.5, color=C1, label=BA + "-conditioned MSE")
        a.semilogy(e, [val(r, "cond_MSE_own_B") if int(float(r.get("cond_n_B") or 0)) >= 20 else np.nan for r in ss],
                   "s", ms=3.5, color=C2, label=BB + "-conditioned MSE")
        a.semilogy(e, [val(r, "local_trace_A") for r in ss], "--", color=C1, lw=1.0, label="local approx., " + BA)
        a.semilogy(e, [val(r, "local_trace_B") for r in ss], "--", color=C2, lw=1.0, label="local approx., " + BB)
        a.set(xlabel=r"$\eta$", title=f"Squared error in $u$, {snr} dB")
        style(a)
        tag(a, t)
    ax[1].set_ylabel("Frequency-pair squared error")
    h, l = ax[1].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=3, frameon=False, fontsize=6.5, bbox_to_anchor=(0.66, -0.005))
    fig.tight_layout(pad=0.4, w_pad=0.8, rect=(0, 0.13, 1, 1))
    fig.savefig(FIG / "fig7_noise_local_global.pdf")
    plt.close(fig)


# ----------------------------------------------------------------------------- supplementary figures
def sup_figs():
    import shutil
    shutil.copy(ROOT / "paper/figures/fig7_phase_resolution.pdf", SUP / "figS1_weak_phase_resolution.pdf")
    shutil.copy(VAL / "stage1_window/window_branches.pdf", SUP / "figS2_window_branches.pdf")
    # noise distributions and 20 dB errors
    trials = list(csv.DictReader((VAL / "stage6_noise/trials.csv").open()))
    roots = {float(k): v for k, v in J(VAL / "stage6_noise/noiseless_roots.json").items()}
    r0 = roots[0.0]
    dvec = np.array(r0["B"]) - np.array(r0["A"])
    dvec /= np.linalg.norm(dvec)
    fig, ax = plt.subplots(1, 4, figsize=(DOUBLE, 2.2))
    for j, snr in enumerate((20, 30, 40)):
        rr = [r for r in trials if int(r["snr"]) == snr and abs(float(r["eta"])) < 1e-12]
        U = np.array([[float(r["u1"]), float(r["u2"])] for r in rr])
        cvals = (U - np.array(r0["A"])) @ dvec
        ax[j].hist(cvals, bins=60, color=C1, edgecolor="white", linewidth=0.3)
        ax[j].set(xlabel=r"coordinate along $\mathcal{A}\to\mathcal{B}$", title=f"{snr} dB, $\\eta=0$")
        style(ax[j])
    ax[0].set_ylabel("records")
    summ = list(csv.DictReader((VAL / "stage6_noise/summary.csv").open()))
    ss = sorted([r for r in summ if int(r["snr"]) == 20], key=lambda r: float(r["eta"]))
    e = [float(r["eta"]) for r in ss]
    ax[3].semilogy(e, [float(r["MSE_to_winner"]) for r in ss], "-", color=INK, label="global MSE")
    ax[3].semilogy(e, [float(r["cond_MSE_own_A"]) for r in ss], "o", ms=3, color=C1, label=BA + "-conditioned")
    ax[3].semilogy(e, [float(r["cond_MSE_own_B"]) for r in ss], "s", ms=3, color=C2, label=BB + "-conditioned")
    ax[3].semilogy(e, [float(r["local_trace_A"]) for r in ss], "--", color=C1, lw=1.0, label="local, " + BA)
    ax[3].semilogy(e, [float(r["local_trace_B"]) for r in ss], "--", color=C2, lw=1.0, label="local, " + BB)
    ax[3].set(xlabel=r"$\eta$", title="20 dB errors", ylim=(0.03, 600))
    ax[3].legend(frameon=False, fontsize=5.2, loc="center", ncol=2, columnspacing=0.8, handlelength=1.6)
    style(ax[3])
    fig.tight_layout(pad=0.4)
    fig.savefig(SUP / "figS4_noise_distributions.pdf")
    plt.close(fig)
    # post-hoc fixed-delta scaling
    ph = J(VAL / "stage2_phase/POSTHOC_fixed_delta_all_switches.json")
    fig, ax = plt.subplots(figsize=(SINGLE, 2.4))
    cols = [C1, C2, C3, INK, MUTED]
    for (key, c) in zip(("0.05", "0.1", "0.2", "0.3", "0.4"), cols):
        rows = [r for r in ph[key]["rows"] if "eps_c" in r]
        z = np.array([r["z"] for r in rows])
        eps = np.array([r["eps_c"] for r in rows])
        reg = np.array([r["P_cls"] == "regular" for r in rows])
        ax.loglog(z, eps, "-", color=c, lw=0.9, label=fr"$\delta/\pi={key}$, $\alpha={ph[key]['fits']['full']['alpha']:.3f}$")
        ax.loglog(z[reg], eps[reg], "o", ms=4, color=c)
        ax.loglog(z[~reg], eps[~reg], "o", ms=4, mfc="white", color=c)
    ax.set(xlabel=r"$z$", ylabel=r"$\epsilon$ at the switch into the $\mathcal{B}$-like fit")
    ax.set_title("POST-HOC sensitivity (not preregistered)", fontsize=7)
    style(ax)
    ax.legend(frameon=False, fontsize=5.8, loc="upper left")
    fig.tight_layout(pad=0.4)
    fig.savefig(SUP / "figS3_posthoc_fixed_delta.pdf")
    plt.close(fig)


if __name__ == "__main__":
    for f in (fig1, fig2, fig3, fig4, fig5, fig6, fig7, sup_figs):
        f()
        print("done", f.__name__)
