"""Small fixed-design complex-Gaussian BIC sanity check (not a theorem)."""
import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize

from three_to_two_tone_stage5_core import grid_minima, refine, signal, times

ROOT = Path(__file__).resolve().parent
N, Z, B, EC = 21, 2.0, 10.0, 0.2481906301722774
S = times(N)
A = np.array([-4.772447353918614, 0.605170320276862])
BB = np.array([-0.057653561926680, 12.46341935680166])
RNG = np.random.default_rng(20260928)
ETAS = (-0.10, 0.0, 0.10)
SNRS = (20, 30, 40)


def vp3(u, y):
    u = np.asarray(u, float)
    v = np.exp(1j * np.outer(S, u))
    coef = np.linalg.lstsq(v, y, rcond=1e-12)[0]
    residual = y-v@coef
    rss = float(np.vdot(residual, residual).real)
    grad = -2*np.real(np.sum(np.conj(residual[:, None]) *
                             (1j*S[:, None]*v)*coef[None, :], axis=0))
    return rss, grad


def fit3(y, u2):
    # Continuous off-grid variable projection from fixed, untuned seed rules.
    frequencies = np.linspace(-np.pi*N, np.pi*N, 257, endpoint=False)
    v = np.exp(1j*np.outer(S, frequencies))
    v2 = np.exp(1j*np.outer(S, u2))
    residual2 = y-v2@np.linalg.lstsq(v2, y, rcond=None)[0]
    peak = float(frequencies[np.argmax(np.abs(v.conj().T@residual2))])
    seeds = [np.array([-Z, Z, B]), np.r_[u2, peak],
             np.r_[A, B], np.r_[BB, -Z], np.r_[BB, Z]]
    best = None
    for seed in seeds:
        result = minimize(lambda u: vp3(u, y), seed, jac=True, method="BFGS",
                          options={"gtol": 1e-7, "maxiter": 180})
        rss = vp3(result.x, y)[0]
        if best is None or rss < best[0]:
            best = rss, np.sort(result.x), float(np.linalg.norm(vp3(result.x, y)[1]))
    return best


def quantile(vals, p):
    return float(np.quantile(np.asarray(vals, float), p))


def run(trials):
    out = []
    for snr in SNRS:
        for eta in ETAS:
            epsilon = EC*(1+eta)
            x = signal(N, Z, B, epsilon)
            sigma2 = np.vdot(x, x).real/N * 10**(-snr/10)
            for trial in range(trials):
                noise = np.sqrt(sigma2/2)*(RNG.standard_normal(N)+1j*RNG.standard_normal(N))
                y = x+noise
                fit_a = refine(A, y, S)
                fit_b = refine(BB, y, S)
                fit_strong = refine([-Z, Z], y, S)
                fit2 = min((fit_a, fit_b, fit_strong), key=lambda f: f["J"])
                audited = trial % 20 == 0
                audit_improved = False
                if audited:
                    grid = grid_minima(y, S, grid_size=128, top=16)
                    if grid and grid[0]["J"] < fit2["J"]:
                        audit_improved = grid[0]["J"] < fit2["J"]-1e-7
                        fit2 = grid[0]
                rss3, u3, grad3 = fit3(y, fit2["u"])
                rss2 = fit2["J"]
                # 2N real observations; 3K real mean parameters for K complex tones.
                dbic = 2*N*np.log(rss3/rss2)+3*np.log(2*N)
                out.append(dict(snr_db=snr, eta=eta, epsilon=epsilon, trial=trial,
                                sigma2=sigma2, rss2=rss2, rss3=rss3,
                                delta_bic_3_minus_2=dbic,
                                selected_order=3 if dbic < 0 else 2,
                                two_tone_audited=int(audited),
                                two_tone_grid_improved=int(audit_improved),
                                two_tone_u1=fit2["u"][0], two_tone_u2=fit2["u"][1],
                                three_tone_u1=u3[0], three_tone_u2=u3[1],
                                three_tone_u3=u3[2], three_tone_gradient_norm=grad3))
            print("finished", snr, eta, flush=True)
    with (ROOT/"order_selection_sanity.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)

    summary=[]
    for snr in SNRS:
        for eta in ETAS:
            rr=[r for r in out if r["snr_db"]==snr and r["eta"]==eta]
            dd=[r["delta_bic_3_minus_2"] for r in rr]
            n3=sum(r["selected_order"]==3 for r in rr)
            summary.append(dict(snr_db=snr,eta=eta,n=len(rr),n_order2=len(rr)-n3,
                                n_order3=n3,fraction_order2=(len(rr)-n3)/len(rr),
                                fraction_order3=n3/len(rr),
                                delta_bic_q05=quantile(dd,.05),
                                delta_bic_q25=quantile(dd,.25),
                                delta_bic_median=quantile(dd,.5),
                                delta_bic_q75=quantile(dd,.75),
                                delta_bic_q95=quantile(dd,.95)))
    with (ROOT/"order_selection_sanity_summary.csv").open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(summary[0]))
        w.writeheader();w.writerows(summary)

    colors={-.1:"#267a4c",0.:"#116ba2",.1:"#b35b12"}
    fig, axes=plt.subplots(1,3,figsize=(8.8,3.2),sharey=True)
    for ax,snr in zip(axes,SNRS):
        for eta in ETAS:
            vals=np.sort([r["delta_bic_3_minus_2"] for r in out if r["snr_db"]==snr and r["eta"]==eta])
            ax.step(vals,np.arange(1,len(vals)+1)/len(vals),where="post",lw=1.5,
                    color=colors[eta],label=rf"$\eta={eta:+.1f}$")
        ax.axvline(0,color="black",ls="--",lw=1)
        ax.set_title(f"{snr} dB")
        ax.set_xlabel(r"$\mathrm{BIC}_3-\mathrm{BIC}_2$")
        ax.grid(alpha=.18)
    axes[0].set_ylabel("Empirical CDF")
    axes[0].legend(fontsize=7,frameon=False)
    fig.tight_layout()
    fig.savefig(ROOT/"order_selection_sanity.pdf")
    plt.close(fig)
    for rec in summary:
        print(rec)
    print("max three-tone gradient",max(r["three_tone_gradient_norm"] for r in out))


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--trials",type=int,default=120)
    run(parser.parse_args().trials)
