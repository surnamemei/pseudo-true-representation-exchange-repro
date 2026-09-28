"""Deterministic predictions of the branch-selection / global-MSE transition theory.

Governing plan: paper/transition_theory/00_THEORY_PLAN.md (sealed, SEAL_1). This script reads only deterministic inputs:
the frozen Stage-6 code (advcore.py, s6_noise.py) and the noiseless branch continuation
(stage6_noise/noiseless_roots.json). It never opens a Monte Carlo outcome file.

Outputs (results/transition_theory/deterministic/): predictions.csv, deterministic.json, checks.json.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import platform
import sys
import time
from pathlib import Path

sys.dont_write_bytecode = True  # keep the frozen validation folder free of new byte-code files
ROOT = Path(__file__).resolve().parents[3]
VAL = ROOT / "validation/adversarial_overnight"
sys.path.insert(0, str(VAL))

import mpmath as mp  # noqa: E402
import numpy as np  # noqa: E402
import scipy  # noqa: E402
from scipy.optimize import brentq  # noqa: E402
from scipy.special import ndtr, ndtri  # noqa: E402

import s6_noise  # noqa: E402  (frozen Stage-6 module: constants, local_cov, noiseless_roots)
from advcore import make_signal, times, refine_from, Objective  # noqa: E402

OUT = ROOT / "results/transition_theory/deterministic"
SPEC = json.loads((ROOT / "paper/transition_theory/transition_spec.json").read_text())
N, Z, BW, EC = s6_noise.N, s6_noise.Z, s6_noise.BW, s6_noise.EC
ETAS, SNRS = list(s6_noise.ETAS), list(s6_noise.SNRS)
PHI_INV_09 = SPEC["theory"]["Phi_inv_0p9"]
ROOTS_FILE = VAL / "stage6_noise/noiseless_roots.json"


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def check_inputs():
    bad = {rel: h for rel, h in SPEC["deterministic_inputs"].items() if sha(ROOT / rel) != h}
    if bad:
        raise SystemExit(f"deterministic inputs differ from the sealed spec: {sorted(bad)}")
    assert (s6_noise.N, s6_noise.Z, s6_noise.BW, s6_noise.EC) == (21, 2.0, 10.0, SPEC["setting"]["eps_c"])
    assert ETAS == SPEC["setting"]["eta_grid"] and SNRS == SPEC["setting"]["snr_db"]
    assert abs(float(ndtri(0.9)) - PHI_INV_09) < 1e-15


# ----------------------------------------------------------------------------- float64 building blocks
S = times(N)
DXDE = np.exp(1j * math.pi) * np.exp(1j * BW * S)          # dx/d eps, exactly as make_signal builds the weak tone


def residual(x, u):
    """Least-squares residual of the two-tone fit at frequencies u (as in s6_noise.local_cov)."""
    v = np.exp(1j * np.outer(S, np.asarray(u, float)))
    a = np.linalg.lstsq(v, x, rcond=None)[0]
    return x - v @ a


def sigma2_of(x, snr):
    return float(np.vdot(x, x).real / N) * 10 ** (-snr / 10)


def continue_branches(eps, steps=40, tol_g=1e-12):
    """Continue A and B from the frozen crossing roots to eps (same scheme as s6_noise.noiseless_roots)."""
    a, b = s6_noise.A_FROZEN, s6_noise.B_FROZEN
    for j in range(1, steps + 1):
        ee = EC + (eps - EC) * j / steps
        qa, qb = refine_from(Objective(make_signal(N, Z, BW, ee)), [a, b], tol_g=tol_g)
        a, b = qa.pair(), qb.pair()
    qa, qb = refine_from(Objective(make_signal(N, Z, BW, eps)), [a, b], tol_g=tol_g)
    return qa, qb


def branch_quantities(eps, uA, uB, JA=None, JB=None):
    x = make_signal(N, Z, BW, eps)
    rA, rB = residual(x, uA), residual(x, uB)
    d = rA - rB
    JA_r, JB_r = float(np.vdot(rA, rA).real), float(np.vdot(rB, rB).real)
    return dict(x=x, rA=rA, rB=rB, d=d, d_norm=float(np.linalg.norm(d)),
                JA=JA_r if JA is None else JA, JB=JB_r if JB is None else JB, JA_resid=JA_r, JB_resid=JB_r,
                x_norm2=float(np.vdot(x, x).real))


def p_exact(Delta0, sigma2, d_norm):
    return float(ndtr(-Delta0 / (math.sqrt(2 * sigma2) * d_norm)))


# ----------------------------------------------------------------------------- 50-digit mpmath building blocks
mp.mp.dps = 50
S_MP = [mp.mpf(t - (N - 1) // 2) / N for t in range(N)]


def x_mp(eps):
    w = mp.expj(mp.pi)
    return [mp.expj(-Z * s) + mp.expj(Z * s) + eps * w * mp.expj(BW * s) for s in S_MP]


def varpro_mp(u, y):
    v = [[mp.expj(uk * s) for uk in u] for s in S_MP]
    G = mp.matrix(2, 2)
    rhs = mp.matrix(2, 1)
    for t in range(N):
        for j in range(2):
            rhs[j] += mp.conj(v[t][j]) * y[t]
            for k in range(2):
                G[j, k] += mp.conj(v[t][j]) * v[t][k]
    c = mp.lu_solve(G, rhs)
    r = [y[t] - v[t][0] * c[0] - v[t][1] * c[1] for t in range(N)]
    return mp.fsum(abs(rt) ** 2 for rt in r), r, c, v


def grad_mp(u, y):
    """dJ/du_k = -2 Im(conj(c_k) sum_t s_t conj(v_tk) r_t) (envelope theorem over the amplitudes)."""
    J, r, c, v = varpro_mp(u, y)
    return [-2 * mp.im(mp.conj(c[k]) * mp.fsum(S_MP[t] * mp.conj(v[t][k]) * r[t] for t in range(N))) for k in range(2)]


def polish_mp(u0, y):
    sol = mp.findroot(lambda a, b: grad_mp([a, b], y), (mp.mpf(u0[0]), mp.mpf(u0[1])), tol=mp.mpf(10) ** -90,
                      maxsteps=60)
    u = [sol[0], sol[1]]
    J, r, c, v = varpro_mp(u, y)
    gn = max(abs(gk) for gk in grad_mp(u, y))
    return u, J, r, gn


# ----------------------------------------------------------------------------- main computation
def main():
    t0 = time.time()
    check_inputs()
    OUT.mkdir(parents=True, exist_ok=True)
    roots = {float(k): v for k, v in json.loads(ROOTS_FILE.read_text()).items()}
    checks = {}

    # check 5: reproduce the frozen continuation (float64)
    rep = s6_noise.noiseless_roots()
    jdiff = max(abs(rep[e]["J_A"] - roots[e]["J_A"]) / abs(roots[e]["J_A"]) for e in ETAS)
    jdiff = max(jdiff, max(abs(rep[e]["J_B"] - roots[e]["J_B"]) / abs(roots[e]["J_B"]) for e in ETAS))
    udiff = max(float(np.max(np.abs(np.array(rep[e][k]) - np.array(roots[e][k])))) for e in ETAS for k in ("A", "B"))
    checks["roots_reproduced"] = dict(max_rel_diff_J=jdiff, max_abs_diff_u=udiff, pass_=bool(jdiff <= 1e-9 and udiff <= 1e-7))

    # ---- crossing quantities (float64)
    rc = roots[0.0]
    bq = branch_quantities(EC, rc["A"], rc["B"], rc["J_A"], rc["J_B"])
    d_c, dn_c = bq["d"], bq["d_norm"]
    g_env = float(2 * np.vdot(d_c, DXDE).real)
    gA_env = float(2 * np.vdot(bq["rA"], DXDE).real)
    gB_env = float(2 * np.vdot(bq["rB"], DXDE).real)

    def gap(eps):
        qa, qb = refine_from(Objective(make_signal(N, Z, BW, eps)), [rc["A"], rc["B"]], tol_g=1e-13)
        return qa.J - qb.J

    def fd(h):
        return (gap(EC + h) - gap(EC - h)) / (2 * h)

    h = 1e-4 * EC
    D1, D2 = fd(h), fd(h / 2)
    g_fd = (4 * D2 - D1) / 3
    checks["g_float64"] = dict(g_envelope=g_env, g_fd_h=D1, g_fd_h2=D2, g_fd_richardson=g_fd,
                               rel_diff=abs(g_fd - g_env) / abs(g_env), pass_=bool(abs(g_fd - g_env) <= 1e-6 * abs(g_env)))

    # ---- 50-digit crossing check
    y0 = x_mp(mp.mpf(EC))
    uA0, JA0, rA0, gnA = polish_mp(rc["A"], y0)
    uB0, JB0, rB0, gnB = polish_mp(rc["B"], y0)
    dmp = [rA0[t] - rB0[t] for t in range(N)]
    dxde_mp = [mp.expj(mp.pi) * mp.expj(BW * s) for s in S_MP]
    g_env_mp = 2 * mp.re(mp.fsum(mp.conj(dmp[t]) * dxde_mp[t] for t in range(N)))
    dn_mp = mp.sqrt(mp.fsum(abs(dt) ** 2 for dt in dmp))
    hmp = mp.mpf(10) ** -12 * mp.mpf(EC)

    def gap_mp(eps):
        y = x_mp(eps)
        return polish_mp(uA0, y)[1] - polish_mp(uB0, y)[1]

    g_fd_mp = (gap_mp(mp.mpf(EC) + hmp) - gap_mp(mp.mpf(EC) - hmp)) / (2 * hmp)
    rel_mp = abs(g_fd_mp - g_env_mp) / abs(g_env_mp)
    # high-precision crossing root by Newton on the gap, using the envelope slope
    e_star = mp.mpf(EC)
    for _ in range(3):
        y = x_mp(e_star)
        ua, ja, ra, _ = polish_mp(uA0, y)
        ub, jb, rb, _ = polish_mp(uB0, y)
        slope = 2 * mp.re(mp.fsum(mp.conj(ra[t] - rb[t]) * dxde_mp[t] for t in range(N)))
        e_star = e_star - (ja - jb) / slope
    checks["g_mpmath"] = dict(dps=mp.mp.dps, grad_norm_A=float(gnA), grad_norm_B=float(gnB),
                              g_envelope=mp.nstr(g_env_mp, 30), g_fd=mp.nstr(g_fd_mp, 30), rel_diff=float(rel_mp),
                              pass_=bool(rel_mp <= mp.mpf(10) ** -20))
    checks["g_float64_vs_mpmath"] = dict(rel_diff=abs(g_env - float(g_env_mp)) / float(g_env_mp))
    checks["d_norm"] = dict(float64=dn_c, mpmath=mp.nstr(dn_mp, 30), rel_diff=abs(dn_c - float(dn_mp)) / float(dn_mp),
                            pass_=bool(abs(dn_c - float(dn_mp)) <= 1e-9 * float(dn_mp)))

    # ---- check 3: synthetic Gaussian variance of the linearized gap (fixed d_c; no optimizer, no outcomes)
    sig2_30 = sigma2_of(bq["x"], 30)
    rng = np.random.default_rng(np.random.SeedSequence([20260928, 99, 1]))
    M, chunk = 2_000_000, 200_000
    tot = tot2 = 0.0
    for _ in range(M // chunk):
        w = math.sqrt(sig2_30 / 2) * (rng.normal(size=(chunk, N)) + 1j * rng.normal(size=(chunk, N)))
        T = 2 * (w @ np.conj(d_c)).real          # 2 Re <d, w> = 2 Re sum conj(d) w
        tot += T.sum()
        tot2 += (T * T).sum()
    mean = tot / M
    var = tot2 / M - mean * mean
    theory = 2 * sig2_30 * dn_c ** 2
    checks["variance_formula"] = dict(M=M, sigma2=sig2_30, sample_mean=mean, sample_var=var, theory_var=theory,
                                      var_ratio_minus_1=var / theory - 1, tol_var=5 * math.sqrt(2 / M),
                                      mean_over_sd=abs(mean) / math.sqrt(theory), tol_mean=5 / math.sqrt(M),
                                      pass_=bool(abs(var / theory - 1) <= 5 * math.sqrt(2 / M)
                                                 and abs(mean) / math.sqrt(theory) <= 5 / math.sqrt(M)))

    # ---- closed-form constants per SNR
    per_snr = {}
    for snr in SNRS:
        s2c = sigma2_of(bq["x"], snr)
        K = g_env * EC / (math.sqrt(2 * s2c) * dn_c)
        per_snr[snr] = dict(sigma2_c=s2c, sigma_c=math.sqrt(s2c), s_Delta_c=math.sqrt(2 * s2c) * dn_c, K=K,
                            W_closed_form=2 * PHI_INV_09 / abs(K))

    # ---- per-eta predictions
    rows = []
    dnorms = {}
    for eta in ETAS:
        R = roots[eta]
        eps = R["eps"]
        q = branch_quantities(eps, R["A"], R["B"], R["J_A"], R["J_B"])
        Delta0 = R["J_A"] - R["J_B"]
        dnorms[eta] = q["d_norm"]
        uA, uB = np.array(R["A"]), np.array(R["B"])
        uW = uA if eta < 0 else uB
        bA, bB = float(np.sum((uA - uW) ** 2)), float(np.sum((uB - uW) ** 2))
        for snr in SNRS:
            s2 = sigma2_of(q["x"], snr)
            covA = s6_noise.local_cov(q["x"], uA, s2)
            covB = s6_noise.local_cov(q["x"], uB, s2)
            trA, trB = float(np.trace(covA[0])), float(np.trace(covB[0]))
            PA = p_exact(Delta0, s2, q["d_norm"])
            PAl = float(ndtr(-per_snr[snr]["K"] * eta))
            within = PA * trA + (1 - PA) * trB
            between = PA * bA + (1 - PA) * bB
            row = dict(snr=snr, eta=eta, eps=eps, sigma2=s2, x_norm2=q["x_norm2"], J_A=R["J_A"], J_B=R["J_B"],
                       Delta0=Delta0, Delta0_resid=q["JA_resid"] - q["JB_resid"], d_norm=q["d_norm"],
                       s_Delta=math.sqrt(2 * s2) * q["d_norm"], z_exact=-Delta0 / (math.sqrt(2 * s2) * q["d_norm"]),
                       P_A=PA, P_A_local=PAl, trC_A=trA, trC_B=trB, H_min_A=covA[1], H_min_B=covB[1],
                       fd_check_A=covA[2], fd_check_B=covB[2], uA1=uA[0], uA2=uA[1], uB1=uB[0], uB2=uB[1],
                       reference=("A" if eta < 0 else "B"), sqdist_A_to_W=bA, sqdist_B_to_W=bB,
                       MSE_pred=within + between, within=within, between=between,
                       MSE_pred_local=PAl * (trA + bA) + (1 - PAl) * (trB + bB))
            if eta == 0.0:
                dAB = float(np.sum((uA - uB) ** 2))
                row["MSE_pred_Aref"] = PA * trA + (1 - PA) * (trB + dAB)
            rows.append(row)

    # check 4: orientation
    ori = {}
    for snr in SNRS:
        ok = all((r["P_A"] > 0.5 if r["eta"] < 0 else r["P_A"] < 0.5) for r in rows if r["snr"] == snr and r["eta"] != 0)
        okl = all((r["P_A_local"] > 0.5 if r["eta"] < 0 else r["P_A_local"] < 0.5) for r in rows if r["snr"] == snr and r["eta"] != 0)
        mono = all(a["P_A"] >= b["P_A"] for a, b in zip([r for r in rows if r["snr"] == snr], [r for r in rows if r["snr"] == snr][1:]))
        ori[snr] = dict(exact=ok, local=okl, exact_monotone_nonincreasing=mono)
    checks["orientation"] = dict(g_positive=bool(g_env > 0), per_snr=ori,
                                 pass_=bool(g_env > 0 and all(v["exact"] and v["local"] for v in ori.values())))
    checks["d_norm_all_eta"] = {str(e): v for e, v in dnorms.items()}
    checks["no_fit"] = dict(fitted_parameters=[], pass_=True)

    # ---- secondary: exact-curve 10-90 width (deterministic continuation; no Monte Carlo)
    cache = {}

    def P_of(eta, snr):
        if eta not in cache:
            eps = EC * (1 + eta)
            qa, qb = continue_branches(eps)
            q = branch_quantities(eps, qa.pair(), qb.pair(), qa.J, qb.J)
            cache[eta] = (qa.J - qb.J, q["d_norm"], q["x"], qa.pair(), qb.pair())
        D0, dn, x, pa, pb = cache[eta]
        return p_exact(D0, sigma2_of(x, snr), dn)

    for snr in SNRS:
        res = {}
        for target, side in ((0.9, -1), (0.1, 1)):
            lo, hi = (0.0, side * 0.02) if side > 0 else (side * 0.02, 0.0)
            f = lambda e: P_of(e, snr) - target
            span = 0.02
            ok = False
            while span <= 0.6:
                a, b = (-span, 0.0) if side < 0 else (0.0, span)
                if f(a) * f(b) < 0:
                    ok = True
                    break
                span *= 1.5
            res[target] = brentq(f, a, b, xtol=1e-12) if ok else float("nan")
        per_snr[snr]["W_exact_curve"] = res[0.1] - res[0.9]
        per_snr[snr]["eta_P090"], per_snr[snr]["eta_P010"] = res[0.9], res[0.1]
    per_snr_ratio = dict(W30_over_W40_closed_form=per_snr[30]["W_closed_form"] / per_snr[40]["W_closed_form"],
                         W30_over_W40_exact_curve=per_snr[30]["W_exact_curve"] / per_snr[40]["W_exact_curve"],
                         W20_over_W30_closed_form=per_snr[20]["W_closed_form"] / per_snr[30]["W_closed_form"])
    eta50_pred = -(rc["J_A"] - rc["J_B"]) / (g_env * EC)

    det = dict(
        plan_seal=json.loads((ROOT / "paper/transition_theory/SEAL_1.json").read_text()),
        eps_c=EC, eps_c_mp_crossing_root=mp.nstr(e_star, 25), eta_of_mp_crossing_root=float(e_star / EC - 1),
        g=g_env, g_A=gA_env, g_B=gB_env, g_mp=mp.nstr(g_env_mp, 25), d_norm_c=dn_c, d_norm_c_mp=mp.nstr(dn_mp, 25),
        J_A_c=rc["J_A"], J_B_c=rc["J_B"], Delta0_c=rc["J_A"] - rc["J_B"], x_norm2_c=bq["x_norm2"],
        u_A_c=rc["A"], u_B_c=rc["B"], sqdist_uA_uB_c=float(np.sum((np.array(rc["A"]) - np.array(rc["B"])) ** 2)),
        eta50_pred_exact=eta50_pred, eta50_pred_local=0.0, per_snr={str(k): v for k, v in per_snr.items()},
        ratios=per_snr_ratio, Phi_inv_0p9=PHI_INV_09,
        environment=dict(python=sys.version.split()[0], numpy=np.__version__, scipy=scipy.__version__,
                         mpmath=mp.__version__, platform=platform.platform()),
        inputs={rel: h for rel, h in SPEC["deterministic_inputs"].items()}, seconds=time.time() - t0)
    checks["all_pass"] = bool(all(v.get("pass_", True) for v in checks.values() if isinstance(v, dict)))
    with (OUT / "predictions.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(dict.fromkeys(k for r in rows for k in r)))
        w.writeheader()
        w.writerows(rows)
    (OUT / "deterministic.json").write_text(json.dumps(det, indent=1) + "\n")
    (OUT / "checks.json").write_text(json.dumps(checks, indent=1, default=float) + "\n")
    print(json.dumps(dict(g=g_env, d_norm_c=dn_c, eps_c=EC, per_snr=det["per_snr"], ratios=per_snr_ratio,
                          checks_all_pass=checks["all_pass"]), indent=1))


if __name__ == "__main__":
    main()
