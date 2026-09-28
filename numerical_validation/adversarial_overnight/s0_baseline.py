"""Stage 0: numerical baseline reproduction (certificate replays run via tools/replay_sandbox.py)."""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from mpmath import mp

from advrun import HERE, PROJECT, log, jdump, jload, status_update, run_configs
from advcore import Objective, make_signal, global_search, refine_from, pair_distance

SD = HERE / "stage0_baseline"
EC_FROZEN = 0.2481906301722774
REF = json.loads((PROJECT / "stage2_reference.json").read_text())


def mp_crossing(dps=50):
    mp.dps = dps
    N = 21
    s = [mp.mpf(t) / N for t in range(-10, 11)]

    def rec(eps):
        return [2 * mp.cos(2 * x) - eps * mp.expj(10 * x) for x in s]

    def J_and_grad(u1, u2, eps):
        x = rec(eps)
        v1 = [mp.expj(u1 * t) for t in s]
        v2 = [mp.expj(u2 * t) for t in s]
        h1 = mp.fsum(mp.conj(a) * b for a, b in zip(v1, x))
        h2 = mp.fsum(mp.conj(a) * b for a, b in zip(v2, x))
        g = mp.fsum(mp.conj(a) * b for a, b in zip(v1, v2))
        det = N * N - abs(g) ** 2
        B1 = (N * h1 - g * h2) / det
        B2 = (N * h2 - mp.conj(g) * h1) / det
        r = [xx - B1 * a - B2 * b for xx, a, b in zip(x, v1, v2)]
        J = mp.fsum(abs(q) ** 2 for q in r)
        g1 = -2 * mp.re(mp.fsum(mp.conj(q) * 1j * t * a * B1 for q, t, a in zip(r, s, v1)))
        g2 = -2 * mp.re(mp.fsum(mp.conj(q) * 1j * t * b * B2 for q, t, b in zip(r, s, v2)))
        return J, g1, g2

    def F(a1, a2, b1, b2, e):
        Ja, ga1, ga2 = J_and_grad(a1, a2, e)
        Jb, gb1, gb2 = J_and_grad(b1, b2, e)
        return [ga1, ga2, gb1, gb2, Ja - Jb]

    x0 = [mp.mpf(REF["A_at_crossing"][0][:25]), mp.mpf(REF["A_at_crossing"][1][:25]),
          mp.mpf(REF["B_at_crossing"][0][:25]), mp.mpf(REF["B_at_crossing"][1][:25]),
          mp.mpf(REF["epsilon_cross"][:25])]
    # perturb the start so the solve is not a copy of the reference
    x0 = [v * (1 + mp.mpf("1e-9")) for v in x0]
    sol = mp.findroot(F, x0, tol=mp.mpf(10) ** (-(dps - 6)), maxsteps=60)
    a1, a2, b1, b2, e = [sol[i] for i in range(5)]
    Ja = J_and_grad(a1, a2, e)[0]
    ref = [mp.mpf(REF["A_at_crossing"][0]), mp.mpf(REF["A_at_crossing"][1]),
           mp.mpf(REF["B_at_crossing"][0]), mp.mpf(REF["B_at_crossing"][1]), mp.mpf(REF["epsilon_cross"])]
    diffs = [abs(u - v) for u, v in zip([a1, a2, b1, b2, e], ref)]
    Jref = mp.mpf(REF["records"][1]["J_A"])
    return dict(dps=dps, eps_c=mp.nstr(e, 45), A=[mp.nstr(a1, 40), mp.nstr(a2, 40)],
                B=[mp.nstr(b1, 40), mp.nstr(b2, 40)], J_star=mp.nstr(Ja, 40),
                max_abs_diff_vs_110digit_reference=mp.nstr(max(diffs), 5),
                eps_c_abs_diff=mp.nstr(diffs[4], 5), J_star_abs_diff=mp.nstr(abs(Ja - Jref), 5),
                pass_B02=bool(max(diffs) < mp.mpf("1e-30") and abs(Ja - Jref) < mp.mpf("1e-30")))


def run(ctx=None):
    SD.mkdir(parents=True, exist_ok=True)
    out = SD / "baseline_numeric.json"
    if out.exists():
        log("Stage 0 numeric already complete", 0)
        return jload(out)
    status_update(0, status="running_numeric")
    cfg = dict(cid="S0_baseline_rect_N21_z2_b10", stage=0, N=21, z=2.0, b=10.0, phi=math.pi, delta=0.0,
               amp_minus=1.0, amp_plus=1.0, window="rect", scale="abs", lo=0.0, hi=1.0)
    res = run_configs(0, SD, [cfg], desc="S0 baseline CA")[cfg["cid"]]
    ex = [r for r in res["switches"] if r.get("type") == "exchange"]
    near = min(ex, key=lambda r: abs(r["eps_c"] - EC_FROZEN)) if ex else None
    rep = dict(ca_verdict=res["verdict"], ca_switches=len(res["switches"]))
    if near:
        rep.update(eps_c_float64=near["eps_c"], eps_c_abs_diff=abs(near["eps_c"] - EC_FROZEN),
                   pass_B01=abs(near["eps_c"] - EC_FROZEN) <= 1e-10, J_star=near["J_star"],
                   A=near["P"], B=near["Q"], third_margin=near["third_margin"], slope=near["slope"],
                   status=near["status"])
    # B0.3/B0.4 at the frozen crossing
    obj = Objective(make_signal(21, 2.0, 10.0, EC_FROZEN))
    A0 = (float(REF["A_at_crossing"][0]), float(REF["A_at_crossing"][1]))
    B0 = (float(REF["B_at_crossing"][0]), float(REF["B_at_crossing"][1]))
    fa, fb = refine_from(obj, [A0, B0], tol_g=1e-12)
    mins = global_search(obj, "R", z=2.0, b=10.0, extra=[A0, B0], rng_key=(0, 7, 7))
    two = mins[:2]
    is_ab = sorted([min(pair_distance(q.pair(), A0, 21), pair_distance(q.pair(), B0, 21)) for q in two])
    others = [q for q in mins if pair_distance(q.pair(), A0, 21) > 0.05 and pair_distance(q.pair(), B0, 21) > 0.05]
    third_margin = others[0].J - 0.5 * (fa.J + fb.J)
    rep.update(hessian_eig_u_A=fa.eig_u, hessian_eig_u_B=fb.eig_u,
               pass_B03=fa.eig_u[0] > 0 and fb.eig_u[0] > 0,
               GS_R_two_lowest_are_AB=bool(is_ab[1] < 1e-6), third_minimum=others[0].as_dict(),
               third_margin_at_frozen_ec=third_margin,
               pass_B04=bool(is_ab[1] < 1e-6 and abs(third_margin - 0.88296) < 1e-5))
    # B0.5 ordering
    order = []
    for e in (EC_FROZEN - 1e-3, EC_FROZEN + 1e-3, 0.213190630172278, 0.283190630172277):
        o = Objective(make_signal(21, 2.0, 10.0, e))
        a, b = refine_from(o, [fa.pair(), fb.pair()], tol_g=1e-12)
        order.append(dict(eps=e, J_A=a.J, J_B=b.J, A=a.pair(), B=b.pair(), A_below_B=a.J < b.J))
    rep["ordering"] = order
    rep["pass_B05"] = bool(order[0]["A_below_B"] and not order[1]["A_below_B"]
                           and order[2]["A_below_B"] and not order[3]["A_below_B"])
    log("mpmath 50-digit joint crossing solve ...", 0)
    rep["mp50"] = mp_crossing(50)
    rep["pass_B02"] = rep["mp50"]["pass_B02"]
    rep["pass_all_numeric"] = all(rep.get(k) for k in ("pass_B01", "pass_B02", "pass_B03", "pass_B04", "pass_B05"))
    rep["evidence_level"] = "GLOBAL_NUMERICAL (reproduction of an interval-certified result)"
    jdump(out, rep)
    status_update(0, status="numeric_done", numeric_pass=rep["pass_all_numeric"])
    log(f"Stage 0 numeric pass={rep['pass_all_numeric']} eps_c={rep.get('eps_c_float64')}", 0)
    return rep
