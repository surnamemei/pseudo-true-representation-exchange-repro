"""Stage 7: targeted numerical-invariance audit (search density, seeds, tolerance, precision, certificate path)."""
from __future__ import annotations

import ast
import csv
import glob
import json
import math
import random
from fractions import Fraction
from pathlib import Path

import numpy as np
from mpmath import mp, iv

from advrun import HERE, PROJECT, log, jdump, jload, status_update, run_configs, pmap
from advcore import Objective, make_signal, weights, global_search

SD = HERE / "stage7_numinv"
VARIANTS = {
    "base": {},
    "grid_coarse": dict(mm=512, mh=129),
    "grid_fine": dict(mm=2048, mh=513),
    "rand_8": dict(n_rand=8),
    "rand_128": dict(n_rand=128),
    "kgrid_16": dict(k_grid=16),
    "kgrid_200": dict(k_grid=200),
    "tol_1e-7": dict(tol_g=1e-7),
    "tol_1e-11": dict(tol_g=1e-11),
}


def base_points():
    pts = {}
    pts["baseline"] = dict(N=21, z=2.0, b=10.0, phi=math.pi, delta=0.0, amp_minus=1.0, amp_plus=1.0, window="rect",
                           scale="abs", lo=0.0, hi=1.0)
    pts["hann"] = dict(pts["baseline"], window="hann")
    pts["kappa0p2_z0p5"] = dict(pts["baseline"], z=0.5, delta=0.2 * 0.5, scale="z2")
    pts["fixeddelta0p2_z0p2"] = dict(pts["baseline"], z=0.2, delta=0.2 * math.pi, scale="z")
    # weakest relative third margin among present crossings of stages 1-4
    best = None
    for f in glob.glob(str(HERE / "stage[1-4]_*" / "analysis" / "*.json")):
        a = jload(f)
        for s in a["switches"]:
            if s.get("type") == "exchange" and s["status"] == "present":
                if best is None or s["third_margin_rel"] < best[0]:
                    best = (s["third_margin_rel"], a["config"], s["eps_c"], f)
    c = dict(best[1])
    for k in ("cid", "cid_int", "stage"):
        c.pop(k, None)
    pts["weakest_margin"] = c
    return pts, dict(third_margin_rel=best[0], eps_c=best[2], source=Path(best[3]).name)


def noisy_labels(args):
    variant, over = args
    from s6_noise import ETAS, EC, N, Z, BW, label
    roots = {float(k): v for k, v in jload(HERE / "stage6_noise" / "noiseless_roots.json").items()}
    ei = ETAS.index(0.0)
    x = make_signal(N, Z, BW, EC)
    sigma2 = float(np.vdot(x, x).real / N) * 10 ** (-30 / 10)
    rng = np.random.default_rng(np.random.SeedSequence([20260928, 6, 30, ei, 0]))
    out = []
    for t in range(50):
        w = math.sqrt(sigma2 / 2) * (rng.normal(size=N) + 1j * rng.normal(size=N))
        o = Objective(x + w)
        best = global_search(o, "P", blind=True, rng_key=(6, 30, ei, 0, t), **over)[0]
        out.append(dict(t=t, J=best.J, label=label(best.pair(), roots[0.0])[0], u1=best.u1, u2=best.u2))
    return variant, out


def mp_crossing_weighted(cfg, P, Q, eps_c, dps=40):
    """Independent mpmath polish of (P, Q, eps_c) for a weighted objective."""
    mp.dps = dps
    n = cfg["N"]
    w = weights(cfg["window"], n)
    ss = [mp.mpf(t) / n - mp.mpf(n - 1) / (2 * n) for t in range(n)]
    wv = [mp.mpf(repr(float(x))) for x in w]
    fac = {"abs": 1, "z2": cfg["z"] ** 2, "z": cfg["z"]}[cfg["scale"]]

    def rec(e):
        return [cfg["amp_minus"] * mp.expj(-(cfg["z"] * s + cfg["delta"])) + cfg["amp_plus"] * mp.expj(cfg["z"] * s + cfg["delta"])
                + e * mp.expj(cfg["phi"]) * mp.expj(cfg["b"] * s) for s in ss]

    def J(u1, u2, e):
        x = rec(e)
        v1 = [mp.expj(u1 * s) for s in ss]
        v2 = [mp.expj(u2 * s) for s in ss]
        G11 = mp.fsum(wv)
        G12 = mp.fsum(ww * mp.conj(a) * b for ww, a, b in zip(wv, v1, v2))
        h1 = mp.fsum(ww * mp.conj(a) * xx for ww, a, xx in zip(wv, v1, x))
        h2 = mp.fsum(ww * mp.conj(b) * xx for ww, b, xx in zip(wv, v2, x))
        det = G11 * G11 - abs(G12) ** 2
        B1 = (G11 * h1 - G12 * h2) / det
        B2 = (G11 * h2 - mp.conj(G12) * h1) / det
        return mp.fsum(ww * abs(xx - B1 * a - B2 * b) ** 2 for ww, xx, a, b in zip(wv, x, v1, v2))

    def F(a1, a2, b1, b2, e):
        return [mp.diff(lambda t: J(t, a2, e), a1), mp.diff(lambda t: J(a1, t, e), a2),
                mp.diff(lambda t: J(t, b2, e), b1), mp.diff(lambda t: J(b1, t, e), b2), J(a1, a2, e) - J(b1, b2, e)]

    x0 = [mp.mpf(P[0]), mp.mpf(P[1]), mp.mpf(Q[0]), mp.mpf(Q[1]), mp.mpf(eps_c)]
    sol = mp.findroot(F, x0, tol=mp.mpf(10) ** (-(dps - 10)), maxsteps=40)
    return float(sol[4]), mp.nstr(sol[4], 30)


def interval_path_checks():
    rnd = random.Random(20260928)
    iv.dps = 50
    res = dict(iv_arith=0, iv_arith_fail=0, iv_trans=0, iv_trans_fail=0, arb_arith=0, arb_arith_fail=0,
               arb_trans=0, arb_trans_fail=0)

    def frac_raw(t):
        sign, man, exp, _ = t
        v = Fraction(man) * Fraction(2) ** exp if exp >= 0 else Fraction(man, 2 ** (-exp))
        return -v if sign else v
    for _ in range(400):
        a = Fraction(rnd.randint(-10 ** 12, 10 ** 12), rnd.randint(1, 10 ** 9))
        b = Fraction(rnd.randint(-10 ** 12, 10 ** 12), rnd.randint(1, 10 ** 9)) or Fraction(1, 3)
        A = iv.mpf([str(a.numerator), str(a.numerator)]) / iv.mpf(a.denominator)
        B = iv.mpf([str(b.numerator), str(b.numerator)]) / iv.mpf(b.denominator)
        for exact, enc in ((a + b, A + B), (a - b, A - B), (a * b, A * B), (a / b, A / B)):
            res["iv_arith"] += 1
            lo_t, hi_t = enc._mpi_
            if not (frac_raw(lo_t) <= exact <= frac_raw(hi_t)):
                res["iv_arith_fail"] += 1
    mp.dps = 120
    for _ in range(200):
        x = mp.mpf(rnd.uniform(-60, 60))
        X = iv.mpf(x)
        for f, g in ((mp.sin, iv.sin), (mp.cos, iv.cos), (mp.exp, iv.exp)):
            if f is mp.exp and x > 40:
                continue
            ref = f(x)
            enc = g(X)
            res["iv_trans"] += 1
            if not (mp.mpf(enc.a) <= ref <= mp.mpf(enc.b)):
                res["iv_trans_fail"] += 1
    try:
        import sys
        sys.path.insert(0, str(PROJECT / "stage12_deps"))
        from flint import arb, ctx, fmpq
        ctx.prec = 200
        for _ in range(300):
            a = Fraction(rnd.randint(-10 ** 12, 10 ** 12), rnd.randint(1, 10 ** 9))
            b = Fraction(rnd.randint(1, 10 ** 12), rnd.randint(1, 10 ** 9))
            A = arb(fmpq(a.numerator, a.denominator))
            B = arb(fmpq(b.numerator, b.denominator))
            for exact, enc in ((a + b, A + B), (a * b, A * B), (a / b, A / B)):
                res["arb_arith"] += 1
                if not enc.overlaps(arb(fmpq(exact.numerator, exact.denominator))):
                    res["arb_arith_fail"] += 1
        mp.dps = 120
        for _ in range(200):
            xv = rnd.uniform(-60, 60)
            X = arb(repr(xv))
            for name in ("sin", "cos", "exp"):
                if name == "exp" and xv > 40:
                    continue
                enc = getattr(X, name)()
                ref = getattr(mp, name)(mp.mpf(repr(xv)))
                res["arb_trans"] += 1
                if not enc.overlaps(arb(mp.nstr(ref, 110))):
                    res["arb_trans_fail"] += 1
        res["arb_available"] = True
    except Exception as exc:
        res["arb_available"] = f"error: {exc!r}"
    mp.dps = 60
    return res


def import_audit():
    files = {
        "stage12_full_arb_replay.py": PROJECT / "stage12_full_arb_replay.py",
        "stage12_arb_backend.py": PROJECT / "stage12_arb_backend.py",
        "stage16_independent_arb.py": PROJECT / "stage16_independent_arb.py",
        "bcert_arb.py (this run)": HERE / "bcert_arb.py",
    }
    primary = {"three_to_two_tone_stage2_interval", "stage14_continuum_global_interval", "stage15_uniform_interval",
               "stage17_uniform_trial", "bcert_primary", "bcert_block", "three_to_two_tone_stage4_global_interval"}
    out = {}
    for name, p in files.items():
        if not p.exists():
            out[name] = "absent"
            continue
        tree = ast.parse(p.read_text(encoding="utf-8"))
        mods = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                mods |= {a.name.split(".")[0] for a in node.names}
            elif isinstance(node, ast.ImportFrom) and node.module:
                mods.add(node.module.split(".")[0])
        out[name] = dict(imports=sorted(mods), imports_primary=sorted(mods & primary), independent=not (mods & primary))
    return out


def reproducibility():
    """Same seeds -> bitwise identical global-search output (two fresh evaluations)."""
    o = Objective(make_signal(21, 2.0, 10.0, 0.2481906301722774))
    a = global_search(o, "P", z=2.0, b=10.0, rng_key=(7, 1, 2))
    b = global_search(Objective(make_signal(21, 2.0, 10.0, 0.2481906301722774)), "P", z=2.0, b=10.0, rng_key=(7, 1, 2))
    same = len(a) == len(b) and all(x.as_dict() == y.as_dict() for x, y in zip(a, b))
    return dict(global_search_bitwise_identical=bool(same), n_minima=len(a))


def run(ctx=None):
    SD.mkdir(parents=True, exist_ok=True)
    done = SD / "DONE.json"
    if done.exists():
        log("Stage 7 already complete", 7)
        return jload(done)
    status_update(7, status="running")
    pts, weakest = base_points()
    cfgs = []
    for pname, p in pts.items():
        for vname, over in VARIANTS.items():
            c = dict(p, cid=f"S7_{pname}_{vname}", stage=7, gs_over=over, point=pname, variant=vname)
            cfgs.append(c)
    res = run_configs(7, SD, cfgs, desc="S7 variants")
    table = []
    for c in cfgs:
        a = res[c["cid"]]
        ex = [s for s in a["switches"] if s.get("type") == "exchange"]
        pres = [s for s in ex if s["status"] == "present"]
        p = pres[0] if pres else (ex[0] if ex else None)
        table.append(dict(point=c["point"], variant=c["variant"], verdict=a["verdict"],
                          eps_c=p["eps_c"] if p else None, status=p["status"] if p else None,
                          third_margin_rel=p["third_margin_rel"] if p else None,
                          n_exchanges=len(ex)))
    changes = []
    for pname in pts:
        base = next(t for t in table if t["point"] == pname and t["variant"] == "base")
        for t in table:
            if t["point"] != pname or t["variant"] == "base":
                continue
            rel = (abs(t["eps_c"] - base["eps_c"]) / abs(base["eps_c"])) if (t["eps_c"] and base["eps_c"]) else None
            t["eps_c_rel_change"] = rel
            if t["verdict"] != base["verdict"] or t["status"] != base["status"] or (rel is not None and rel > 1e-6):
                changes.append(dict(point=pname, variant=t["variant"], base=base["verdict"], now=t["verdict"],
                                    rel_change=rel))
    with (SD / "variant_table.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(dict.fromkeys(k for t in table for k in t)))
        w.writeheader()
        w.writerows(table)
    # precision: mpmath polish of three crossings
    prec = {}
    for pname in ("baseline", "hann", "kappa0p2_z0p5"):
        a = res[f"S7_{pname}_base"]
        s = next((x for x in a["switches"] if x.get("type") == "exchange" and x["status"] == "present"), None)
        if s:
            try:
                f64 = s["eps_c"]
                mpv, mps = mp_crossing_weighted(pts[pname], (s["P"]["u1"], s["P"]["u2"]), (s["Q"]["u1"], s["Q"]["u2"]),
                                                s["eps_c"])
                prec[pname] = dict(float64=f64, mpmath40=mps, rel_diff=abs(f64 - mpv) / abs(mpv))
            except Exception as exc:
                prec[pname] = dict(error=repr(exc))
    # noisy labels under variants
    nres = pmap(noisy_labels, list(VARIANTS.items()), desc="S7 noisy", stage=7)
    base_n = dict(nres)["base"]
    noisy = {}
    for v, rows in nres:
        diff = sum(1 for x, y in zip(rows, base_n) if x["label"] != y["label"] or abs(x["J"] - y["J"]) > 1e-8)
        noisy[v] = dict(label_or_J_changes=diff, labels={k: sum(r["label"] == k for r in rows) for k in ("A", "B", "OUT")})
    out = dict(points=list(pts), weakest_margin_point=weakest, variant_changes=changes, n_variant_runs=len(cfgs),
               precision=prec, noisy=noisy, interval_path=interval_path_checks(), replay_independence=import_audit(),
               reproducibility=reproducibility())
    out["classification_changes"] = [c for c in changes if c["base"] != c["now"]]
    out["any_classification_change"] = bool(out["classification_changes"]) or any(v["label_or_J_changes"] for v in noisy.values())
    jdump(done, out)
    write_report(out, table)
    status_update(7, status="done", any_classification_change=out["any_classification_change"])
    log(f"Stage 7 any classification change: {out['any_classification_change']}", 7)
    return out


def write_report(out, table):
    L = ["# Numerical invariance report (Stage 7)", "",
         "Evidence level of every row: GLOBAL_NUMERICAL (search-protocol variations) unless stated.", "",
         f"Weakest relative third-margin qualifying crossing in Stages 1-4: {out['weakest_margin_point']}", "",
         "## Search-protocol variations (one factor at a time from production P)", "",
         "| point | variant | verdict | eps_c | rel. change vs base | status |", "|---|---|---|---|---|---|"]
    for t in table:
        L.append(f"| {t['point']} | {t['variant']} | {t['verdict']} | {t['eps_c']} | {t.get('eps_c_rel_change')} | {t['status']} |")
    L += ["", f"Classification changes: {out['classification_changes'] or 'none'}",
          f"Other flagged changes (eps_c rel > 1e-6 or status change): {out['variant_changes'] or 'none'}", "",
          "## Floating precision", "", json.dumps(out["precision"], indent=1), "",
          "## Noisy near-crossing records (50 records, eta=0, 30 dB)", "", json.dumps(out["noisy"], indent=1), "",
          "## Certificate arithmetic path", "", json.dumps(out["interval_path"], indent=1), "",
          "## Replay independence (static import audit)", "", json.dumps(out["replay_independence"], indent=1), "",
          "## Deterministic reproducibility", "", json.dumps(out["reproducibility"], indent=1), ""]
    (HERE / "numerical_invariance_report.md").write_text("\n".join(L), encoding="utf-8")
