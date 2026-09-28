"""Pre-registered crossing analysis (CA) on top of advcore.

A configuration is a dict with keys
  cid, stage, N, z, b, phi, delta, amp_minus, amp_plus, window, scale ('abs'|'z2'|'z'), lo, hi
The scan parameter p is converted to epsilon = p * scale_factor.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.optimize import brentq

from advcore import (Objective, make_signal, weights, global_search, refine_from, merge_minima,
                     pair_distance, Minimum, level_params)

NPTS = 101
TOPK_STORE = 12


def eps_factor(cfg):
    sc = cfg.get("scale", "abs")
    if sc == "abs":
        return 1.0
    if sc == "z2":
        return cfg["z"] ** 2
    if sc == "z":
        return cfg["z"]
    raise ValueError(sc)


def build_obj(cfg, eps, noise=None):
    y = make_signal(cfg["N"], cfg["z"], cfg["b"], eps, cfg.get("phi", math.pi), cfg.get("delta", 0.0),
                    cfg.get("amp_minus", 1.0), cfg.get("amp_plus", 1.0))
    if noise is not None:
        y = y + noise
    return Objective(y, weights(cfg.get("window", "rect"), cfg["N"]))


def scan_grid(cfg, extended=False):
    lo, hi = cfg["lo"], cfg["hi"]
    if extended:
        lo, hi = hi, 3 * hi
    return np.linspace(lo, hi, NPTS)


def scan_task(args):
    """Worker: GS(level) at one scan point. Returns compact minima dicts."""
    cfg, k, p, level, extended = args
    eps = p * eps_factor(cfg)
    obj = build_obj(cfg, eps)
    mins = global_search(obj, level, z=cfg["z"], b=cfg["b"],
                         rng_key=(cfg["stage"], cfg["cid_int"], k + (1000 if extended else 0)),
                         **cfg.get("gs_over", {}))
    return dict(k=k, p=float(p), eps=float(eps), E=obj.E,
                minima=[q.as_dict() for q in mins[:TOPK_STORE]], n_minima=len(mins))


def _mk(d):
    return Minimum(d["m"], d["h"], d["J"], d["u1"], d["u2"], d["sep"], d["grad"],
                   (d["eig_mh_min"], d["eig_mh_max"]), (d["eig_u_min"], d["eig_u_max"]),
                   d["amp_max"], d["cls"])


def continuation_merge(cfg, pts):
    """Forward/backward continuation seeds from neighbouring scan points (P0 step 2)."""
    n = cfg["N"]
    lists = [[_mk(d) for d in pt["minima"]] for pt in pts]
    for direction in (1, -1):
        rng = range(1, len(pts)) if direction == 1 else range(len(pts) - 2, -1, -1)
        for k in rng:
            prev = lists[k - direction]
            seeds = [q.pair() for q in prev[:8]]
            if not seeds:
                continue
            obj = build_obj(cfg, pts[k]["eps"])
            ref = [q for q in refine_from(obj, seeds, tol_g=1e-10)
                   if q.eig_mh[0] >= -1e-9 * max(1.0, abs(q.J)) and q.grad <= 1e-7 * max(obj.E, 1.0)]
            lists[k] = merge_minima(n, [lists[k], ref])[:TOPK_STORE]
    return lists


class Branch:
    """Warm-started local continuation of one minimum as epsilon varies."""

    def __init__(self, cfg, eps0, pair):
        self.cfg = cfg
        self.cache = {float(eps0): tuple(pair)}
        self.fold = False

    def at(self, eps, tol_g=1e-11):
        keys = np.array(list(self.cache))
        e0 = float(keys[np.argmin(np.abs(keys - eps))])
        start = self.cache[e0]
        # sub-step if far
        nsub = max(1, int(math.ceil(abs(eps - e0) / (0.004 * max(1.0, self.cfg["hi"] * eps_factor(self.cfg))))))
        cur = start
        for j in range(1, nsub + 1):
            e = e0 + (eps - e0) * j / nsub
            obj = build_obj(self.cfg, e)
            q = refine_from(obj, [cur], tol_g=tol_g)[0]
            if pair_distance(q.pair(), cur, self.cfg["N"]) > 0.5:
                self.fold = True
            cur = q.pair()
        self.cache[float(eps)] = cur
        return q


def analyse_config(cfg, pts, validate_level="R"):
    """Crossing analysis for one configuration given its scan points (sorted by k)."""
    n = cfg["N"]
    pts = sorted(pts, key=lambda d: d["k"])
    lists = continuation_merge(cfg, pts)
    eps = np.array([pt["eps"] for pt in pts])
    best = [L[0] if L else None for L in lists]
    switches = []
    for k in range(len(pts) - 1):
        a, b = best[k], best[k + 1]
        if a is None or b is None:
            continue
        if pair_distance(a.pair(), b.pair(), n) > 1.0:
            switches.append(k)
    results = []
    for k in switches:
        results.append(solve_switch(cfg, eps[k], eps[k + 1], best[k], best[k + 1], validate_level, k))
    return lists, results


def solve_switch(cfg, e0, e1, P0, Q0, validate_level, k):
    n = cfg["N"]
    rec = dict(k=int(k), eps_lo=float(e0), eps_hi=float(e1), P_lo=P0.as_dict(), Q_hi=Q0.as_dict())
    bp = Branch(cfg, e0, P0.pair())
    bq = Branch(cfg, e1, Q0.pair())
    try:
        p_at_1 = bp.at(e1)
        q_at_0 = bq.at(e0)
    except Exception as exc:  # pragma: no cover
        rec.update(type="error", error=str(exc))
        return rec
    rec["P_fold"] = bool(bp.fold)
    rec["Q_fold"] = bool(bq.fold)
    if bp.fold or bq.fold:
        rec["type"] = "fold"
        return rec
    g0 = P0.J - q_at_0.J
    g1 = p_at_1.J - Q0.J
    rec["gap_lo"], rec["gap_hi"] = float(g0), float(g1)
    if not (g0 <= 0 <= g1 or g0 >= 0 >= g1):
        rec["type"] = "no_sign_change"
        return rec

    def gap(e):
        return bp.at(e).J - bq.at(e).J

    try:
        ec = brentq(gap, e0, e1, xtol=1e-14, rtol=1e-15, maxiter=200)
    except Exception as exc:
        rec.update(type="brent_failed", error=str(exc))
        return rec
    if bp.fold or bq.fold:
        rec["type"] = "fold"
        return rec
    P = bp.at(ec, tol_g=1e-12)
    Q = bq.at(ec, tol_g=1e-12)
    d = 1e-6 * max(1.0, abs(ec))
    slope = (gap(ec + d) - gap(ec - d)) / (2 * d)
    Jstar = 0.5 * (P.J + Q.J)
    obj = build_obj(cfg, ec)
    extra = [P.pair(), Q.pair()]
    vm = global_search(obj, validate_level, z=cfg["z"], b=cfg["b"], extra=extra,
                       rng_key=(cfg["stage"], cfg["cid_int"], 50000 + k))
    tol = 1e-9 * max(1.0, abs(Jstar))
    lower_other = [q for q in vm if q.J < Jstar - tol and pair_distance(q.pair(), P.pair(), n) > 0.05
                   and pair_distance(q.pair(), Q.pair(), n) > 0.05]
    foundP = any(pair_distance(q.pair(), P.pair(), n) < 1e-3 for q in vm)
    foundQ = any(pair_distance(q.pair(), Q.pair(), n) < 1e-3 for q in vm)
    others = [q for q in vm if pair_distance(q.pair(), P.pair(), n) > 0.05
              and pair_distance(q.pair(), Q.pair(), n) > 0.05]
    third = others[0] if others else None
    coal = [q for q in vm if q.cls == "coalescent"]
    third_margin = (third.J - Jstar) if third else float("inf")
    coal_margin = (min(q.J for q in coal) - Jstar) if coal else float("inf")
    dist = float(pair_distance(P.pair(), Q.pair(), n))
    c1 = dist >= 1.0
    regP = P.sep >= 1.0 and P.eig_u[0] > 1e-8 and P.amp_max < 1e3
    regQ = Q.sep >= 1.0 and Q.eig_u[0] > 1e-8 and Q.amp_max < 1e3
    c2 = regP and regQ
    c3_orig = abs(slope) >= 1e-3
    c4_orig = third_margin >= max(1e-4, 0.02 * abs(Jstar))
    # Amendment A1: scale-consistent thresholds (see AMENDMENTS.md)
    E = obj.E
    c3 = (abs(ec) * abs(slope) / max(abs(Jstar), 1e-300) >= 1e-3) and (abs(slope) >= 1e-10 * E / max(abs(ec), 1e-300))
    c4 = third_margin >= max(0.02 * abs(Jstar), 1e-10 * E)
    c5 = (not lower_other) and foundP and foundQ
    if c1 and c2 and c3_orig and c4_orig and c5:
        status_orig = "present"
    elif c1 and c2 and c3_orig:
        status_orig = "not_global" if lower_other else "inconclusive"
    else:
        status_orig = "not_qualifying"
    if c1 and c2 and c3 and c4 and c5:
        status = "present"
    elif c1 and c2 and c3:
        status = "not_global" if lower_other else "inconclusive"
    else:
        status = "not_qualifying"
    near_coal = (min(P.sep, Q.sep) < 1.0 or max(P.amp_max, Q.amp_max) >= 1e3
                 or coal_margin < 0.01 * abs(Jstar))
    rec.update(type="exchange", eps_c=float(ec), J_star=float(Jstar), J_P=float(P.J), J_Q=float(Q.J),
               P=P.as_dict(), Q=Q.as_dict(), slope=float(slope), branch_distance=dist,
               third=(third.as_dict() if third else None), third_margin=float(third_margin),
               third_margin_rel=float(third_margin / abs(Jstar)) if Jstar else float("inf"),
               coalescent_margin=float(coal_margin), n_validation_minima=len(vm),
               validation_lower_other=[q.as_dict() for q in lower_other[:3]],
               validation_found_P=bool(foundP), validation_found_Q=bool(foundQ),
               C1=bool(c1), C2=bool(c2), C3=bool(c3), C4=bool(c4), C5=bool(c5), status=status,
               C3_orig=bool(c3_orig), C4_orig=bool(c4_orig), status_orig=status_orig, E=float(E),
               transversality=float(abs(ec) * abs(slope) / max(abs(Jstar), 1e-300)),
               approaches_coalescence=bool(near_coal))
    return rec


def classify(results):
    """Configuration-level verdict from switch records (pre-registered)."""
    ex = [r for r in results if r.get("type") == "exchange"]
    if any(r["status"] == "present" for r in ex):
        return "present"
    if any(r["status"] in ("inconclusive",) for r in ex):
        return "inconclusive"
    return "absent"
