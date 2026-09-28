"""Stage 3B two-level block certificate (uses the directed-interval kernels of bcert_primary).

Block = b in [b0, b0 + 0.01], split into 100 fine sub-tiles of width 1e-4.
  local : parametric Krawczyk per sub-tile (mean-value form in b) + convexity of each root box
  near  : per sub-tile, windows W_A, W_B (hull of the block's root boxes +- 0.05) minus that
          sub-tile's root box are excluded (cost / gradient / concave / monotone)
  far   : [-100,100] minus the windows, adaptive (v, b)-bisection with hull lambda and
          incumbent enclosures over the b sub-range
  plus coalescent and |v|>=100 tail margins per sub-tile, exact rational coverage checks.
"""
from __future__ import annotations

import math
import time
from decimal import Decimal, getcontext
from fractions import Fraction

from mpmath import iv, mp

from bcert_primary import (TileEval, candidate, krawczyk, lo, hi, excl0)

BLOCK_W = Decimal("0.01")
SUB_N = 100
NEAR_M = 0.1


def sub_edges(b0: str):
    getcontext().prec = 40
    B0 = Decimal(b0)
    return [str(B0 + BLOCK_W * Decimal(k) / Decimal(SUB_N)) for k in range(SUB_N + 1)]


def local_subtile(e0, e1):
    T = iv.mpf([e0, e1])
    bc = (float(e0) + float(e1)) / 2
    xc, dx = candidate(bc)
    half = (float(e1) - float(e0)) / 2
    for fac in (4.0, 6.0, 3.0):
        rad = [abs(dx[i]) * half * fac + 1e-10 * (1 + abs(xc[i])) for i in range(3)]
        loc, X, K = krawczyk(T, xc, rad, bc)
        if loc["ok"]:
            break
    loc["center"] = [repr(float(v)) for v in xc]
    loc["radius"] = [repr(float(v)) for v in rad]
    return loc, T, X, K


def predicate(te, U, a, c):
    z = te.calc(a, c)
    if z is None:
        return None
    R, Rv, Rvv = z
    if lo(R) > hi(U):
        return "cost"
    if excl0(Rv):
        return "gradient"
    if hi(Rvv) < 0:
        return "concave"
    if lo(Rvv) > 0:
        za, zc = te.calc(a, a), te.calc(c, c)
        if za and zc and ((lo(za[1]) > 0 and lo(zc[1]) > 0) or (hi(za[1]) < 0 and hi(zc[1]) < 0)):
            return "monotone"
    return None


def exclude_cells(te, U, cells, min_width=1e-10, budget=20000):
    leaves, un = [], []
    stack = list(reversed(cells))
    visits = 0
    while stack:
        a, c = stack.pop()
        visits += 1
        kind = predicate(te, U, a, c)
        if kind:
            leaves.append((a, c, kind))
            continue
        if c - a < min_width or visits > budget:
            un.append((a, c))
            if visits > budget:
                un.extend(stack)
                break
            continue
        m = (a + c) / 2
        stack.extend([(m, c), (a, m)])
    return leaves, un, visits


def chain_ok(segs, lo_v, hi_v):
    segs = sorted(segs)
    return (Fraction(segs[0][0]) == Fraction(lo_v) and Fraction(segs[-1][1]) == Fraction(hi_v)
            and all(Fraction(p[1]) == Fraction(q[0]) for p, q in zip(segs, segs[1:])))


def certify_block(b0: str, budget_far=2_500_000, verbose=False):
    t0 = time.time()
    edges = sub_edges(b0)
    subs = []
    for i in range(SUB_N):
        loc, T, X, K = local_subtile(edges[i], edges[i + 1])
        rec = dict(i=i, b=[edges[i], edges[i + 1]], local=loc)
        if loc["ok"]:
            L = iv.mpf([max(lo(K[2]), lo(X[2])), min(hi(K[2]), hi(X[2]))])
            boxA = (math.nextafter(float(lo(X[0])), -math.inf), math.nextafter(float(hi(X[0])), math.inf))
            boxB = (math.nextafter(float(lo(X[1])), -math.inf), math.nextafter(float(hi(X[1])), math.inf))
            te = TileEval(T, L)
            zA, zB = te.calc(*boxA), te.calc(*boxB)
            # incumbent: R at one fixed point bounds the A-minimum from above (as in the frozen Stage-14 code)
            vc = float(loc["center"][0])
            zc = te.calc(vc, vc)
            rec.update(L=L, boxA=boxA, boxB=boxB, U=zc[0] if zc else None, te=te, U_point=vc,
                       rootA_convex=bool(zA and lo(zA[2]) > 0), rootB_convex=bool(zB and lo(zB[2]) > 0))
        subs.append(rec)
    loc_ok = all(r["local"]["ok"] and r.get("rootA_convex") and r.get("rootB_convex") for r in subs)
    if verbose:
        print(f"  local done {time.time()-t0:.1f}s ok={loc_ok}", flush=True)
    out = dict(block=[edges[0], edges[-1]], n_sub=SUB_N, local_all_ok=bool(loc_ok),
               local_fail=[r["b"] for r in subs if not (r["local"]["ok"] and r.get("rootA_convex") and r.get("rootB_convex"))][:10])
    if not loc_ok:
        out.update(ok=False, reason="local", seconds=time.time() - t0)
        return out, None
    WA = (min(r["boxA"][0] for r in subs) - NEAR_M, max(r["boxA"][1] for r in subs) + NEAR_M)
    WB = (min(r["boxB"][0] for r in subs) - NEAR_M, max(r["boxB"][1] for r in subs) + NEAR_M)
    out.update(window_A=WA, window_B=WB)
    coal_min = tail_min = None
    for r in subs:
        cm = lo(r["te"].coalescent()) - hi(r["U"])
        tm = lo(r["te"].tail()) - hi(r["U"])
        coal_min = cm if coal_min is None else min(coal_min, cm)
        tail_min = tm if tail_min is None else min(tail_min, tm)
    out.update(coalescent_margin_min=mp.nstr(coal_min, 10), tail_margin_min=mp.nstr(tail_min, 10))
    cache = {}

    def range_eval(i0, i1):
        if (i0, i1) not in cache:
            T = iv.mpf([edges[i0], edges[i1]])
            L = iv.mpf([min(lo(subs[k]["L"]) for k in range(i0, i1)), max(hi(subs[k]["L"]) for k in range(i0, i1))])
            U = iv.mpf([min(lo(subs[k]["U"]) for k in range(i0, i1)), max(hi(subs[k]["U"]) for k in range(i0, i1))])
            cache[(i0, i1)] = (TileEval(T, L), U)
        return cache[(i0, i1)]

    def in_window(a, c):
        return (a >= WA[0] and c <= WA[1]) or (a >= WB[0] and c <= WB[1])

    cuts = sorted(set([-100.0, 100.0] + [k / 4 for k in range(-400, 401)] + [WA[0], WA[1], WB[0], WB[1]]))
    cells = [(a, c) for a, c in zip(cuts, cuts[1:])]
    stack = [(a, c, 0, SUB_N) for a, c in reversed(cells)]
    leaves, far_un, vis = [], [], 0
    tl = time.time()
    while stack:
        a, c, i0, i1 = stack.pop()
        vis += 1
        if verbose and vis % 20000 == 0:
            import collections
            print(f"  visits {vis} leaves {len(leaves)} stack {len(stack)} ranges {len(cache)} t={time.time()-tl:.0f}s "
                  f"last=({a:.4f},{c:.4f},{i0},{i1}) kinds={collections.Counter(k for *_, k in leaves)}", flush=True)
        te, U = range_eval(i0, i1)
        if in_window(a, c):
            z = te.calc(a, c)
            kind = "convex_window" if (z is not None and lo(z[2]) > 0) else None
        else:
            kind = predicate(te, U, a, c)
        if kind:
            leaves.append((a, c, i0, i1, kind))
            continue
        if vis > budget_far:
            far_un.append((a, c, i0, i1))
            far_un.extend(stack)
            break
        if (c - a > 2.0 ** -8 or i1 - i0 == 1) and c - a >= 1e-9:
            m = (a + c) / 2
            stack.extend([(m, c, i0, i1), (a, m, i0, i1)])
        elif i1 - i0 > 1:
            m = (i0 + i1) // 2
            stack.extend([(a, c, m, i1), (a, c, i0, m)])
        else:
            far_un.append((a, c, i0, i1))
    # root boxes must lie strictly inside their windows (uniqueness by strict convexity on the window)
    boxes_inside = all(WA[0] < r["boxA"][0] and r["boxA"][1] < WA[1] and WB[0] < r["boxB"][0] and r["boxB"][1] < WB[1]
                       for r in subs)
    cov_ok = not far_un and all(
        chain_ok([(a, c) for a, c, i0, i1, _ in leaves if i0 <= k < i1], -100.0, 100.0) for k in range(SUB_N))
    kinds = {}
    for *_, k in leaves:
        kinds[k] = kinds.get(k, 0) + 1
    out.update(visits=vis, leaves=len(leaves), unresolved=len(far_un), kinds=kinds, boxes_inside_windows=bool(boxes_inside),
               coverage_exact=bool(cov_ok), b_ranges_used=len(cache))
    out["ok"] = bool(loc_ok and boxes_inside and not far_un and cov_ok and coal_min > 0 and tail_min > 0)
    out["seconds"] = time.time() - t0
    m = [r["local"]["margins"] for r in subs]
    out["min_margins"] = dict(
        Rvv_A=min(float(x["Rvv_A"]) for x in m), Rvv_B=min(float(x["Rvv_B"]) for x in m),
        slope=min(float(x["slope"]) for x in m), beta_A=min(float(x["beta_A"]) for x in m),
        beta_B=min(float(x["beta_B"]) for x in m),
        gram=min(min(float(x["gram_A"]), float(x["gram_B"])) for x in m),
        contraction_max=max(r["local"]["contraction"] for r in subs),
        coalescent=float(coal_min), tail=float(tail_min),
        lambda_range=[float(min(lo(r["L"]) for r in subs)), float(max(hi(r["L"]) for r in subs))])
    replay = dict(edges=edges,
                  subs=[dict(i=r["i"], b=r["b"], center=r["local"]["center"], radius=r["local"]["radius"],
                             boxA=r["boxA"], boxB=r["boxB"]) for r in subs],
                  window_A=WA, window_B=WB, leaves=leaves)
    return out, replay
