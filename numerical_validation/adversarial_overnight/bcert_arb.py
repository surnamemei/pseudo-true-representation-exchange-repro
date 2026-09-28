"""Stage 3B independent replay with python-flint / Arb ball arithmetic.

Does not import bcert_primary / bcert_block. Inputs are only the *geometry* written by the
primary run (b edges, Krawczyk centres/radii, root boxes, windows, leaf cells and their
claimed predicate). Every enclosure, preconditioner, lambda enclosure, incumbent and margin
is recomputed here from scratch in Arb; stored numbers are never used as validity inputs.
"""
from __future__ import annotations

import math
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "stage12_deps"))
from flint import arb, ctx, fmpq  # noqa: E402

ctx.prec = 256
MU2 = arb(1) / 12
MU4 = arb(1) / 80
VV2 = MU4 - MU2 * MU2
TWO = arb(2)
PI = arb.pi()


def ball(a: float, c: float):
    A, C = arb(a), arb(c)
    return (A + C) / 2 + arb(0, 1) * ((C - A) / 2)   # exact endpoints, outward


def dball(d0: str, d1: str):
    A = arb(fmpq(*Fraction(d0).as_integer_ratio()))
    C = arb(fmpq(*Fraction(d1).as_integer_ratio()))
    return (A + C) / 2 + arb(0, 1) * ((C - A) / 2)


def hull(xs):
    lo_ = min(x.lower() for x in xs)
    hi_ = max(x.upper() for x in xs)
    return (lo_ + hi_) / 2 + arb(0, 1) * ((hi_ - lo_) / 2)


def neg_def(x):
    return x < 0  # arb comparison is True only if certainly


def pos_def(x):
    return x > 0


def excludes0(x):
    return x > 0 or x < 0


def maxabs(x):
    return max(abs(x.lower()), abs(x.upper()))


# ----------------------------------------------------------------------------- kernel
SER = 40
_fact = [math.factorial(k) for k in range(2 * SER + 6)]


def Dk(x, k=0):
    if maxabs(x) <= 1:
        out = arb(0)
        for m in range((k + 1) // 2, SER):
            out += (-1) ** m * x ** (2 * m - k) / (arb(2 * m + 1) * _fact[2 * m - k] * TWO ** (2 * m))
        return out + arb(0, arb("1e-45"))
    if x.contains(0):
        return None
    out = 2 * (x / 2).sin() / x
    for j in range(1, k + 1):
        out = (2 * (x / 2 + j * PI / 2).sin() / TWO ** j - j * out) / x
    return out


# ----------------------------------------------------------------------------- jets in (v, lam, b)
class Jet:
    __slots__ = ("v", "g", "H")

    def __init__(self, v, g, H):
        self.v, self.g, self.H = v, g, H

    @staticmethod
    def const(x):
        z = arb(0)
        return Jet(arb(x) if not isinstance(x, arb) else x, [z, z, z], [[z] * 3 for _ in range(3)])

    @staticmethod
    def var(x, k):
        j = Jet.const(x)
        j.g = [arb(int(i == k)) for i in range(3)]
        return j

    def _c(self, o):
        return o if isinstance(o, Jet) else Jet.const(o)

    def __add__(a, b):
        b = a._c(b)
        return Jet(a.v + b.v, [x + y for x, y in zip(a.g, b.g)], [[a.H[i][j] + b.H[i][j] for j in range(3)] for i in range(3)])
    __radd__ = __add__

    def __neg__(a):
        return Jet(-a.v, [-x for x in a.g], [[-x for x in r] for r in a.H])

    def __sub__(a, b):
        return a + (-a._c(b))

    def __rsub__(a, b):
        return (-a) + b

    def __mul__(a, b):
        b = a._c(b)
        return Jet(a.v * b.v, [a.g[i] * b.v + a.v * b.g[i] for i in range(3)],
                   [[a.H[i][j] * b.v + a.v * b.H[i][j] + a.g[i] * b.g[j] + b.g[i] * a.g[j] for j in range(3)] for i in range(3)])
    __rmul__ = __mul__

    def inv(a):
        r = 1 / a.v
        return Jet(r, [-x * r * r for x in a.g], [[2 * a.g[i] * a.g[j] * r ** 3 - a.H[i][j] * r * r for j in range(3)] for i in range(3)])

    def __truediv__(a, b):
        return a * a._c(b).inv()


def Dj(x: Jet, k=0):
    d0, d1, d2 = Dk(x.v, k), Dk(x.v, k + 1), Dk(x.v, k + 2)
    return Jet(d0, [d1 * g for g in x.g], [[d1 * x.H[i][j] + d2 * x.g[i] * x.g[j] for j in range(3)] for i in range(3)])


def Rjet(v, lam, b):
    v, lam, b = Jet.var(v, 0), Jet.var(lam, 1), Jet.var(b, 2)
    d, dp, d2 = Dj(v), Dj(v, 1), Dj(v, 2)
    q0 = Jet.const(-MU2) - lam * Dj(b)
    q1 = 2 * lam * Dj(b, 1)
    A = 1 - d * d - dp * dp / MU2
    Bn = d2 - lam * Dj(b - v) - d * q0 + dp * q1 / (2 * MU2)
    C = Jet.const(MU4) - 2 * lam * Dj(b, 2) + lam * lam - q0 * q0 - q1 * q1 / (4 * MU2)
    return C - Bn * Bn / A, A, Bn


def arb_krawczyk(T, bc, center, radius):
    xc = [arb(float(c)) for c in center]
    rad = [arb(float(r)) for r in radius]
    X = [xc[i] + arb(0, 1) * rad[i] for i in range(3)]
    Bc = arb(bc)

    def FJ(xA, xB, lam, b):
        RA, AA, BA = Rjet(xA, lam, b)
        RB, AB, BB = Rjet(xB, lam, b)
        F = [RA.g[0], RB.g[0], RA.v - RB.v]
        Jx = [[RA.H[0][0], arb(0), RA.H[0][1]], [arb(0), RB.H[0][0], RB.H[0][1]], [RA.g[0], -RB.g[0], RA.g[1] - RB.g[1]]]
        Fb = [RA.H[0][2], RB.H[0][2], RA.g[2] - RB.g[2]]
        return F, Jx, Fb, (AA, BA, AB, BB)
    F0, J0, _, _ = FJ(xc[0], xc[1], xc[2], Bc)
    _, _, Fb, _ = FJ(xc[0], xc[1], xc[2], T)
    _, JX, _, (AA, BA, AB, BB) = FJ(X[0], X[1], X[2], T)
    mid = np.array([[float(J0[i][j].mid()) for j in range(3)] for i in range(3)])
    Ci = np.linalg.inv(mid)
    C = [[arb(float(Ci[i, j])) for j in range(3)] for i in range(3)]
    Fx = [F0[i] + Fb[i] * (T - Bc) for i in range(3)]
    Y = [X[i] - xc[i] for i in range(3)]
    K = []
    for i in range(3):
        acc = xc[i]
        for k in range(3):
            acc -= C[i][k] * Fx[k]
        for j in range(3):
            tij = arb(int(i == j)) - sum((C[i][k] * JX[k][j] for k in range(3)), arb(0))
            acc += tij * Y[j]
        K.append(acc)
    inside = all(K[i].lower() > X[i].lower() and K[i].upper() < X[i].upper() for i in range(3))
    marg = dict(Rvv_A=JX[0][0], Rvv_B=JX[1][1], slope=JX[2][2], gram_A=AA.v, gram_B=AB.v, beta_A=BA.v / AA.v,
                beta_B=-(BB.v / AB.v))
    ok = inside and all(pos_def(v) for v in marg.values())
    L = K[2].intersection(X[2]) if inside else X[2]
    return ok, L, {k: float(v.lower()) for k, v in marg.items()}, K


class Eval:
    def __init__(self, T, L):
        self.T, self.L = T, L
        self.q0 = -MU2 - L * Dk(T)
        self.q1 = 2 * L * Dk(T, 1)
        self.C = MU4 - 2 * L * Dk(T, 2) + L * L - self.q0 * self.q0 - self.q1 * self.q1 / (4 * MU2)
        K = 46
        dc = [arb(0)] * (K + 3)
        for m in range((K + 3) // 2):
            k = 2 * m
            if k < len(dc):
                dc[k] = arb((-1) ** m) / (_fact[k + 1] * TWO ** k)
        dp = [(k + 1) * dc[k + 1] for k in range(K + 1)]
        db = [(-1) ** k * Dk(T, k) / _fact[k] for k in range(K + 1)]
        ac, bc = [], []
        for k in range(K + 1):
            ac.append(arb(int(k == 0)) - sum((dc[j] * dc[k - j] + dp[j] * dp[k - j] / MU2 for j in range(k + 1)), arb(0)))
            bc.append((k + 2) * (k + 1) * dc[k + 2] - L * db[k] - dc[k] * self.q0 + dp[k] * self.q1 / (2 * MU2))
        self.ac, self.bc = ac[4:45], bc[2:43]

    @staticmethod
    def poly(coeff, x, r):
        out = arb(0)
        for k in reversed(range(r, len(coeff))):
            out = out * x + coeff[k] * (_fact[k] // _fact[k - r])
        return out + arb(0, arb("1e-30"))

    def calc(self, x):
        L, T = self.L, self.T
        if maxabs(x) <= 1:
            g, g1, g2 = [self.poly(self.ac, x, j) for j in range(3)]
            h, h1, h2 = [self.poly(self.bc, x, j) for j in range(3)]
        else:
            if x.contains(0):
                return None
            d, d1, d2, d3, d4 = [Dk(x, j) for j in range(5)]
            g = 1 - d * d - d1 * d1 / MU2
            g1 = -2 * d * d1 - 2 * d1 * d2 / MU2
            g2 = -2 * (d1 * d1 + d * d2) - 2 * (d2 * d2 + d1 * d3) / MU2
            w = T - x
            e0, e1, e2 = Dk(w), Dk(w, 1), Dk(w, 2)
            if e0 is None:
                return None
            h = d2 - L * e0 - d * self.q0 + d1 * self.q1 / (2 * MU2)
            h1 = d3 + L * e1 - d1 * self.q0 + d2 * self.q1 / (2 * MU2)
            h2 = d4 - L * e2 - d2 * self.q0 + d3 * self.q1 / (2 * MU2)
        if not pos_def(g):
            return None
        R = self.C - h * h / g
        Rv = -2 * h * h1 / g + h * h * g1 / (g * g)
        Rvv = -2 * (h1 * h1 + h * h2) / g + 4 * h * h1 * g1 / (g * g) + h * h * g2 / (g * g) - 2 * h * h * g1 * g1 / (g ** 3)
        return R, Rv, Rvv

    def coalescent(self):
        q2 = -MU4 + self.L * Dk(self.T, 2) - MU2 * self.q0
        return self.C - q2 * q2 / VV2

    def tail(self):
        V = arb(100)
        dd = 2 / V
        dpb = (1 + dd) / V
        d2b = 1 / (2 * V) + 2 * dpb / V
        dbw = 2 / (V - self.T.upper())
        AA = 1 - dd * dd - dpb * dpb / MU2
        BB = d2b + maxabs(self.L) * dbw + dd * maxabs(self.q0) + dpb * maxabs(self.q1) / (2 * MU2)
        return self.C - BB * BB / AA


def check_leaf(ev, U, a, c, kind):
    z = ev.calc(ball(a, c))
    if z is None:
        return False
    R, Rv, Rvv = z
    if kind == "cost":
        return R.lower() > U.upper()
    if kind == "gradient":
        return excludes0(Rv)
    if kind == "concave":
        return neg_def(Rvv)
    if kind in ("convex_window", "root_convex"):
        return pos_def(Rvv)
    if kind == "monotone":
        za, zc = ev.calc(arb(a)), ev.calc(arb(c))
        if not (za and zc and pos_def(Rvv)):
            return False
        return (pos_def(za[1]) and pos_def(zc[1])) or (neg_def(za[1]) and neg_def(zc[1]))
    return False


def chain_ok(segs, lo_v, hi_v):
    segs = sorted(segs)
    return (Fraction(segs[0][0]) == Fraction(lo_v) and Fraction(segs[-1][1]) == Fraction(hi_v)
            and all(Fraction(p[1]) == Fraction(q[0]) for p, q in zip(segs, segs[1:])))


def replay_block(rep):
    edges = rep["edges"]
    n = len(edges) - 1
    subs = []
    fails = []
    for s in rep["subs"]:
        i = s["i"]
        T = dball(edges[i], edges[i + 1])
        bc = (float(edges[i]) + float(edges[i + 1])) / 2
        ok, L, marg, K = arb_krawczyk(T, bc, s["center"], s["radius"])
        ev = Eval(T, L)
        vc = float(s["center"][0])
        zc = ev.calc(arb(vc))
        U = zc[0] if zc else None
        # every root of the joint system in X lies in K (Krawczyk); require K's v-projections inside the
        # primary's root boxes and the root boxes strictly inside the convexity windows
        boxA, boxB = s["boxA"], s["boxB"]
        box_ok = bool(ok and arb(boxA[0]) <= K[0].lower() and K[0].upper() <= arb(boxA[1]) and arb(boxB[0]) <= K[1].lower()
                      and K[1].upper() <= arb(boxB[1]) and rep["window_A"][0] < boxA[0] and boxA[1] < rep["window_A"][1]
                      and rep["window_B"][0] < boxB[0] and boxB[1] < rep["window_B"][1])
        cz = ev.coalescent()
        tz = ev.tail()
        cm = (cz - U).lower() if U is not None else None
        tm = (tz - U).lower() if U is not None else None
        good = bool(ok and U is not None and box_ok and cm > 0 and tm > 0)
        if not good:
            fails.append(dict(i=i, krawczyk=ok, box_ok=bool(box_ok)))
        subs.append(dict(T=T, L=L, U=U, marg=marg, coal=float(cm) if cm is not None else None,
                         tail=float(tm) if tm is not None else None))
    cache = {}

    def range_eval(i0, i1):
        if (i0, i1) not in cache:
            T = dball(edges[i0], edges[i1])
            L = hull([subs[k]["L"] for k in range(i0, i1)])
            U = hull([subs[k]["U"] for k in range(i0, i1)])
            cache[(i0, i1)] = (Eval(T, L), U)
        return cache[(i0, i1)]
    bad = 0
    refined = 0
    WA, WB = rep["window_A"], rep["window_B"]
    FAR = ("cost", "gradient", "concave", "monotone")

    def verify(a, c, i0, i1, inw, depth):
        nonlocal refined
        ev, U = range_eval(i0, i1)
        kinds = ("convex_window",) if inw else FAR
        if any(check_leaf(ev, U, a, c, k) for k in kinds):
            return True
        if depth >= 8:
            return False
        refined += 1
        if i1 - i0 > 1 and depth % 2 == 1:
            m = (i0 + i1) // 2
            return verify(a, c, i0, m, inw, depth + 1) and verify(a, c, m, i1, inw, depth + 1)
        m = (a + c) / 2
        return verify(a, m, i0, i1, inw, depth + 1) and verify(m, c, i0, i1, inw, depth + 1)

    for a, c, i0, i1, kind in rep["leaves"]:
        inw = (a >= WA[0] and c <= WA[1]) or (a >= WB[0] and c <= WB[1])
        if inw != (kind == "convex_window"):
            bad += 1
            continue
        ev, U = range_eval(i0, i1)
        if check_leaf(ev, U, a, c, kind):
            continue
        if not verify(a, c, i0, i1, inw, 0):
            bad += 1
    cov = all(chain_ok([(a, c) for a, c, i0, i1, _ in rep["leaves"] if i0 <= k < i1], -100.0, 100.0) for k in range(n))
    bcov = all(Fraction(edges[k + 1]) == Fraction(edges[k]) + (Fraction(edges[-1]) - Fraction(edges[0])) / n for k in range(n))
    out = dict(block=[edges[0], edges[-1]], subtiles=n, subtile_failures=fails[:10], n_subtile_failures=len(fails),
               leaves=len(rep["leaves"]), leaf_failures=bad, replay_refinements=refined, coverage_exact=bool(cov),
               b_edges_exact=bool(bcov),
               min_coalescent=min(s["coal"] for s in subs if s["coal"] is not None),
               min_tail=min(s["tail"] for s in subs if s["tail"] is not None),
               min_margins={k: min(s["marg"][k] for s in subs) for k in subs[0]["marg"]})
    out["ok"] = bool(not fails and bad == 0 and cov and bcov)
    return out
