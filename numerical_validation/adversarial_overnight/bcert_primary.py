"""Stage 3B primary: b-parameterised directed-interval certificate of the continuum A/B crossing.

Architecture follows the frozen Stage-14 continuum certificate (mpmath.iv, libmpi directed
rounding): joint Krawczyk isolation of (v_A, v_B, lambda) plus full-line exclusion on
[-100,100] by cost / gradient / concavity / monotone-derivative predicates, a coalescent
bound and an analytic |v|>=100 tail bound -- here uniformly for b in a tile T = [b0, b1].
The only structural additions are (i) b enters every enclosure as an interval, and
(ii) the Krawczyk residual uses a mean-value form in b, F(x_c,T) in F(x_c,b_c)+dF/db(x_c,T)(T-b_c).
Numerical candidates (float64) never decide acceptance.
"""
from __future__ import annotations

import json
import math
import time
from fractions import Fraction

import numpy as np
from mpmath import iv, mp

iv.dps = 50
mp.dps = 60
M2 = iv.mpf(1) / 12
M4 = iv.mpf(1) / 80
V2 = M4 - M2 * M2
TAIL_V = 100


def lo(x):
    return mp.mpf(x.a)


def hi(x):
    return mp.mpf(x.b)


def excl0(x):
    return lo(x) > 0 or hi(x) < 0


def s(x):
    return [mp.nstr(lo(x), 25), mp.nstr(hi(x), 25)]


# ----------------------------------------------------------------------------- kernel d(x)=2 sin(x/2)/x
SER = 40


def D(x, k=0):
    a, b = abs(lo(x)), abs(hi(x))
    if max(a, b) <= 1:
        out = iv.mpf(0)
        for m in range((k + 1) // 2, SER):
            out += (-1) ** m * x ** (2 * m - k) / (iv.mpf(2 * m + 1) * math.factorial(2 * m - k) * iv.mpf(2) ** (2 * m))
        return out + iv.mpf(["-1e-45", "1e-45"])
    out = 2 * iv.sin(x / 2) / x
    for j in range(1, k + 1):
        out = (2 * iv.sin(x / 2 + j * iv.pi / 2) / iv.mpf(2) ** j - j * out) / x
    return out


# ----------------------------------------------------------------------------- 3-variable jets (v, lam, b)
class J3:
    __slots__ = ("v", "g", "H")

    def __init__(self, v, g, H):
        self.v, self.g, self.H = v, g, H

    @staticmethod
    def c(x):
        z = iv.mpf(0)
        return J3(x if isinstance(x, iv.mpf) else iv.mpf(x), [z, z, z], [[z, z, z], [z, z, z], [z, z, z]])

    @staticmethod
    def var(x, k):
        j = J3.c(x)
        j.g = [iv.mpf(int(i == k)) for i in range(3)]
        return j

    def __add__(a, b):
        b = b if isinstance(b, J3) else J3.c(b)
        return J3(a.v + b.v, [a.g[i] + b.g[i] for i in range(3)],
                  [[a.H[i][j] + b.H[i][j] for j in range(3)] for i in range(3)])
    __radd__ = __add__

    def __neg__(a):
        return J3(-a.v, [-x for x in a.g], [[-x for x in r] for r in a.H])

    def __sub__(a, b):
        return a + (-(b if isinstance(b, J3) else J3.c(b)))

    def __rsub__(a, b):
        return (-a) + b

    def __mul__(a, b):
        if not isinstance(b, J3):
            b = b if isinstance(b, iv.mpf) else iv.mpf(b)
            return J3(a.v * b, [x * b for x in a.g], [[x * b for x in r] for r in a.H])
        return J3(a.v * b.v, [a.g[i] * b.v + a.v * b.g[i] for i in range(3)],
                  [[a.H[i][j] * b.v + a.v * b.H[i][j] + a.g[i] * b.g[j] + b.g[i] * a.g[j] for j in range(3)]
                   for i in range(3)])
    __rmul__ = __mul__

    def recip(a):
        r = 1 / a.v
        return J3(r, [-x * r * r for x in a.g],
                  [[2 * a.g[i] * a.g[j] * r * r * r - a.H[i][j] * r * r for j in range(3)] for i in range(3)])

    def __truediv__(a, b):
        if not isinstance(b, J3):
            b = b if isinstance(b, iv.mpf) else iv.mpf(b)
            return J3(a.v / b, [x / b for x in a.g], [[x / b for x in r] for r in a.H])
        return a * b.recip()


def Dj(x: J3, k=0):
    d0, d1, d2 = D(x.v, k), D(x.v, k + 1), D(x.v, k + 2)
    return J3(d0, [d1 * g for g in x.g], [[d1 * x.H[i][j] + d2 * x.g[i] * x.g[j] for j in range(3)] for i in range(3)])


def tangent_jet(v, lam, b):
    v, lam, b = J3.var(v, 0), J3.var(lam, 1), J3.var(b, 2)
    d, dp, d2 = Dj(v), Dj(v, 1), Dj(v, 2)
    q0 = J3.c(-M2) - lam * Dj(b)
    q1 = 2 * lam * Dj(b, 1)
    A = 1 - d * d - dp * dp / M2
    Bn = d2 - lam * Dj(b - v) - d * q0 + dp * q1 / (2 * M2)
    C = J3.c(M4) - 2 * lam * Dj(b, 2) + lam * lam - q0 * q0 - q1 * q1 / (4 * M2)
    R = C - Bn * Bn / A
    return R, A, Bn


def FJ(xA, xB, lam, b):
    RA, AA, BA = tangent_jet(xA, lam, b)
    RB, AB, BB = tangent_jet(xB, lam, b)
    F = [RA.g[0], RB.g[0], RA.v - RB.v]
    Jx = [[RA.H[0][0], iv.mpf(0), RA.H[0][1]],
          [iv.mpf(0), RB.H[0][0], RB.H[0][1]],
          [RA.g[0], -RB.g[0], RA.g[1] - RB.g[1]]]
    Fb = [RA.H[0][2], RB.H[0][2], RA.g[2] - RB.g[2]]
    return F, Jx, Fb, (RA, AA, BA, RB, AB, BB)


# ----------------------------------------------------------------------------- numerical candidates
def candidate(bc):
    from advtangent import TangentProblem
    anchor = (-4.060030216180106892789, 12.42390445039620145, 0.0665069703923955354559)
    x = np.array(anchor)
    steps = max(1, int(math.ceil(abs(bc - 10.0) / 0.02)))
    for j in range(1, steps + 1):
        bb = 10.0 + (bc - 10.0) * j / steps
        x = TangentProblem(21, bb, math.pi, 0.0, continuum=True).crossing(*x)
    h = 1e-4
    xp = TangentProblem(21, bc + h, math.pi, 0.0, continuum=True).crossing(*x)
    xm = TangentProblem(21, bc - h, math.pi, 0.0, continuum=True).crossing(*x)
    return x, (xp - xm) / (2 * h)


# ----------------------------------------------------------------------------- local (Krawczyk) part
def krawczyk(T, xc, rad, bc_float):
    bc = iv.mpf(bc_float)
    x0 = [iv.mpf(float(v)) for v in xc]
    X = [x0[i] + iv.mpf([-rad[i], rad[i]]) for i in range(3)]
    F0, J0, _, _ = FJ(x0[0], x0[1], x0[2], bc)
    _, _, Fb, _ = FJ(x0[0], x0[1], x0[2], T)
    _, JX, _, parts = FJ(X[0], X[1], X[2], T)
    mid = np.array([[float(mp.mpf(J0[i][j].mid)) for j in range(3)] for i in range(3)])
    Ci = np.linalg.inv(mid)
    C = [[iv.mpf(float(Ci[i, j])) for j in range(3)] for i in range(3)]
    Fx0 = [F0[i] + Fb[i] * (T - bc) for i in range(3)]
    Tm = [[iv.mpf(int(i == j)) - sum(C[i][k] * JX[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
    Y = [X[i] - x0[i] for i in range(3)]
    K = [x0[i] - sum(C[i][k] * Fx0[k] for k in range(3)) + sum(Tm[i][j] * Y[j] for j in range(3)) for i in range(3)]
    inside = all(lo(K[i]) > lo(X[i]) and hi(K[i]) < hi(X[i]) for i in range(3))
    contraction = max(hi(sum(abs(Tm[i][j]) * iv.mpf(rad[j]) for j in range(3)) / iv.mpf(rad[i])) for i in range(3))
    RA, AA, BA, RB, AB, BB = parts
    betaA, betaB = BA.v / AA.v, BB.v / AB.v
    m = dict(Rvv_A=lo(JX[0][0]), Rvv_B=lo(JX[1][1]), slope=lo(JX[2][2]), gram_A=lo(AA.v), gram_B=lo(AB.v),
             beta_A=lo(betaA), beta_B=-hi(betaB))
    ok = bool(inside and contraction < 1 and all(v > 0 for v in m.values()))
    return dict(ok=ok, inside=bool(inside), contraction=float(contraction), X=[s(x) for x in X], K=[s(x) for x in K],
                margins={k: mp.nstr(v, 12) for k, v in m.items()}), X, K


# ----------------------------------------------------------------------------- global part
class TileEval:
    def __init__(self, T, L):
        self.T, self.L = T, L
        self.q0 = -M2 - L * D(T)
        self.q1 = 2 * L * D(T, 1)
        self.C = M4 - 2 * L * D(T, 2) + L * L - self.q0 * self.q0 - self.q1 * self.q1 / (4 * M2)
        K = 46
        dc = [iv.mpf(0)] * (K + 3)
        for mm in range((K + 3) // 2):
            k = 2 * mm
            if k < len(dc):
                dc[k] = iv.mpf((-1) ** mm) / (math.factorial(k + 1) * iv.mpf(2) ** k)
        dp = [(k + 1) * dc[k + 1] for k in range(K + 1)]
        db = [(-1) ** k * D(T, k) / math.factorial(k) for k in range(K + 1)]
        ac, bc = [], []
        for k in range(K + 1):
            ac.append(iv.mpf(int(k == 0)) - sum(dc[j] * dc[k - j] + dp[j] * dp[k - j] / M2 for j in range(k + 1)))
            bc.append((k + 2) * (k + 1) * dc[k + 2] - L * db[k] - dc[k] * self.q0 + dp[k] * self.q1 / (2 * M2))
        self.ac, self.bc = ac[4:45], bc[2:43]

    @staticmethod
    def poly(coeff, x, r):
        out = iv.mpf(0)
        for k in reversed(range(r, len(coeff))):
            out = out * x + coeff[k] * (math.factorial(k) // math.factorial(k - r))
        return out + iv.mpf(["-1e-30", "1e-30"])

    def calc(self, a, c):
        x = iv.mpf([a, c])
        L, T = self.L, self.T
        if max(abs(a), abs(c)) <= 1:
            g, g1, g2 = [self.poly(self.ac, x, j) for j in range(3)]
            h, h1, h2 = [self.poly(self.bc, x, j) for j in range(3)]
        else:
            if a <= 0 <= c:
                return None
            d, d1, d2, d3, d4 = [D(x, j) for j in range(5)]
            g = 1 - d * d - d1 * d1 / M2
            g1 = -2 * d * d1 - 2 * d1 * d2 / M2
            g2 = -2 * (d1 * d1 + d * d2) - 2 * (d2 * d2 + d1 * d3) / M2
            w = T - x
            h = d2 - L * D(w) - d * self.q0 + d1 * self.q1 / (2 * M2)
            h1 = d3 + L * D(w, 1) - d1 * self.q0 + d2 * self.q1 / (2 * M2)
            h2 = d4 - L * D(w, 2) - d2 * self.q0 + d3 * self.q1 / (2 * M2)
        if lo(g) <= 0:
            return None
        R = self.C - h * h / g
        Rv = -2 * h * h1 / g + h * h * g1 / (g * g)
        Rvv = -2 * (h1 * h1 + h * h2) / g + 4 * h * h1 * g1 / (g * g) + h * h * g2 / (g * g) - 2 * h * h * g1 * g1 / (g * g * g)
        return R, Rv, Rvv

    def coalescent(self):
        q2 = -M4 + self.L * D(self.T, 2) - M2 * self.q0
        return self.C - q2 * q2 / V2

    def tail(self):
        V = iv.mpf(TAIL_V)
        dd = 2 / V
        dpb = (1 + dd) / V
        d2b = 1 / (2 * V) + 2 * dpb / V
        dbw = 2 / (V - self.T.b)
        AA = 1 - dd * dd - dpb * dpb / M2
        BB = d2b + abs(self.L).b * dbw + dd * abs(self.q0).b + dpb * abs(self.q1).b / (2 * M2.a)
        return self.C - BB * BB / AA


def global_exclusion(T, L, boxA, boxB, U, min_width=1e-9, budget=400000):
    te = TileEval(T, L)
    fixed = sorted([boxA, boxB])
    cuts = sorted(set([-100.0, 100.0] + [k / 4 for k in range(-400, 401)] + [fixed[0][0], fixed[0][1], fixed[1][0], fixed[1][1]]))
    cells = [(a, c) for a, c in zip(cuts, cuts[1:]) if not any(a >= f[0] and c <= f[1] for f in fixed)]
    leaves, un = [], []
    stack = list(reversed(cells))
    visits = 0
    cost_margin = None
    while stack:
        a, c = stack.pop()
        visits += 1
        z = te.calc(a, c)
        kind, margin = None, None
        if z is not None:
            R, Rv, Rvv = z
            if lo(R) > hi(U):
                kind, margin = "cost", lo(R) - hi(U)
                cost_margin = margin if cost_margin is None else min(cost_margin, margin)
            elif excl0(Rv):
                kind, margin = "gradient", min(abs(lo(Rv)), abs(hi(Rv)))
            elif hi(Rvv) < 0:
                kind, margin = "concave", -hi(Rvv)
            elif lo(Rvv) > 0:
                za, zc = te.calc(a, a), te.calc(c, c)
                if za and zc and ((lo(za[1]) > 0 and lo(zc[1]) > 0) or (hi(za[1]) < 0 and hi(zc[1]) < 0)):
                    kind, margin = "monotone", min(abs(lo(za[1])), abs(lo(zc[1])))
        if kind:
            leaves.append((a, c, kind, mp.nstr(margin, 8)))
            continue
        if c - a < min_width or visits > budget:
            un.append((a, c))
            if visits > budget:
                un.extend(stack)
                break
            continue
        m = (a + c) / 2
        stack.extend([(m, c), (a, m)])
    coal = te.coalescent()
    tail = te.tail()
    # root boxes: strict convexity in v over the whole box (uniqueness of the stationary point)
    roots = []
    for name, (a, c) in (("A", boxA), ("B", boxB)):
        z = te.calc(a, c)
        roots.append(dict(name=name, lo=a, hi=c, Rvv_lo=mp.nstr(lo(z[2]), 10) if z else None,
                          ok=bool(z and lo(z[2]) > 0)))
        leaves.append((a, c, "root_" + name, mp.nstr(lo(z[2]), 8) if z else "nan"))
    leaves.sort(key=lambda r: r[0])
    # exact coverage (Fractions of the binary endpoints)
    cov = Fraction(leaves[0][0]) == Fraction(-100) and Fraction(leaves[-1][1]) == Fraction(100) and all(
        Fraction(p[1]) == Fraction(q[0]) for p, q in zip(leaves, leaves[1:]))
    return dict(visits=visits, leaves=len(leaves), unresolved=len(un), unresolved_cells=un[:20],
                coverage_exact=bool(cov), coalescent_margin=mp.nstr(lo(coal) - hi(U), 12),
                tail_margin=mp.nstr(lo(tail) - hi(U), 12), cost_cell_margin_min=mp.nstr(cost_margin, 8) if cost_margin else None,
                roots=roots, U=s(U),
                ok=bool(not un and cov and lo(coal) > hi(U) and lo(tail) > hi(U) and all(r["ok"] for r in roots))), leaves


def certify_tile(b0: str, b1: str, rad_factor=1.5, write_leaves=None):
    t0 = time.time()
    T = iv.mpf([b0, b1])
    bc = (float(b0) + float(b1)) / 2
    xc, dx = candidate(bc)
    half = (float(b1) - float(b0)) / 2
    rec = dict(b=[b0, b1], candidate=list(map(float, xc)), dx_db=list(map(float, dx)))
    for attempt, fac in enumerate((rad_factor, 3.0, 6.0)):
        rad = [abs(dx[i]) * half * fac + 1e-9 * (1 + abs(xc[i])) for i in range(3)]
        loc, X, K = krawczyk(T, xc, rad, bc)
        if loc["ok"]:
            break
    rec["local"] = loc
    rec["local_attempts"] = attempt + 1
    if not loc["ok"]:
        rec.update(ok=False, reason="krawczyk", seconds=time.time() - t0)
        return rec, None
    L = iv.mpf([max(lo(K[2]), lo(X[2])), min(hi(K[2]), hi(X[2]))])
    boxA = (float(mp.mpf(X[0].a)), float(mp.mpf(X[0].b)))
    boxB = (float(mp.mpf(X[1].a)), float(mp.mpf(X[1].b)))
    # outward float boxes (floats chosen to contain the interval boxes)
    boxA = (math.nextafter(boxA[0], -math.inf), math.nextafter(boxA[1], math.inf))
    boxB = (math.nextafter(boxB[0], -math.inf), math.nextafter(boxB[1], math.inf))
    te = TileEval(T, L)
    zA = te.calc(*boxA)
    U = zA[0] if zA else iv.mpf([-1, 1])
    glob, leaves = global_exclusion(T, L, boxA, boxB, U)
    rec["global"] = glob
    rec["lambda_enclosure"] = s(L)
    rec["ok"] = bool(loc["ok"] and glob["ok"])
    rec["seconds"] = time.time() - t0
    return rec, leaves


if __name__ == "__main__":
    import sys
    r, lv = certify_tile(sys.argv[1], sys.argv[2])
    print(json.dumps({k: v for k, v in r.items() if k != "global"}, indent=1, default=str))
    print(json.dumps({k: v for k, v in r.get("global", {}).items() if k != "unresolved_cells"}, indent=1, default=str))
