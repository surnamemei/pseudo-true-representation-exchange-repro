"""Independent tangent-level (small-z) and continuum losses by direct realified least squares.

Finite N:   R_N(v) = min_{c in R^5} || q - [1, i, i s, e^{ivs}, i e^{ivs}] c ||^2   over samples s_t = t/N.
Continuum:  same with Gauss-Legendre quadrature on [-1/2, 1/2] (weights included).
Coalescent: basis [1, i, s, i s, s^2] (complex 1, complex s, real s^2).
Target:     q(s) = -(s + kappa)^2 + lam * e^{i phi} e^{i b s}.
These are numerical evaluations (float64 scans, mpmath polishing), not interval proofs.
"""
from __future__ import annotations

import math

import numpy as np
from mpmath import mp
from scipy.optimize import minimize_scalar


class TangentProblem:
    def __init__(self, N=21, b=10.0, phi=math.pi, kappa=0.0, continuum=False, gl_nodes=256):
        self.N, self.b, self.phi, self.kappa, self.continuum = N, float(b), float(phi), float(kappa), continuum
        if continuum:
            x, w = np.polynomial.legendre.leggauss(gl_nodes)
            self.s = 0.5 * x
            self.w = 0.5 * w
        else:
            self.s = (np.arange(N) - (N - 1) / 2) / N
            self.w = np.ones(N)
        self.sw = np.sqrt(self.w)
        s = self.s
        self.q_strong = -(s + self.kappa) ** 2
        self.q_weak = np.exp(1j * self.phi) * np.exp(1j * self.b * s)
        self.base_cols = np.stack([np.ones_like(s), 1j * np.ones_like(s), 1j * s], axis=1)
        self.coal_cols = np.stack([np.ones_like(s) + 0j, 1j * np.ones_like(s), s + 0j, 1j * s, s * s + 0j], axis=1)
        self.period = 2 * math.pi * N if not continuum else np.inf

    def target(self, lam):
        return self.q_strong + lam * self.q_weak

    @staticmethod
    def _real(M):
        return np.concatenate([M.real, M.imag], axis=-2)

    def R(self, v, lam):
        """Vectorised regular-chart tangent loss over an array of v."""
        v = np.atleast_1d(np.asarray(v, float))
        s = self.s
        ev = np.exp(1j * np.outer(v, s))                        # (P, n)
        cols = np.concatenate([np.broadcast_to(self.base_cols, (len(v),) + self.base_cols.shape),
                               ev[:, :, None], (1j * ev)[:, :, None]], axis=2)   # (P, n, 5)
        cols = cols * self.sw[None, :, None]
        q = self.target(lam) * self.sw
        A = self._real(cols)                                    # (P, 2n, 5)
        qr = np.concatenate([q.real, q.imag])
        G = np.einsum("pki,pkj->pij", A, A)
        rhs = np.einsum("pki,k->pi", A, qr)
        try:
            c = np.linalg.solve(G, rhs[..., None])[..., 0]
        except np.linalg.LinAlgError:
            c = np.stack([np.linalg.lstsq(G[i], rhs[i], rcond=None)[0] for i in range(len(v))])
        res = qr[None, :] - np.einsum("pki,pi->pk", A, c)
        return np.sum(res * res, axis=1), c

    def coal(self, lam):
        cols = self.coal_cols * self.sw[:, None]
        q = self.target(lam) * self.sw
        A = np.concatenate([cols.real, cols.imag], axis=0)
        qr = np.concatenate([q.real, q.imag])
        c, *_ = np.linalg.lstsq(A, qr, rcond=None)
        r = qr - A @ c
        return float(r @ r)

    def beta(self, v, lam):
        _, c = self.R([v], lam)
        return complex(c[0, 3], c[0, 4])

    # ---- derivatives by central differences (float64)
    def Rv(self, v, lam, h=2e-5):
        return (self.R([v + h], lam)[0][0] - self.R([v - h], lam)[0][0]) / (2 * h)

    def Rvv(self, v, lam, h=2e-4):
        r = self.R([v - h, v, v + h], lam)[0]
        return (r[0] - 2 * r[1] + r[2]) / h ** 2

    def local_min(self, v0, lam, width=0.3):
        res = minimize_scalar(lambda x: float(self.R([x], lam)[0][0]), bounds=(v0 - width, v0 + width),
                              method="bounded", options={"xatol": 1e-13})
        return float(res.x), float(res.fun)

    # ---- crossing solve in float64 (then optional mpmath polish)
    def crossing(self, vA, vB, lam, iters=60):
        x = np.array([vA, vB, lam], float)
        for _ in range(iters):
            F = self._F(x)
            Jm = np.zeros((3, 3))
            for j, h in enumerate((1e-6, 1e-6, 1e-8)):
                e = np.zeros(3); e[j] = h
                Jm[:, j] = (self._F(x + e) - self._F(x - e)) / (2 * h)
            dx = np.linalg.solve(Jm, -F)
            x = x + dx
            if np.max(np.abs(dx)) < 1e-13:
                break
        return x

    def _F(self, x):
        vA, vB, lam = x
        r = self.R([vA, vB], lam)[0]
        return np.array([self.Rv(vA, lam), self.Rv(vB, lam), r[0] - r[1]])

    # ---- full-period / full-line scan of regular minima
    def scan_minima(self, lam, lo=None, hi=None, npts=8192, exclude=0.05):
        if lo is None:
            lo, hi = (-self.period / 2, self.period / 2) if not self.continuum else (-100.0, 100.0)
        grid = np.linspace(lo, hi, npts)
        vals = np.concatenate([self.R(grid[i:i + 2048], lam)[0] for i in range(0, len(grid), 2048)])
        vals = np.where(np.abs(grid) < exclude, np.inf, vals)
        idx = [i for i in range(1, len(grid) - 1) if vals[i] <= vals[i - 1] and vals[i] <= vals[i + 1]
               and np.isfinite(vals[i])]
        mins = []
        step = grid[1] - grid[0]
        for i in idx:
            x, f = self.local_min(grid[i], lam, width=1.5 * step)
            if abs(x) > exclude:
                mins.append((x, f))
        # explicit probes close to the central chart
        for a, b in ((exclude, 0.6), (-0.6, -exclude)):
            res = minimize_scalar(lambda x: float(self.R([x], lam)[0][0]), bounds=(a, b), method="bounded",
                                  options={"xatol": 1e-12})
            if abs(res.x - a) > 1e-6 and abs(res.x - b) > 1e-6:
                mins.append((float(res.x), float(res.fun)))
        mins.sort(key=lambda t: t[0])
        ded = []
        for x, f in mins:
            if not ded or abs(x - ded[-1][0]) > 1e-4:
                ded.append((x, f))
            elif f < ded[-1][1]:
                ded[-1] = (x, f)
        return sorted(ded, key=lambda t: t[1])


def analyse_crossing(tp: TangentProblem, vA, vB, lam, scan_pts=8192):
    """Crossing + curvatures + coefficients + slope + margins (numerical)."""
    vA, vB, lam = tp.crossing(vA, vB, lam)
    rA, rB = tp.R([vA, vB], lam)[0]
    Rstar = 0.5 * (rA + rB)
    cA, cB = tp.Rvv(vA, lam), tp.Rvv(vB, lam)
    bA, bB = tp.beta(vA, lam), tp.beta(vB, lam)
    d = 1e-6
    slope = ((tp.R([vA], lam + d)[0][0] - tp.R([vB], lam + d)[0][0])
             - (tp.R([vA], lam - d)[0][0] - tp.R([vB], lam - d)[0][0])) / (2 * d)
    # envelope: optimise v at lam +/- d as well (first-order identical); use fixed-v difference
    coal = tp.coal(lam)
    mins = tp.scan_minima(lam, npts=scan_pts)
    others = [(x, f) for x, f in mins if abs(x - vA) > 0.05 and abs(x - vB) > 0.05]
    third = others[0] if others else (float("nan"), float("inf"))
    best = mins[0] if mins else (float("nan"), float("nan"))
    return dict(lam=float(lam), vA=float(vA), vB=float(vB), R_star=float(Rstar), R_diff=float(rA - rB),
                Rvv_A=float(cA), Rvv_B=float(cB), beta_A_re=bA.real, beta_A_im=bA.imag,
                beta_B_re=bB.real, beta_B_im=bB.imag, beta_A_abs=abs(bA), beta_B_abs=abs(bB),
                slope=float(slope), coalescent_cost=float(coal), coalescent_margin=float(coal - Rstar),
                third_v=float(third[0]), third_margin=float(third[1] - Rstar), n_regular_minima=len(mins),
                best_scan_v=float(best[0]), best_scan_is_AB=bool(min(abs(best[0] - vA), abs(best[0] - vB)) < 0.05),
                separation=float(abs(vB - vA)))


# ----------------------------------------------------------------------------- continuum tail bound
def continuum_tail_margin(lam, b, Rstar, V=100.0):
    mu2 = 1.0 / 12
    mu4 = 1.0 / 80

    def d(x, k=0):
        x = float(x)
        if abs(x) < 1e-3:
            return [1 - x * x / 24, -x / 12, -1 / 12][k]
        s, c = math.sin(x / 2), math.cos(x / 2)
        v0 = 2 * s / x
        if k == 0:
            return v0
        v1 = c / x - 2 * s / x ** 2
        if k == 1:
            return v1
        return -s / (2 * x) - 2 * c / x ** 2 + 4 * s / x ** 3

    q0 = -mu2 - lam * d(b)
    q1 = 2 * lam * d(b, 1)
    base = mu4 - 2 * lam * d(b, 2) + lam * lam - q0 * q0 - q1 * q1 / (4 * mu2)
    dmax = 2 / V
    dpmax = 1 / V + 2 / V ** 2
    ddmax = 0.5 / V + 2 / V ** 2 + 4 / V ** 3
    bmax = ddmax + abs(lam) * 2 / (V - abs(b)) + dmax * abs(q0) + dpmax * abs(q1) / (2 * mu2)
    amin = 1 - dmax * dmax - dpmax * dpmax / mu2
    return base - bmax * bmax / amin - Rstar, base
