"""Core numerics for the adversarial overnight validation (new code path).

Two-tone variable-projection objective in smooth (m, h) coordinates:
    u1 = m - h, u2 = m + h,
    span{e^{i u1 s}, e^{i u2 s}} = span{e^{ims} cos(hs), e^{ims} s sinc(hs)}.
For symmetric weights the two basis vectors are w-orthogonal, so
    J(m,h) = E - |C|^2/a - |S|^2/b,
which is smooth everywhere, even in h, 2*pi*N periodic in m, and equals the
coalescent closure at h = 0.  Everything here is numerical (float64); nothing
in this module is an interval certificate.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

TWO_PI = 2.0 * math.pi
SEED_ROOT = 20260928


# ----------------------------------------------------------------------------- signals
def times(n: int) -> np.ndarray:
    return (np.arange(n) - (n - 1) / 2) / n


def make_signal(n: int, z: float, b: float, eps: float, phi: float = math.pi, delta: float = 0.0,
                amp_minus: float = 1.0, amp_plus: float = 1.0) -> np.ndarray:
    """x_t = A- e^{-i(z s+delta)} + A+ e^{+i(z s+delta)} + eps e^{i phi} e^{i b s}."""
    s = times(n)
    return (amp_minus * np.exp(-1j * (z * s + delta)) + amp_plus * np.exp(1j * (z * s + delta))
            + eps * np.exp(1j * phi) * np.exp(1j * b * s))


def weights(name: str, n: int) -> np.ndarray:
    """Symmetric LS weights normalised to sum N (pre-registered definitions)."""
    name = name.lower()
    if name in ("rect", "rectangular"):
        w = np.ones(n)
    elif name == "hann":
        w = np.hanning(n)
    elif name == "hamming":
        w = np.hamming(n)
    elif name == "dpss":
        from scipy.signal.windows import dpss
        w = dpss(n, 2.5, sym=True)
        w = w * np.sign(w[n // 2])
    elif name == "hann2":
        w = np.hanning(n) ** 2
    elif name.startswith("homotopy:"):
        theta = float(name.split(":")[1])
        hw = np.hanning(n)
        hw = n * hw / hw.sum()
        w = (1 - theta) * np.ones(n) + theta * hw
    else:
        raise ValueError(name)
    w = np.asarray(w, float)
    return n * w / w.sum()


# ----------------------------------------------------------------------------- sinc helpers
def _sinc_derivs(x: np.ndarray):
    """sinc(x)=sin x/x and its first two derivatives, accurate near 0."""
    x = np.asarray(x, float)
    small = np.abs(x) < 0.2
    xs = np.where(small, 1.0, x)
    s0 = np.sin(xs) / xs
    s1 = (xs * np.cos(xs) - np.sin(xs)) / xs ** 2
    s2 = -s0 - 2.0 * s1 / xs
    x2 = x * x
    t0 = 1 - x2 / 6 + x2 ** 2 / 120 - x2 ** 3 / 5040 + x2 ** 4 / 362880 - x2 ** 5 / 39916800
    t1 = x * (-1 / 3 + x2 / 30 - x2 ** 2 / 840 + x2 ** 3 / 45360 - x2 ** 4 / 3991680)
    t2 = -1 / 3 + x2 / 10 - x2 ** 2 / 168 + x2 ** 3 / 6480 - x2 ** 4 / 443520 + x2 ** 5 * 132 / 6227020800
    return np.where(small, t0, s0), np.where(small, t1, s1), np.where(small, t2, s2)


# ----------------------------------------------------------------------------- objective
class Objective:
    """Weighted two-tone VarPro objective for one record y."""

    def __init__(self, y: np.ndarray, w: np.ndarray | None = None):
        self.y = np.asarray(y, complex)
        self.n = len(self.y)
        self.s = times(self.n)
        self.w = np.ones(self.n) if w is None else np.asarray(w, float)
        self.E = float(np.sum(self.w * np.abs(self.y) ** 2))
        self.wy = self.w * self.y
        self.period = TWO_PI * self.n

    # ---- point evaluation (vectorised over P points)
    def evaluate(self, m, h, order: int = 2):
        m = np.atleast_1d(np.asarray(m, float))
        h = np.atleast_1d(np.asarray(h, float))
        s = self.s[None, :]
        hs = h[:, None] * s
        ph = np.exp(-1j * m[:, None] * s)
        base = self.wy[None, :] * ph                      # w y e^{-ims}
        c = np.cos(hs)
        sn = np.sin(hs)
        s0, s1, s2 = _sinc_derivs(hs)
        sig = s * s0                                      # s sinc(hs) = sin(hs)/h
        wv = self.w[None, :]
        C = np.sum(base * c, axis=1)
        S = np.sum(base * sig, axis=1)
        a = np.sum(wv * c * c, axis=1)
        bb = np.sum(wv * sig * sig, axis=1)
        J = self.E - np.abs(C) ** 2 / a - np.abs(S) ** 2 / bb
        if order == 0:
            return J
        c1 = -s * sn
        c2 = -s * s * c
        sg1 = s * s * s1
        sg2 = s * s * s * s2
        mi = -1j * s
        Cm = np.sum(base * c * mi, axis=1)
        Sm = np.sum(base * sig * mi, axis=1)
        Ch = np.sum(base * c1, axis=1)
        Sh = np.sum(base * sg1, axis=1)
        a_h = np.sum(wv * 2 * c * c1, axis=1)
        b_h = np.sum(wv * 2 * sig * sg1, axis=1)

        def parts(X, Xm, Xh, nrm, nrm_h):
            F = np.abs(X) ** 2
            Fm = 2 * np.real(np.conj(X) * Xm)
            Fh = 2 * np.real(np.conj(X) * Xh)
            return F, Fm, Fh, Fm / nrm, Fh / nrm - F * nrm_h / nrm ** 2

        FC, FCm, FCh, gCm, gCh = parts(C, Cm, Ch, a, a_h)
        FS, FSm, FSh, gSm, gSh = parts(S, Sm, Sh, bb, b_h)
        grad = -np.stack([gCm + gSm, gCh + gSh], axis=1)
        if order == 1:
            return J, grad
        Cmm = np.sum(base * c * (-s * s), axis=1)
        Smm = np.sum(base * sig * (-s * s), axis=1)
        Chh = np.sum(base * c2, axis=1)
        Shh = np.sum(base * sg2, axis=1)
        Cmh = np.sum(base * c1 * mi, axis=1)
        Smh = np.sum(base * sg1 * mi, axis=1)
        a_hh = np.sum(wv * 2 * (c1 * c1 + c * c2), axis=1)
        b_hh = np.sum(wv * 2 * (sg1 * sg1 + sig * sg2), axis=1)

        def hess(X, Xm, Xh, Xmm, Xhh, Xmh, F, Fm, Fh, nrm, nrm_h, nrm_hh):
            Fmm = 2 * np.real(np.conj(X) * Xmm) + 2 * np.abs(Xm) ** 2
            Fhh = 2 * np.real(np.conj(X) * Xhh) + 2 * np.abs(Xh) ** 2
            Fmh = 2 * np.real(np.conj(X) * Xmh) + 2 * np.real(np.conj(Xh) * Xm)
            gmm = Fmm / nrm
            gmh = Fmh / nrm - Fm * nrm_h / nrm ** 2
            ghh = (Fhh / nrm - 2 * Fh * nrm_h / nrm ** 2 - F * nrm_hh / nrm ** 2
                   + 2 * F * nrm_h ** 2 / nrm ** 3)
            return gmm, gmh, ghh

        Cmm_, Cmh_, Chh_ = hess(C, Cm, Ch, Cmm, Chh, Cmh, FC, FCm, FCh, a, a_h, a_hh)
        Smm_, Smh_, Shh_ = hess(S, Sm, Sh, Smm, Shh, Smh, FS, FSm, FSh, bb, b_h, b_hh)
        H = -np.stack([np.stack([Cmm_ + Smm_, Cmh_ + Smh_], axis=1),
                       np.stack([Cmh_ + Smh_, Chh_ + Shh_], axis=1)], axis=1)
        return J, grad, H

    def J(self, m, h):
        return self.evaluate(m, h, 0)

    # ---- amplitudes / residual in the ordinary (u1,u2) chart
    def amplitudes(self, m: float, h: float):
        s = self.s
        c = np.exp(1j * m * s) * np.cos(h * s)
        sig = s * _sinc_derivs(h * s)[0]
        d = np.exp(1j * m * s) * sig
        gc = np.sum(self.w * np.conj(c) * self.y) / np.sum(self.w * np.abs(c) ** 2)
        gd = np.sum(self.w * np.conj(d) * self.y) / np.sum(self.w * np.abs(d) ** 2)
        resid = self.y - gc * c - gd * d
        if abs(h) > 1e-12:
            b2 = gc / 2 + gd / (2j * h)
            b1 = gc / 2 - gd / (2j * h)
        else:
            b1 = b2 = complex("inf")
        return np.array([b1, b2]), resid, (gc, gd)

    # ---- grid
    def grid(self, mm: int, mh: int):
        n = self.n
        mg = -math.pi * n + TWO_PI * n * np.arange(mm) / mm
        hg = np.linspace(0.0, math.pi * n / 2, mh)
        s = self.s
        hs = hg[:, None] * s[None, :]
        c = np.cos(hs)
        sig = s[None, :] * _sinc_derivs(hs)[0]
        a = (self.w[None, :] * c * c).sum(1)
        bb = (self.w[None, :] * sig * sig).sum(1)
        ph = np.exp(-1j * np.outer(s, mg))           # (n, mm)
        C = (c * self.wy[None, :]) @ ph               # (mh, mm)
        S = (sig * self.wy[None, :]) @ ph
        Jg = self.E - np.abs(C) ** 2 / a[:, None] - np.abs(S) ** 2 / bb[:, None]
        return mg, hg, Jg


# ----------------------------------------------------------------------------- geometry
def canon(m, h, n):
    """Canonical (m,h) with h in [0, pi N/2], m in [-pi N, pi N)."""
    period = TWO_PI * n
    m = np.asarray(m, float).copy()
    h = np.abs(np.asarray(h, float)).copy()
    h = np.mod(h, period)
    h = np.where(h > period / 2, period - h, h)
    flip = h > period / 4
    m = np.where(flip, m + period / 2, m)
    h = np.where(flip, period / 2 - h, h)
    m = np.mod(m + period / 2, period) - period / 2
    return m, h


def to_u(m, h, n):
    period = TWO_PI * n
    u1 = np.mod(np.asarray(m) - np.asarray(h) + period / 2, period) - period / 2
    u2 = np.mod(np.asarray(m) + np.asarray(h) + period / 2, period) - period / 2
    lo = np.minimum(u1, u2)
    hi = np.maximum(u1, u2)
    return lo, hi


def circ(a, b, n):
    period = TWO_PI * n
    d = np.mod(np.asarray(a) - np.asarray(b), period)
    return np.minimum(d, period - d)


def pair_distance(p, q, n):
    """Max-norm distance between unordered frequency pairs on the torus."""
    a1, a2 = p
    b1, b2 = q
    d1 = np.maximum(circ(a1, b1, n), circ(a2, b2, n))
    d2 = np.maximum(circ(a1, b2, n), circ(a2, b1, n))
    return np.minimum(d1, d2)


def separation(u1, u2, n):
    return circ(u1, u2, n)


# ----------------------------------------------------------------------------- refinement
def newton_refine(obj: Objective, m0, h0, tol_g: float = 1e-10, max_iter: int = 400,
                  max_step: float = 2.0):
    """Vectorised damped (absolute-eigenvalue) Newton minimisation in (m,h)."""
    x = np.stack([np.asarray(m0, float), np.asarray(h0, float)], axis=1).copy()
    P = len(x)
    active = np.ones(P, bool)
    iters = np.zeros(P, int)
    J, g, H = obj.evaluate(x[:, 0], x[:, 1])
    floor = 5e-14 * max(obj.E, 1e-300)
    last_step = np.full(P, np.inf)
    for it in range(max_iter):
        gn = np.linalg.norm(g, axis=1)
        done = (gn <= max(tol_g * max(obj.E, 1.0) * 1e-2, floor)) | (last_step <= 1e-11 * (1 + np.abs(x).max(1)))
        active &= ~done
        if not active.any():
            break
        idx = np.where(active)[0]
        Hs = 0.5 * (H[idx] + np.transpose(H[idx], (0, 2, 1)))
        ev, V = np.linalg.eigh(Hs)
        scale = np.maximum(np.abs(ev).max(1), 1e-300)
        ev_mod = np.maximum(np.abs(ev), 1e-10 * scale[:, None] + 1e-300)
        gl = np.einsum("pij,pi->pj", V, g[idx])
        step = -np.einsum("pij,pj->pi", V, gl / ev_mod)
        sn = np.linalg.norm(step, axis=1)
        fac = np.minimum(1.0, max_step / np.maximum(sn, 1e-300))
        step *= fac[:, None]
        alpha = np.ones(len(idx))
        accepted = np.zeros(len(idx), bool)
        slope = np.einsum("pi,pi->p", g[idx], step)
        newJ = J[idx].copy()
        newx = x[idx].copy()
        for _ in range(40):
            todo = ~accepted
            if not todo.any():
                break
            trial = x[idx[todo]] + alpha[todo, None] * step[todo]
            Jt = obj.evaluate(trial[:, 0], trial[:, 1], 0)
            ok = Jt <= J[idx[todo]] + 1e-4 * alpha[todo] * slope[todo] + 4e-16 * max(obj.E, 1.0)
            sub = np.where(todo)[0]
            acc = sub[ok]
            newJ[acc] = Jt[ok]
            newx[acc] = trial[ok]
            accepted[acc] = True
            alpha[sub[~ok]] *= 0.5
        moved = np.linalg.norm(newx - x[idx], axis=1)
        last_step[idx] = np.where(accepted, moved, 0.0)
        # points whose line search failed are at numerical stationarity
        x[idx] = newx
        iters[idx] += 1
        Jn, gnw, Hn = obj.evaluate(x[idx, 0], x[idx, 1])
        J[idx], g[idx], H[idx] = Jn, gnw, Hn
        stuck = ~accepted
        active[idx[stuck]] = False
    return x, J, g, H, iters


# ----------------------------------------------------------------------------- minima records
@dataclass
class Minimum:
    m: float
    h: float
    J: float
    u1: float
    u2: float
    sep: float
    grad: float
    eig_mh: tuple
    eig_u: tuple
    amp_max: float
    cls: str  # 'regular' or 'coalescent'

    def pair(self):
        return (self.u1, self.u2)

    def as_dict(self):
        return dict(m=self.m, h=self.h, J=self.J, u1=self.u1, u2=self.u2, sep=self.sep,
                    grad=self.grad, eig_mh_min=self.eig_mh[0], eig_mh_max=self.eig_mh[1],
                    eig_u_min=self.eig_u[0], eig_u_max=self.eig_u[1], amp_max=self.amp_max,
                    cls=self.cls)


T_U = np.array([[0.5, 0.5], [-0.5, 0.5]])  # rows (m,h), cols (u1,u2)


def make_minima(obj: Objective, x, J, g, H, dedup: float = 1e-3, sep_reg: float = 1.0):
    n = obj.n
    out: list[Minimum] = []
    order = np.argsort(J)
    for k in order:
        Hs = 0.5 * (H[k] + H[k].T)
        ev = np.linalg.eigvalsh(Hs)
        if ev[0] < -1e-9 * max(1.0, abs(J[k])):
            continue
        if not np.isfinite(J[k]) or np.linalg.norm(g[k]) > 1e-7 * max(obj.E, 1.0):
            continue
        m, h = canon(x[k, 0], x[k, 1], n)
        m, h = float(m), float(h)
        u1, u2 = to_u(m, h, n)
        u1, u2 = float(u1), float(u2)
        if any(pair_distance((u1, u2), q.pair(), n) < dedup for q in out):
            continue
        Hu = T_U.T @ Hs @ T_U
        evu = np.linalg.eigvalsh(Hu)
        sep = float(separation(u1, u2, n))
        if sep >= 1e-8:
            amps = obj.amplitudes(m, h)[0]
            amax = float(np.max(np.abs(amps)))
        else:
            amax = float("inf")
        out.append(Minimum(m, h, float(J[k]), u1, u2, sep, float(np.linalg.norm(g[k])),
                           (float(ev[0]), float(ev[1])), (float(evu[0]), float(evu[1])), amax,
                           "regular" if sep >= sep_reg else "coalescent"))
    return out


# ----------------------------------------------------------------------------- global search
LEVELS = {
    "P": dict(mm=1024, mh=257, k_grid=64, n_sat=48, n_rand=32, tol_g=1e-10),
    "R": dict(mm=2048, mh=513, k_grid=200, n_sat=96, n_rand=128, tol_g=1e-11),
    "P2": dict(mm=2048, mh=513, k_grid=128, n_sat=48, n_rand=64, tol_g=1e-10),
}


def level_params(level, n: int, **over):
    p = dict(LEVELS[level]) if isinstance(level, str) else dict(level)
    scale = n / 21.0
    p["mm"] = int(2 * math.ceil(p["mm"] * scale / 2))
    p["mh"] = int(math.ceil((p["mh"] - 1) * scale)) + 1
    p.update(over)
    return p


def grid_minima(obj: Objective, mm: int, mh: int, k: int):
    mg, hg, Jg = obj.grid(mm, mh)
    # pad: periodic in m; mirror at h=0; identification at h=pi N/2
    top = np.roll(Jg[-2:-1, :], mm // 2, axis=1)
    bot = Jg[1:2, :]
    P = np.vstack([bot, Jg, top])
    P = np.hstack([P[:, -1:], P, P[:, :1]])
    core = P[1:-1, 1:-1]
    ismin = np.ones_like(core, bool)
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            if di == 0 and dj == 0:
                continue
            ismin &= core <= P[1 + di:P.shape[0] - 1 + di, 1 + dj:P.shape[1] - 1 + dj]
    ii, jj = np.nonzero(ismin)
    vals = core[ii, jj]
    order = np.argsort(vals)[:k]
    return mg[jj[order]], hg[ii[order]], vals[order], (mg, hg, Jg)


def rng_for(*key):
    return np.random.default_rng(np.random.SeedSequence([SEED_ROOT] + [int(abs(k)) for k in key]))


def seeds_structured(n: int, z: float | None, b: float | None, n_sat: int, extra=()):
    ms, hs = [], []
    if z is not None:
        ms.append(0.0); hs.append(abs(z))
        for pair in ((-z, z), (-z, b), (z, b)) if b is not None else ((-z, z),):
            ms.append((pair[0] + pair[1]) / 2); hs.append(abs(pair[1] - pair[0]) / 2)
    period = TWO_PI * n
    for j in range(n_sat):
        v = -period / 2 + period * (j + 0.5) / n_sat
        ms.append(v / 2); hs.append(abs(v) / 2)
    for (u1, u2) in extra:
        ms.append((u1 + u2) / 2); hs.append(abs(u2 - u1) / 2)
    return np.array(ms, float), np.array(hs, float)


def global_search(obj: Objective, level="P", z=None, b=None, extra=(), rng_key=(0,), blind=False,
                  return_grid=False, **over):
    """Pre-registered global search protocol GS(level). Returns sorted list of Minimum."""
    p = level_params(level, obj.n, **over)
    gm, gh, gv, grid = grid_minima(obj, p["mm"], p["mh"], p["k_grid"])
    sm, sh = seeds_structured(obj.n, None if blind else z, None if blind else b, p["n_sat"], extra)
    rng = rng_for(*rng_key)
    period = TWO_PI * obj.n
    rm = rng.uniform(-period / 2, period / 2, p["n_rand"])
    rh = rng.uniform(0, period / 4, p["n_rand"])
    m0 = np.concatenate([gm, sm, rm])
    h0 = np.concatenate([gh, sh, rh])
    x, J, g, H, it = newton_refine(obj, m0, h0, tol_g=p["tol_g"])
    mins = make_minima(obj, x, J, g, H)
    if return_grid:
        return mins, grid
    return mins


def refine_from(obj: Objective, pairs, tol_g=1e-11):
    """Local refinement from given (u1,u2) seeds; returns list aligned with seeds."""
    pairs = np.asarray(pairs, float).reshape(-1, 2)
    m0 = (pairs[:, 0] + pairs[:, 1]) / 2
    h0 = (pairs[:, 1] - pairs[:, 0]) / 2
    x, J, g, H, it = newton_refine(obj, m0, h0, tol_g=tol_g, max_step=0.5)
    res = []
    for k in range(len(pairs)):
        m, h = canon(x[k, 0], x[k, 1], obj.n)
        u1, u2 = to_u(m, h, obj.n)
        Hs = 0.5 * (H[k] + H[k].T)
        ev = np.linalg.eigvalsh(Hs)
        evu = np.linalg.eigvalsh(T_U.T @ Hs @ T_U)
        sep = float(separation(u1, u2, obj.n))
        amax = float(np.max(np.abs(obj.amplitudes(float(m), float(h))[0]))) if sep > 1e-8 else float("inf")
        res.append(Minimum(float(m), float(h), float(J[k]), float(u1), float(u2), sep,
                           float(np.linalg.norm(g[k])), (float(ev[0]), float(ev[1])),
                           (float(evu[0]), float(evu[1])), amax,
                           "regular" if sep >= 1.0 else "coalescent"))
    return res


def merge_minima(n: int, lists, dedup=1e-3):
    allm = sorted([q for L in lists for q in L], key=lambda q: q.J)
    out = []
    for q in allm:
        if not any(pair_distance(q.pair(), r.pair(), n) < dedup for r in out):
            out.append(q)
    return out
