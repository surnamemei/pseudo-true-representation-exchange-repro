"""Stage 8 (OPTIONAL_EXPLORATORY): real-sinusoid analogue of the finite-spacing exchange.

x_t = cos((w0-D)t) + cos((w0+D)t) - eps cos((w0+w3)t), t=-10..10, D=2/21, w3=10/21,
fitted by two real sinusoids (4 real amplitudes, 2 frequencies in (0, pi)).
"""
from __future__ import annotations

import csv
import math

import numpy as np
from scipy.optimize import minimize, brentq

from advrun import HERE, log, jdump, jload, status_update, pmap

SD = HERE / "stage8_real"
N = 21
T = np.arange(N) - (N - 1) / 2
DEL, W3 = 2.0 / 21, 10.0 / 21
W0S = {"pi_over_2": math.pi / 2, "6pi_over_21": 6 * math.pi / 21}


def record(w0, eps):
    return np.cos((w0 - DEL) * T) + np.cos((w0 + DEL) * T) - eps * np.cos((w0 + W3) * T)


def J(nu, x):
    cols = np.column_stack([np.cos(nu[0] * T), np.sin(nu[0] * T), np.cos(nu[1] * T), np.sin(nu[1] * T)])
    c, *_ = np.linalg.lstsq(cols, x, rcond=None)
    r = x - cols @ c
    return float(r @ r), c


def grid_minima(x, m=420, top=60):
    g = np.linspace(1e-3, math.pi - 1e-3, m)
    C = np.cos(np.outer(g, T))
    S = np.sin(np.outer(g, T))
    B = np.stack([C, S], axis=2)                      # (m, N, 2)
    i, j = np.triu_indices(m, 1)
    A = np.concatenate([B[i], B[j]], axis=2)          # (P, N, 4)
    G = np.einsum("pki,pkj->pij", A, A) + 1e-12 * np.eye(4)
    h = np.einsum("pki,k->pi", A, x)
    sol = np.linalg.solve(G, h[..., None])[..., 0]
    vals = x @ x - np.einsum("pi,pi->p", h, sol)
    Mv = np.full((m, m), np.inf)
    Mv[i, j] = vals
    loc = []
    for a, b in zip(i, j):
        v = Mv[a, b]
        nb = Mv[max(a - 1, 0):a + 2, max(b - 1, 0):b + 2]
        if v <= nb.min() + 1e-15:
            loc.append((v, g[a], g[b]))
    loc.sort()
    return loc[:top]


def refine(seed, x):
    f = lambda p: J(np.sort(np.clip(p, 1e-6, math.pi - 1e-6)), x)[0]
    r = minimize(f, seed, method="Nelder-Mead", options=dict(xatol=1e-11, fatol=1e-15, maxiter=4000))
    nu = np.sort(np.clip(r.x, 1e-6, math.pi - 1e-6))
    val = J(nu, x)[0]
    h = 1e-5
    H = np.zeros((2, 2))
    for a in range(2):
        for b in range(2):
            ea, eb = np.eye(2)[a] * h, np.eye(2)[b] * h
            H[a, b] = (J(nu + ea + eb, x)[0] - J(nu + ea - eb, x)[0] - J(nu - ea + eb, x)[0] + J(nu - ea - eb, x)[0]) / (4 * h * h)
    return nu, val, float(np.linalg.eigvalsh(H)[0])


def search(x, extra=()):
    seeds = [(a, b) for _, a, b in grid_minima(x)] + list(extra)
    out = []
    for s in seeds:
        nu, val, eig = refine(np.array(s), x)
        if eig > -1e-9 and not any(np.max(np.abs(nu - q[0])) < 1e-4 for q in out):
            out.append((nu, val, eig))
    out.sort(key=lambda q: q[1])
    return out


def point_task(args):
    w0, eps = args
    x = record(w0, eps)
    mins = search(x, extra=[(w0 - DEL, w0 + DEL), (w0 - DEL, w0 + W3), (w0 + DEL, w0 + W3)])
    return dict(eps=eps, minima=[dict(nu=q[0].tolist(), J=q[1], eig=q[2]) for q in mins[:5]])


def run(ctx=None):
    SD.mkdir(parents=True, exist_ok=True)
    done = SD / "DONE.json"
    if done.exists():
        return jload(done)
    status_update(8, status="running")
    eps_grid = np.linspace(0, 1.0, 101)
    out = {}
    for name, w0 in W0S.items():
        res = pmap(point_task, [(w0, float(e)) for e in eps_grid], desc=f"S8 {name}", stage=8)
        res.sort(key=lambda r: r["eps"])
        switches = []
        for a, b in zip(res, res[1:]):
            pa, pb = np.array(a["minima"][0]["nu"]), np.array(b["minima"][0]["nu"])
            if np.max(np.abs(pa - pb)) * N > 1.0:
                # continue both branches and solve the equal-cost point
                def gap(e):
                    x = record(w0, e)
                    return refine(pa, x)[1] - refine(pb, x)[1]
                try:
                    ec = brentq(gap, a["eps"], b["eps"], xtol=1e-12)
                    x = record(w0, ec)
                    P, Q = refine(pa, x), refine(pb, x)
                    mins = search(x, extra=[tuple(P[0]), tuple(Q[0])])
                    Js = 0.5 * (P[1] + Q[1])
                    others = [m for m in mins if np.max(np.abs(m[0] - P[0])) * N > 0.05 and np.max(np.abs(m[0] - Q[0])) * N > 0.05]
                    third = others[0][1] - Js if others else float("inf")
                    lower = [m for m in others if m[1] < Js - 1e-9]
                    sepP, sepQ = abs(P[0][1] - P[0][0]) * N, abs(Q[0][1] - Q[0][0]) * N
                    switches.append(dict(eps_c=ec, J_star=Js, P_nu=P[0].tolist(), Q_nu=Q[0].tolist(), P_eig=P[2], Q_eig=Q[2],
                                         P_sep_u=sepP, Q_sep_u=sepQ, third_margin=third,
                                         third_margin_rel=third / Js if Js > 0 else None, global_ok=not lower,
                                         qualifying=bool(not lower and min(sepP, sepQ) >= 1.0 and min(P[2], Q[2]) > 0
                                                         and third >= 0.02 * Js)))
                except Exception as exc:
                    switches.append(dict(error=repr(exc), eps_lo=a["eps"], eps_hi=b["eps"]))
        out[name] = dict(w0=w0, switches=switches, any_qualifying=any(s.get("qualifying") for s in switches))
    out["evidence_level"] = "OPTIONAL_EXPLORATORY (numerical grid + Nelder-Mead; not certified)"
    jdump(done, out)
    status_update(8, status="done")
    log(f"Stage 8: {[(k, v['any_qualifying']) for k, v in out.items() if isinstance(v, dict)]}", 8)
    return out
