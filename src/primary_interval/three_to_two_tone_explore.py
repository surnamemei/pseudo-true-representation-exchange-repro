import numpy as np
from scipy.optimize import minimize
from scipy.ndimage import minimum_filter

N = 21
t = np.arange(N) - (N - 1) / 2


def atom(w):
    return np.exp(1j * w * t)


def record(d, b, eps, phase):
    return atom(-d / N) + atom(d / N) + eps * np.exp(1j * phase) * atom(b / N)


def cost(ww, x):
    w = np.sort(np.asarray(ww))
    v = np.column_stack((atom(w[0]), atom(w[1])))
    if abs(w[1] - w[0]) < 1e-7:
        v = np.column_stack((atom(w[0]), 1j * t * atom(w[0])))
    c = np.linalg.lstsq(v, x, rcond=None)[0]
    e = x - v @ c
    return float(np.vdot(e, e).real)


def grid_cost(x, grid):
    v = np.exp(1j * np.outer(grid, t))
    h = v.conj() @ x
    diff = grid[:, None] - grid[None, :]
    den_d = np.sin(diff / 2)
    with np.errstate(divide='ignore', invalid='ignore'):
        g = np.sin(N * diff / 2) / den_d
    g[np.abs(den_d) < 1e-12] = N
    h1, h2 = h[:, None], h[None, :]
    den = N * N - g * g
    cap = (N * (abs(h1) ** 2 + abs(h2) ** 2) - 2 * g * (h1.conj() * h2).real) / np.maximum(den, 1e-12)
    out = np.vdot(x, x).real - cap
    # Row is the smaller frequency and column is the larger one.
    out[np.tril_indices(len(grid), 0)] = np.inf
    return out


def minima(d, b, eps, phase, m=256):
    x = record(d, b, eps, phase)
    grid = np.linspace(-np.pi, np.pi, m, endpoint=False)
    vals = grid_cost(x, grid)
    filt = minimum_filter(vals, size=3, mode='nearest')
    inds = np.argwhere((vals <= filt + 1e-10) & np.isfinite(vals))
    inds = sorted(inds, key=lambda ij: vals[tuple(ij)])[:30]
    found = []
    for ij in inds:
        init = grid[ij]
        res = minimize(cost, init, args=(x,), method='Nelder-Mead', options={'xatol': 1e-12, 'fatol': 1e-12, 'maxiter': 1000})
        w = np.sort((res.x + np.pi) % (2 * np.pi) - np.pi)
        c = cost(w, x)
        if not any(np.linalg.norm(w - a[1]) < 1e-3 for a in found):
            found.append((c, w))
    return sorted(found, key=lambda e: e[0])


if __name__ == '__main__':
    for d, b, phase in [(2, 10, 0), (2.5, 10, 0), (3, 10, 0), (2, 10, np.pi / 2), (2, 10, np.pi)]:
        for eps in [0.05, 0.15, 0.25, 0.4, 0.6, 0.9]:
            a = minima(d, b, eps, phase)
            print('d,b,ph,ep', d, b, round(phase, 3), eps, [(round(q, 5), np.round(w * N, 3)) for q, w in a[:3]], flush=True)
