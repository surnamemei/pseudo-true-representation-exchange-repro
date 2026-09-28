"""Variable-projection two-tone LS for Stage 5, with centered time and u=Nw."""
import numpy as np
from scipy.optimize import minimize
from scipy.ndimage import minimum_filter


def times(n):
    return (np.arange(n) - (n - 1) / 2) / n


def signal(n, z, b, epsilon, phase=np.pi, ratio=1.0, strong_shift=0.0):
    s = times(n)
    return (np.exp(-1j*z*s) + ratio*np.exp(1j*(z+strong_shift)*s)
            + epsilon*np.exp(1j*phase)*np.exp(1j*b*s))


def objective(u, y, s, want_detail=False):
    u = np.asarray(u, dtype=float)
    v = np.exp(1j*np.outer(s, u))
    if abs(u[1]-u[0]) < 1e-6:
        # Stable confluent objective; the gradient is not used near coalescence.
        v[:, 1] = 1j*s*v[:, 0]
    c = np.linalg.lstsq(v, y, rcond=None)[0]
    residual = y-v@c
    val = float(np.vdot(residual, residual).real)
    if want_detail:
        return val, c, residual, v
    if abs(u[1]-u[0]) < 1e-6:
        return val, np.array([0.0, 0.0])
    grad = -2*np.real(np.sum(np.conj(residual[:, None])*(1j*s[:, None]*v)*c[None, :], axis=0))
    return val, grad


def refine(seed, y, s, maxiter=250):
    def fun(u):
        return objective(u, y, s)
    result = minimize(fun, np.asarray(seed, dtype=float), jac=True, method="BFGS",
                      options={"gtol": 2e-10, "maxiter": maxiter})
    u = np.sort(result.x)
    val, c, residual, v = objective(u, y, s, True)
    # Frequency torus: the canonical u range is [-pi*N, pi*N).
    h = np.empty((2, 2))
    step = 2e-4
    for k in range(2):
        du = np.eye(2)[k]*step
        h[:, k] = (objective(u+du, y, s)[1]-objective(u-du, y, s)[1])/(2*step)
    h = (h+h.T)/2
    return dict(u=u, J=val, coefficients=c, residual=residual,
                hessian_eigs=np.linalg.eigvalsh(h),
                gram_cond=np.linalg.cond(v.conj().T@v),
                converged=np.linalg.norm(objective(u, y, s)[1]) < 1e-6)


def grid_minima(y, s, grid_size=128, top=16):
    n = len(s)
    grid = np.linspace(-np.pi*n, np.pi*n, grid_size, endpoint=False)
    v = np.exp(1j*np.outer(grid, s))
    h = v.conj()@y
    g = v.conj()@v.T
    g = np.real(g)
    den = n*n-g*g
    h1, h2 = h[:, None], h[None, :]
    cap = (n*(np.abs(h1)**2+np.abs(h2)**2)-2*g*np.real(np.conj(h1)*h2))/np.maximum(den,1e-10)
    vals = np.vdot(y,y).real-cap
    vals[np.tril_indices(grid_size, 0)] = np.inf
    loc = minimum_filter(vals, size=3, mode="nearest")
    inds = np.argwhere(np.isfinite(vals)&(vals <= loc+1e-9))
    inds = sorted(inds, key=lambda ij: vals[tuple(ij)])[:top]
    roots = []
    for ij in inds:
        fit = refine(grid[ij], y, s)
        if not any(np.linalg.norm(fit["u"]-a["u"])<.02 for a in roots):
            roots.append(fit)
    return sorted(roots, key=lambda a:a["J"])
