"""Stage-20 tangent diagnostics for strong-pair relative phase.

Numerical scans are candidate searches, not interval certificates.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from mpmath import mp
from scipy.optimize import minimize_scalar

from stage14_continuum import D as Dc, MU2, R as Rc

ROOT = Path(__file__).resolve().parent
N = 21
B = 10
mp.dps = 50
S = [mp.mpf(t) / N for t in range(-10, 11)]
S2 = mp.fsum(s * s for s in S)
S4 = mp.fsum(s**4 for s in S)
V2 = S4 - S2**2 / N


def D(x, deriv=0):
    if deriv == 0:
        return mp.fsum(mp.cos(x * s) for s in S)
    if deriv == 1:
        return -mp.fsum(s * mp.sin(x * s) for s in S)
    if deriv == 2:
        return -mp.fsum(s * s * mp.cos(x * s) for s in S)
    raise ValueError(deriv)


def P(v):
    d = D(v)
    g = N - d * d / N
    return S2 - D(v, 1) ** 2 / g


def R0(v, lam):
    d, dp = D(v), D(v, 1)
    q0 = -S2 - lam * D(B)
    q1 = 2 * lam * D(B, 1)
    a = N - d * d / N - dp * dp / S2
    bb = D(v, 2) - lam * D(B - v) - d * q0 / N + dp * q1 / (2 * S2)
    c = S4 - 2 * lam * D(B, 2) + N * lam * lam - q0 * q0 / N - q1 * q1 / (4 * S2)
    return c - bb * bb / a


def H(v):
    d, dp = D(v), D(v, 1)
    q0 = -D(B)
    q1 = 2 * D(B, 1)
    a = N - d * d / N - dp * dp / S2
    bb = -D(B - v) - d * q0 / N + dp * q1 / (2 * S2)
    c = N - q0 * q0 / N - q1 * q1 / (4 * S2)
    return c - bb * bb / a


def K_coal():
    return N - D(B) ** 2 / N - D(B, 1) ** 2 / S2 - (D(B, 2) + S2 * D(B) / N) ** 2 / V2


def Rk(v, lam, kap):
    return R0(v, lam) + 4 * kap * kap * P(v)


def solve_kappa(kap, continuum=False, seed=None):
    kap = mp.mpf(kap)
    if continuum:
        def PP(v):
            return MU2 - Dc(v, 1) ** 2 / (1 - Dc(v) ** 2)
        def RR(v, lam):
            return Rc(v, lam) + 4 * kap * kap * PP(v)
        if seed is None:
            seed = (-4.06003, 12.42390, 0.06650697)
    else:
        PP, RR = P, lambda v, lam: Rk(v, lam, kap)
        if seed is None:
            seed = (-4.15344, 12.46809, 0.06598041)
    va, vb, lam = mp.findroot(
        (lambda a, b, l: mp.diff(lambda v: RR(v, l), a),
         lambda a, b, l: mp.diff(lambda v: RR(v, l), b),
         lambda a, b, l: RR(a, l) - RR(b, l)),
        seed, tol=mp.mpf("1e-32"), maxsteps=40,
    )
    curvature_a = mp.diff(lambda v: RR(v, lam), va, 2)
    curvature_b = mp.diff(lambda v: RR(v, lam), vb, 2)
    slope = mp.diff(lambda l: RR(va, l) - RR(vb, l), lam)
    rstar = RR(va, lam)
    if continuum:
        co = Rc(mp.mpf("2"), lam)  # overwritten below by exact coalescent formula
        from stage14_continuum import coal
        co = coal(lam)
    else:
        # The shifted linear term belongs to the confluent {1,s} class.
        q0 = -S2 - lam * D(B)
        q2 = -S4 + lam * D(B, 2) - S2 * q0 / N
        c = S4 - 2 * lam * D(B, 2) + N * lam * lam - q0 * q0 / N - (2 * lam * D(B, 1)) ** 2 / (4 * S2)
        co = c - q2 * q2 / V2
    return dict(kappa=float(kap),lambda_c=float(lam),v_A=float(va),v_B=float(vb),
                R_star=float(rstar),Rvv_A=float(curvature_a),Rvv_B=float(curvature_b),
                slope=float(slope),separation=float(vb-va),
                coalescent_margin=float(co-rstar),
                P_A=float(PP(va)),P_B=float(PP(vb)))


def finite_scan(lam, kap, grid_size=4097):
    # Full normalized frequency period, with near-zero and endpoints excluded.
    grid = np.linspace(-np.pi * N, np.pi * N, grid_size)
    vals = np.full(len(grid), np.inf)
    for i, v in enumerate(grid):
        if abs(v) > 0.015:
            vals[i] = float(Rk(mp.mpf(float(v)), mp.mpf(lam), mp.mpf(kap)))
    ids = [i for i in range(1, len(grid)-1) if vals[i] < vals[i-1] and vals[i] < vals[i+1]]
    out=[]
    for i in ids:
        a,b = grid[i-1],grid[i+1]
        res = minimize_scalar(lambda x:float(Rk(mp.mpf(x),mp.mpf(lam),mp.mpf(kap))),bounds=(a,b),method="bounded",options={"xatol":1e-8})
        if abs(res.x) > 0.12:  # reject artificial minima at the excised singular grid cell
            out.append((res.x,res.fun))
    return sorted(out,key=lambda x:x[1])


def fixed_scan():
    # P(v)+mu H(v) is independent of delta after mu=(lambda/(2 sin delta))^2.
    grid=np.linspace(-np.pi*N,np.pi*N,6001)
    p=np.array([float(P(mp.mpf(float(v)))) if abs(v)>.02 else np.nan for v in grid])
    h=np.array([float(H(mp.mpf(float(v)))) if abs(v)>.02 else np.nan for v in grid])
    results=[]
    for mu in np.geomspace(1e-5,1e5,81):
        vals=p+mu*h
        ids=[i for i in range(1,len(grid)-1) if np.isfinite(vals[i]) and vals[i]<vals[i-1] and vals[i]<vals[i+1]]
        roots=[]
        for i in ids:
            res=minimize_scalar(lambda x:float(P(mp.mpf(x))+mp.mpf(mu)*H(mp.mpf(x))),bounds=(grid[i-1],grid[i+1]),method="bounded")
            if abs(res.x)>0.08:
                roots.append((res.x,res.fun))
        # The global branch can lie much closer to zero than one grid step.
        for low,high in ((1e-9,3),(-3,-1e-9)):
            res=minimize_scalar(lambda x:float(P(mp.mpf(x))+mp.mpf(mu)*H(mp.mpf(x))),bounds=(low,high),method="bounded",options={"xatol":1e-10})
            if abs(res.x)<2.95 and not any(abs(res.x-v)<.01 for v,_ in roots):
                roots.append((res.x,res.fun))
        roots.sort(key=lambda q:q[1])
        co=float(mu*K_coal())
        results.append(dict(mu=float(mu),lambda_over_2sin_delta=float(np.sqrt(mu)),best_v=roots[0][0] if roots else None,
                            best_cost=roots[0][1] if roots else None,second_v=roots[1][0] if len(roots)>1 else None,
                            second_gap=roots[1][1]-roots[0][1] if len(roots)>1 else None,
                            coalescent_gap=co-roots[0][1] if roots else None,local_minima=len(roots)))
    return results


def main():
    kappas=[0,-.05,.05,-.10,.10,-.20,.20,-.30,.30,-.50,.50]
    rows=[solve_kappa(k) for k in kappas]
    # The fixed-kappa root is numerically the same under sign reversal.
    for row in rows:
        best=finite_scan(row['lambda_c'],row['kappa'],2049)
        row['sampled_local_minima']=len(best)
        row['sampled_third_margin']=min((cost-row['R_star'] for v,cost in best if abs(v-row['v_A'])>.1 and abs(v-row['v_B'])>.1),default=float('nan'))
        row['sampled_best_v']=best[0][0]
        row['sampled_global_AB']=bool(abs(best[0][0]-row['v_A'])<.1 or abs(best[0][0]-row['v_B'])<.1)
    with (ROOT/'kappa_crossing_map.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    cont=[solve_kappa(k,True) for k in [0,.05,.1,.2,.3,.5]]
    fixed=fixed_scan()
    (ROOT/'stage20_phase_diagnostics.json').write_text(json.dumps(dict(kappa_finite=rows,kappa_continuum=cont,fixed_scan=fixed,
                        K_coalescent=float(K_coal()),S2=float(S2)),indent=2),encoding='utf-8')
    print(json.dumps(dict(kappa_finite=rows,kappa_continuum=cont,
                          fixed_transitions=[(r['mu'],r['best_v'],r['coalescent_gap']) for r in fixed[::8]],
                          K_coalescent=float(K_coal())),indent=2))


if __name__=='__main__': main()
