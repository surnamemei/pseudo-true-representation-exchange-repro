"""Consistency check (not a new experiment): evaluate the manuscript's closed-form tangent loss
R_{N,chi}(v,lam) = R_{N,0}(v,lam) + 4 chi^2 P_N(v)   (eq. phase-chi / S23)
at N=21, b=10, phi=pi, and solve the A/B equal-cost crossing lambda_21(chi) by local Newton/root finding
from the certified chi=0 roots. Compare with (i) the envelope coefficient k_21=1.9283886651 and
(ii) the Stage-2 numerical tangent values."""
import mpmath as mp

mp.mp.dps = 40
N = 21
b = mp.mpf(10)
ts = [mp.mpf(t) for t in range(-(N - 1) // 2, (N - 1) // 2 + 1)]
s = [t / N for t in ts]
S2 = mp.fsum(x ** 2 for x in s)
S4 = mp.fsum(x ** 4 for x in s)


def D(q):
    return mp.fsum(mp.cos(q * x) for x in s)


def D1(q):
    return -mp.fsum(x * mp.sin(q * x) for x in s)


def D2(q):
    return -mp.fsum(x ** 2 * mp.cos(q * x) for x in s)


Db, D1b, D2b = D(b), D1(b), D2(b)


def R0(v, lam):
    q0 = -S2 - lam * Db
    q1 = 2 * lam * D1b
    q3 = D2(v) - lam * D(b - v)
    G = N - D(v) ** 2 / N - D1(v) ** 2 / S2
    H = q3 - D(v) * q0 / N + D1(v) * q1 / (2 * S2)
    return S4 - 2 * lam * D2b + N * lam ** 2 - q0 ** 2 / N - q1 ** 2 / (4 * S2) - H ** 2 / G


def P(v):
    return S2 - D1(v) ** 2 / (N - D(v) ** 2 / N)


def R(v, lam, chi):
    return R0(v, lam) + 4 * chi ** 2 * P(v)


def crossing(chi, vA, vB, lam):
    f = lambda a, bb, l: [mp.diff(lambda x: R(x, l, chi), a), mp.diff(lambda x: R(x, l, chi), bb), R(a, l, chi) - R(bb, l, chi)]
    return mp.findroot(f, (vA, vB, lam))


vA, vB, lam = mp.mpf("-4.15344048664877"), mp.mpf("12.4680948443009"), mp.mpf("0.0659804111269914")
sol = crossing(mp.mpf(0), vA, vB, lam)
lam0 = sol[2]
print("chi=0:", [mp.nstr(x, 16) for x in sol])
# envelope coefficient
Gam = mp.diff(lambda l: R(sol[0], l, 0) - R(sol[1], l, 0), lam0)
k = 4 * (P(sol[1]) - P(sol[0])) / Gam
print("Gamma =", mp.nstr(Gam, 12), " P_A =", mp.nstr(P(sol[0]), 12), " P_B =", mp.nstr(P(sol[1]), 12), " k =", mp.nstr(k, 12))
stage2 = {0.05: 0.0709172933824949, 0.1: 0.08551102766366643, 0.2: 0.1322256614331775}
prev = sol
for chi in (0.01, 0.02, 0.05, 0.1, 0.2):
    c = mp.mpf(chi)
    prev = crossing(c, prev[0], prev[1], prev[2])
    shift = prev[2] - lam0
    line = f"chi={chi}: vA={mp.nstr(prev[0], 10)} vB={mp.nstr(prev[1], 10)} lam={mp.nstr(prev[2], 14)} shift/chi^2={mp.nstr(shift / c ** 2, 10)}"
    if chi in stage2:
        line += f"  stage2 lam={stage2[chi]:.14f}  diff={float(prev[2]) - stage2[chi]:.3e}"
    print(line)
