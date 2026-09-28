"""Stage-15 data reduction from derived coefficients and frozen Stage-14 rows."""
import csv
import json
from pathlib import Path
from mpmath import mp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
mp.dps = 70
coeff = json.loads((ROOT / 'stage14_asymptotic_coefficients.json').read_text())
linf = mp.mpf(coeff['lambda_infinity'])
a = mp.mpf(coeff['a_1_over_N2'])
b = mp.mpf(coeff['b_1_over_N4'])
with (ROOT / 'lambdaN_asymptotics.csv').open(newline='', encoding='utf-8') as f:
    existing = list(csv.DictReader(f))

def fmt(x, n=35):
    return mp.nstr(x, n)

large, residual = [], []
for old in existing:
    N = int(old['N'])
    actual = mp.mpf(old['lambda_N'])
    p0 = linf
    p2 = linf + a / N**2
    p4 = p2 + b / N**4
    res2 = actual - p2
    res4 = actual - p4
    large.append({'N': N, 'lambda_N': fmt(actual, 42), 'origin': old['origin'],
                  'lambda_infinity': fmt(p0, 42), 'N_minus_2_prediction': fmt(p2, 42),
                  'N_minus_4_prediction': fmt(p4, 42),
                  'residual_after_N_minus_2': fmt(res2, 35),
                  'N4_times_residual': fmt(res2 * N**4, 35)})
    residual.append({'N': N, 'origin': old['origin'],
                     'N2_times_lambda_difference': fmt((actual-p0)*N**2, 35),
                     'N4_times_residual_after_N_minus_2': fmt(res2*N**4, 35),
                     'N6_times_residual_after_N_minus_4': fmt(res4*N**6, 35),
                     'absolute_residual_after_N_minus_2': fmt(abs(res2), 35),
                     'absolute_residual_after_N_minus_4': fmt(abs(res4), 35)})

for name, rows in [('lambdaN_largeN_data.csv', large),
                   ('lambdaN_residual_scaling.csv', residual)]:
    with (ROOT / name).open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)

xs = [1/int(r['N'])**2 for r in large]
ys = [float(r['lambda_N']) for r in large]
res = [float(r['N4_times_residual_after_N_minus_2']) for r in residual]
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.2, 3.4))
grid = [max(xs)*i/300 for i in range(301)]
ax1.plot(grid, [float(linf+a*x+b*x*x) for x in grid], color='black',
         label=r'$\lambda_\infty+a/N^2+b/N^4$')
ax1.scatter(xs[:5], ys[:5], label='certified finite $N$', color='#125e91', s=28)
ax1.scatter(xs[5:], ys[5:], label='numerical validation', color='#c45a20', s=28)
ax1.set(xlabel=r'$N^{-2}$', ylabel=r'$\lambda_N$')
ax1.legend(fontsize=7); ax1.grid(alpha=.2)
ax2.axhline(float(b), color='black', lw=1, label=r'$b=-0.28314327\ldots$')
ax2.scatter(xs[:5], res[:5], color='#125e91', s=28)
ax2.scatter(xs[5:], res[5:], color='#c45a20', s=28)
ax2.set(xlabel=r'$N^{-2}$', ylabel=r'$N^4(\lambda_N-\lambda_\infty-a/N^2)$')
ax2.legend(fontsize=7); ax2.grid(alpha=.2)
fig.tight_layout()
fig.savefig(ROOT / 'largeN_theory_vs_data.pdf')
plt.close(fig)
