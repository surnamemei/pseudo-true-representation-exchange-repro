"""Stage 21 plot from archived Stage-18 continuation; no new continuation solve."""
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
LAM = 0.06598041112699136
C4 = -0.00101774929746980
CERT = 0.2481906301722774

with (ROOT / "p_bridge_numerical.csv").open(newline="", encoding="utf-8") as f:
    old = list(csv.DictReader(f))

rows = []
for rec in old:
    z = float(rec["z"])
    observed = float(rec["epsilon_cross"])
    quadratic = LAM*z*z
    quartic = quadratic+C4*z**4
    rows.append(dict(z=z, numerical_epsilon=observed,
                     quadratic_epsilon=quadratic, quartic_epsilon=quartic,
                     quadratic_relative_error=abs(quadratic-observed)/observed,
                     quartic_relative_error=abs(quartic-observed)/observed,
                     continuation_status="numerical_stationary_equal_cost",
                     globality_status="independently_certified_only_at_z_2" if z == 2 else "not_certified"))

with (ROOT / "z_continuation_vs_asymptotics.csv").open("w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)

fig, ax = plt.subplots(figsize=(6.6, 3.8))
z = np.array([r["z"] for r in rows])
ax.plot(z, [r["numerical_epsilon"] for r in rows], color="#116ba2", lw=2.2,
        label="Numerical A/B equal-cost continuation")
ax.plot(z, [r["quadratic_epsilon"] for r in rows], color="#b35b12", ls="--", lw=1.7,
        label=r"Quadratic $\lambda_{21}z^2$")
ax.plot(z, [r["quartic_epsilon"] for r in rows], color="#267a4c", ls=":", lw=2.1,
        label=r"Quartic $\lambda_{21}z^2+c_{21}z^4$")
ax.plot([2], [CERT], marker="s", linestyle="none", color="#5d267d", markersize=6,
        label=r"Independently certified $z=2$ crossing")
ax.set(xlabel=r"Normalized spacing $z=N\Delta$", ylabel=r"Crossing amplitude $\epsilon_c$")
ax.set_xlim(0, 2.10)
ax.set_ylim(0, .30)
ax.grid(alpha=.2)
ax.legend(loc="upper left", fontsize=7.5, frameon=False)
fig.tight_layout()
fig.savefig(ROOT / "z_continuation_vs_asymptotics.pdf")
plt.close(fig)

for key in ("quadratic_relative_error", "quartic_relative_error"):
    maximum = max(rows, key=lambda r: r[key])
    print(key, "max", maximum[key], "at z", maximum["z"], "z2", rows[-1][key])
print("rows", len(rows))
