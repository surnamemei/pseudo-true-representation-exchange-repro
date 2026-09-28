"""POST_HOC_SENSITIVITY figure: fixed-delta global switch into the B-like fit vs z."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
d = json.load(open("stage2_phase/POSTHOC_fixed_delta_all_switches.json"))
fig, ax = plt.subplots(figsize=(5.6, 3.8))
for key in ("0.05", "0.1", "0.2", "0.3", "0.4"):
    o = d[key]
    rows = [r for r in o["rows"] if "eps_c" in r]
    z = np.array([r["z"] for r in rows]); e = np.array([r["eps_c"] for r in rows])
    reg = np.array([r["P_cls"] == "regular" for r in rows])
    line, = ax.loglog(z, e, "-", lw=1, label=f"delta/pi={key}: alpha={o['fits']['full']['alpha']:.3f}")
    ax.loglog(z[reg], e[reg], "o", color=line.get_color(), ms=5)
    ax.loglog(z[~reg], e[~reg], "o", mfc="none", color=line.get_color(), ms=5)
zz = np.array([0.05, 0.8])
ax.loglog(zz, 0.3 * zz, "k:", lw=0.8, label="slope 1")
ax.loglog(zz, 0.3 * zz ** 2, "k--", lw=0.8, label="slope 2")
ax.set(xlabel="z", ylabel="epsilon_c of switch into B-like fit",
       title="POST_HOC_SENSITIVITY (filled: A regular; open: A near-coalescent)")
ax.title.set_fontsize(8)
ax.legend(fontsize=6)
fig.tight_layout()
fig.savefig("stage2_phase/POSTHOC_fixed_delta_scaling.png", dpi=130)
print("ok")
