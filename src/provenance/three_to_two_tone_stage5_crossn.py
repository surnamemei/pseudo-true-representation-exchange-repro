"""Stage 5 finite-N tangent laws and numerical full-record crossings.

The tangent results inherit Stage 3's local theorem. Full-record crossings
are numerical stationary equal-cost roots, not global certificates.
"""
import csv
from pathlib import Path
from mpmath import mp
from three_to_two_tone_stage3_tangent import solve
from three_to_two_tone_stage3_full_crossings import crossing

ROOT = Path(__file__).resolve().parent
NS = (11, 15, 21, 31, 41, 81, 161, 321)
ZS = ("0.15", "0.25", "0.4", "0.6", "1", "1.5", "2")

rows = []
for n in NS:
    tangent = solve(N=n)
    lam = mp.mpf(tangent["lambda_N"])
    quartic = mp.mpf(tangent["lambda_quartic"])
    print("tangent", n, mp.nstr(lam, 16), mp.nstr(quartic, 12), flush=True)
    base = dict(N=n, lambda_N=float(lam), c_N=float(quartic),
                tangent_v_A=float(tangent["v_A"]), tangent_v_B=float(tangent["v_B"]),
                tangent_Rvv_A=float(tangent["Rvv_A"]),
                tangent_Rvv_B=float(tangent["Rvv_B"]),
                tangent_transversality=float(tangent["delta_R_lambda"]),
                tangent_Gram_det_A=float(tangent["A_Gram_det"]),
                tangent_Gram_det_B=float(tangent["B_Gram_det"]))
    if n not in (11, 15, 21, 31, 41):
        rows.append(dict(base, z="", epsilon_cross="", epsilon_quadratic="",
                         epsilon_quartic="", rel_err_quartic="", u_A1="", u_A2="",
                         u_B1="", u_B2="", objective="", status="tangent_only"))
        continue
    seed = None
    for z in ZS:
        d = mp.mpf(z)
        tangent_seed = (mp.mpf(tangent["v_A"]), mp.mpf(tangent["A_coeff"][2])*d*d,
                        mp.mpf(tangent["B_coeff"][2])*d*d, mp.mpf(tangent["v_B"]),
                        lam*d*d)
        try:
            # A continuation Newton solve can jump to an unrelated equal-cost
            # stationary pair. First use the local theorem's labelled seed.
            root = crossing(z, N=n, seed=tangent_seed)
            e = mp.mpf(root["epsilon_cross"])
            pred = lam*d*d + quartic*d**4
            if not (mp.mpf("0.75")*pred < e < mp.mpf("1.25")*pred):
                raise ValueError("crossing root left the local A/B family")
            row = dict(base, z=float(d), epsilon_cross=float(e),
                       epsilon_quadratic=float(lam*d*d), epsilon_quartic=float(pred),
                       rel_err_quartic=float(abs(pred-e)/e),
                       u_A1=float(root["u_A1"]), u_A2=float(root["u_A2"]),
                       u_B1=float(root["u_B1"]), u_B2=float(root["u_B2"]),
                       objective=float(root["J_A"]), status="numerical_stationary_crossing")
            seed = tuple(mp.mpf(root[k]) for k in
                         ("u_A1", "u_A2", "u_B1", "u_B2", "epsilon_cross"))
            print("cross", n, z, mp.nstr(e, 15), flush=True)
        except Exception as exc:
            row = dict(base, z=float(d), epsilon_cross="", epsilon_quadratic=float(lam*d*d),
                       epsilon_quartic=float(lam*d*d+quartic*d**4),
                       rel_err_quartic="", u_A1="", u_A2="", u_B1="", u_B2="",
                       objective="", status=f"root_failed:{type(exc).__name__}")
            print("failure", n, z, str(exc)[:100], flush=True)
        rows.append(row)

with (ROOT / "crossN_scaling.csv").open("w", newline="") as file:
    writer = csv.DictWriter(file, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
