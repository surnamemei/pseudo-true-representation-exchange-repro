"""A-priori size of the leading term neglected by the first-order theory (diagnostic only; not used in any prediction
or classification).

Re-optimizing branch j after adding noise lowers its cost by q_j(w) = b_j^T H_j^{-1} b_j, b_j = Re(M_j^H w), where M_j is
the model Jacobian (6 real parameters) and H_j the Hessian of s6_noise.local_cov. Hence E q_j = (sigma^2/2) tr(H_j^{-1} G_j),
G_j = Re(M_j^H M_j). The mean noisy gap is Delta0 - E[q_A - q_B], so the first-order theory omits a crossing shift
    delta_eta = E[q_A - q_B] / (g eps_c),
reported here per SNR relative to the predicted width. Reads only deterministic inputs.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "validation/adversarial_overnight"))

import numpy as np  # noqa: E402

import s6_noise  # noqa: E402
from advcore import make_signal, times  # noqa: E402

DET = ROOT / "results/transition_theory/deterministic"
N, Z, BW, EC = s6_noise.N, s6_noise.Z, s6_noise.BW, s6_noise.EC


def trace_HinvG(x, u):
    """Same Jacobian and Hessian construction as s6_noise.local_cov."""
    s = times(N)
    v = np.exp(1j * np.outer(s, u))
    a = np.linalg.lstsq(v, x, rcond=None)[0]
    r = x - v @ a
    jac = np.column_stack([v[:, 0], 1j * v[:, 0], v[:, 1], 1j * v[:, 1], 1j * s * a[0] * v[:, 0], 1j * s * a[1] * v[:, 1]])
    second = np.zeros((N, 6, 6), complex)
    for k in range(2):
        second[:, 2 * k, 4 + k] = second[:, 4 + k, 2 * k] = 1j * s * v[:, k]
        second[:, 2 * k + 1, 4 + k] = second[:, 4 + k, 2 * k + 1] = -s * v[:, k]
        second[:, 4 + k, 4 + k] = -s * s * a[k] * v[:, k]
    G = (jac.conj().T @ jac).real
    H = G - np.einsum("n,nij->ij", np.conj(r), second).real
    return float(np.trace(np.linalg.solve(H, G)))


def main():
    det = json.loads((DET / "deterministic.json").read_text())
    x = make_signal(N, Z, BW, EC)
    tA, tB = trace_HinvG(x, np.array(det["u_A_c"])), trace_HinvG(x, np.array(det["u_B_c"]))
    out = dict(note="diagnostic only; not used in any prediction or classification", tr_HinvG_A=tA, tr_HinvG_B=tB,
               correctly_specified_value=6.0, per_snr={})
    for snr, p in det["per_snr"].items():
        s2 = p["sigma2_c"]
        dq = s2 / 2 * (tA - tB)
        shift = dq / (det["g"] * EC)
        out["per_snr"][snr] = dict(E_qA_minus_qB=dq, crossing_shift_eta=shift, shift_over_W=shift / p["W_closed_form"])
    (DET / "second_order_diagnostic.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
