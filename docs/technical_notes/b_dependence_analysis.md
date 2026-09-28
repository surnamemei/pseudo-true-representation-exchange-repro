# Continuum critical coefficient versus weak-tone location

**Result.** The critical coefficient varies substantially with the normalized weak-tone location `b`. At the nine requested integer locations `b=6,...,14`, a high-precision stationary/equal-cost solve found two positive-curvature branches with nonzero satellite coefficients and a nonzero transverse slope. Numerical competitor scans support A/B globality at the sampled points, but **only b=10 has the existing directed full-line global certificate**. The strict b=10 margins and the implicit-function theorem nevertheless establish an unspecified open interval of b around 10 with a global noncoalescent crossing. No full certified b atlas is claimed.

## Definition and mechanism

At weak phase `phi=pi`, define `d(q)=2 sin(q/2)/q`, `mu2=1/12`, `mu4=1/80`, and

`q0=-mu2-lambda*d(b)`, `q1=2*lambda*d'(b)`,

`A=1-d(v)^2-d'(v)^2/mu2`,

`B=d''(v)-lambda*d(b-v)-d(v)*q0+d'(v)*q1/(2*mu2)`,

`R_inf(v,lambda;b)=mu4-2*lambda*d''(b)+lambda^2-q0^2-q1^2/(4*mu2)-B^2/A`.

For a fixed b, `lambda_inf(b)` is the positive solution of `R_v(v_A)=R_v(v_B)=0` and `R(v_A)=R(v_B)` when the two stationary points are the global competitors. The normalized small-spacing exponent is two because the aligned strong pair satisfies `2 cos(zs)=2-z^2*s^2+O(z^4)`. Its coefficient depends on b because `d(b)`, `d'(b)`, `d''(b)`, and `d(b-v)` change how the omitted tone projects onto the competing finite-window tangent subspaces. The b=6 numerical coefficient is 30.8% above b=10; the sampled coefficient has a shallow minimum near b=9, then rises again.

## Requested map

The high-precision root solve used 65 decimal digits. Displayed decimals below are rounded; full fields, including common cost, both curvatures, both nonzero satellite amplitudes, slope, nearest sampled regular competitor, and tail lower margin, are in [the CSV](lambda_inf_vs_b.csv).

| b | lambda_inf(b) | v_A | v_B | Coalescent cost gap | Nearest scanned regular gap | Assessment |
|---:|---:|---:|---:|---:|---:|---|
| 6 | 0.08696379 | -6.68490 | 10.94836 | 5.58e-4 | 6.46e-4 | A, numerical |
| 7 | 0.07496827 | -6.14923 | 11.28880 | 7.11e-4 | 9.93e-4 | A, numerical |
| 8 | 0.06881440 | -5.54274 | 11.64594 | 8.38e-4 | 1.49e-3 | A, numerical |
| 9 | 0.06638643 | -4.85083 | 12.02335 | 8.99e-4 | 2.19e-3 | A, numerical |
| 10 | 0.06650697 | -4.06003 | 12.42390 | 8.51e-4 | 3.16e-3 | A, **previously certified** |
| 11 | 0.06832370 | -3.16251 | 12.84923 | 6.70e-4 | 4.42e-3 | A, numerical |
| 12 | 0.07099941 | -2.16285 | 13.29923 | 3.84e-4 | 5.92e-3 | A, numerical |
| 13 | 0.07362000 | -1.08481 | 13.77217 | 1.11e-4 | 7.43e-3 | A, numerical |
| 14 | 0.07534750 | +0.02678 | 14.26582 | 7.07e-8 | 8.56e-3 | A at this point, **near central chart** |

The labels use the user's A/B/C/D/E classification: A means a sampled two-branch global crossing after the regular, coalescent, and tail checks described below. They are **numerical classifications** except the independent b=10 certificate. No sampled integer b had a lower third regular competitor (B), disappeared branch (C), lower coalescent competitor (D), or absent crossing (E). The b=14 A label must not be read as a connected noncoalescent theorem from b=10.

## Central-chart contact near b=14

Continuation of the A stationary root toward larger b gives `v_A=-0.0000710` at b=13.9758 and `+0.0000400` at b=13.9759. A 110-digit root solve and scalar bracketing locate its **numerical** central-chart contact at `b≈13.9758640010`. There the regular A coordinate reaches `v=0`, and its cost matches the confluent endpoint. This is a chart contact, not a collision of A and B. At b=14 the regular root exists on the other side, but its satellite coefficient is about `3.64e3` and its coalescent margin is only `7.07e-8`. Therefore a connected strict-margin/noncoalescent claim cannot be carried through the contact without a separate confluent analysis. No such theorem is attempted here.

## Numerical exclusion protocol and limits

For each fixed-b crossing, the script scans the regular scalar loss on `[-100,100]` with 8,001 points, refines detected local minima, and separately probes near the removable `v=0` endpoint. It compares the nearest other regular minimum and the exact limiting coalescent cost with the A/B cost. For `|v|>=100` and `6<=b<=14`, simple sinc derivative inequalities bound the projection gain: `|d(v)|<=2/|v|`, `|d'(v)|<=1/|v|+2/|v|^2`, `|d''(v)|<=1/(2|v|)+2/|v|^2+4/|v|^3`, and `|d(b-v)|<=2/(|v|-b)`. All nine resulting floating-evaluated tail lower margins are positive (minimum `7.15e-4` at b=6). The root solves are high precision; the frequency scan, local-minimum completeness, and tail-parameter substitution are **not directed interval certificates**. The b=10 full-line certificate remains the sole validated global anchor.

The [figure](lambda_inf_vs_b.pdf) connects only the well-separated sampled b=6 through 13 points. It marks certified b=10 separately and displays near-coalescent b=14 as an unconnected numerical point. Regenerate the CSV and figure with `python stage22_b_map.py`.
