# Numerical invariance report (Stage 7)

Evidence level of every row: GLOBAL_NUMERICAL (search-protocol variations) unless stated.

Weakest relative third-margin qualifying crossing in Stages 1-4: {'third_margin_rel': 0.027642372131077647, 'eps_c': 0.009700815431107382, 'source': 'S2B_d0p05_z0p1.json'}

## Search-protocol variations (one factor at a time from production P)

| point | variant | verdict | eps_c | rel. change vs base | status |
|---|---|---|---|---|---|
| baseline | base | present | 0.2481906301722768 | None | present |
| baseline | grid_coarse | present | 0.24819063017227916 | 9.505692965487243e-15 | present |
| baseline | grid_fine | present | 0.24819063017227389 | 1.1742326604425418e-14 | present |
| baseline | rand_8 | present | 0.2481906301722768 | 0.0 | present |
| baseline | rand_128 | present | 0.2481906301722768 | 0.0 | present |
| baseline | kgrid_16 | present | 0.24819063017227777 | 3.914108868141806e-15 | present |
| baseline | kgrid_200 | present | 0.24819063017227747 | 2.68396036672581e-15 | present |
| baseline | tol_1e-7 | present | 0.24819063017227982 | 1.2189653332213054e-14 | present |
| baseline | tol_1e-11 | present | 0.2481906301722766 | 7.828217736283612e-16 | present |
| hann | base | present | 0.1471361491225018 | None | present |
| hann | grid_coarse | present | 0.14713614912249326 | 5.810072739158248e-14 | present |
| hann | grid_fine | present | 0.14713614912249068 | 7.564412884423563e-14 | present |
| hann | rand_8 | present | 0.1471361491225018 | 0.0 | present |
| hann | rand_128 | present | 0.14713614912249498 | 4.640512642314705e-14 | present |
| hann | kgrid_16 | present | 0.1471361491225018 | 0.0 | present |
| hann | kgrid_200 | present | 0.1471361491225018 | 0.0 | present |
| hann | tol_1e-7 | present | 0.14713614912248987 | 8.111465187785866e-14 | present |
| hann | tol_1e-11 | present | 0.1471361491225018 | 0.0 | present |
| kappa0p2_z0p5 | base | present | 0.03280025016219302 | None | present |
| kappa0p2_z0p5 | grid_coarse | present | 0.032800250162197474 | 1.3581511922257227e-13 | present |
| kappa0p2_z0p5 | grid_fine | present | 0.03280025016216711 | 7.899278118023129e-13 | present |
| kappa0p2_z0p5 | rand_8 | present | 0.03280025016219302 | 0.0 | present |
| kappa0p2_z0p5 | rand_128 | present | 0.03280025016219302 | 0.0 | present |
| kappa0p2_z0p5 | kgrid_16 | present | 0.03280025016219302 | 0.0 | present |
| kappa0p2_z0p5 | kgrid_200 | present | 0.03280025016219302 | 0.0 | present |
| kappa0p2_z0p5 | tol_1e-7 | present | 0.03280025016219302 | 0.0 | present |
| kappa0p2_z0p5 | tol_1e-11 | present | 0.032800250162164056 | 8.830098249766615e-13 | present |
| fixeddelta0p2_z0p2 | base | absent | 0.07329997538670686 | None | not_qualifying |
| fixeddelta0p2_z0p2 | grid_coarse | absent | 0.07329997538670686 | 0.0 | not_qualifying |
| fixeddelta0p2_z0p2 | grid_fine | absent | 0.07329997538670686 | 0.0 | not_qualifying |
| fixeddelta0p2_z0p2 | rand_8 | absent | 0.07329997538670686 | 0.0 | not_qualifying |
| fixeddelta0p2_z0p2 | rand_128 | absent | 0.07329997538670686 | 0.0 | not_qualifying |
| fixeddelta0p2_z0p2 | kgrid_16 | absent | 0.07329997538670686 | 0.0 | not_qualifying |
| fixeddelta0p2_z0p2 | kgrid_200 | absent | 0.07329997538670686 | 0.0 | not_qualifying |
| fixeddelta0p2_z0p2 | tol_1e-7 | absent | 0.07329997538669654 | 1.408605401917561e-13 | not_qualifying |
| fixeddelta0p2_z0p2 | tol_1e-11 | absent | 0.07329997538670686 | 0.0 | not_qualifying |
| weakest_margin | base | present | 0.009700815431107382 | None | present |
| weakest_margin | grid_coarse | present | 0.00970081543111876 | 1.172896377602046e-12 | present |
| weakest_margin | grid_fine | present | 0.009700815431117057 | 9.972927424739457e-13 | present |
| weakest_margin | rand_8 | present | 0.009700815431107382 | 0.0 | present |
| weakest_margin | rand_128 | present | 0.009700815431107382 | 0.0 | present |
| weakest_margin | kgrid_16 | present | 0.009700815431107382 | 0.0 | present |
| weakest_margin | kgrid_200 | present | 0.009700815431107382 | 0.0 | present |
| weakest_margin | tol_1e-7 | present | 0.009700815431174407 | 6.909162576827298e-12 | present |
| weakest_margin | tol_1e-11 | present | 0.009700815431102582 | 4.948016887978138e-13 | present |

Classification changes: none
Other flagged changes (eps_c rel > 1e-6 or status change): none

## Floating precision

{
 "baseline": {
  "float64": 0.2481906301722768,
  "mpmath40": "0.248190630172277398203609851355",
  "rel_diff": 2.4602970028319865e-15
 },
 "hann": {
  "float64": 0.1471361491225018,
  "mpmath40": "0.147136149122493769936895695746",
  "rel_diff": 5.470523033623325e-14
 },
 "kappa0p2_z0p5": {
  "float64": 0.03280025016219302,
  "mpmath40": "0.032800250162156159549674493471",
  "rel_diff": 1.1237537590515804e-12
 }
}

## Noisy near-crossing records (50 records, eta=0, 30 dB)

{
 "base": {
  "label_or_J_changes": 0,
  "labels": {
   "A": 31,
   "B": 19,
   "OUT": 0
  }
 },
 "grid_coarse": {
  "label_or_J_changes": 0,
  "labels": {
   "A": 31,
   "B": 19,
   "OUT": 0
  }
 },
 "grid_fine": {
  "label_or_J_changes": 0,
  "labels": {
   "A": 31,
   "B": 19,
   "OUT": 0
  }
 },
 "rand_8": {
  "label_or_J_changes": 0,
  "labels": {
   "A": 31,
   "B": 19,
   "OUT": 0
  }
 },
 "rand_128": {
  "label_or_J_changes": 0,
  "labels": {
   "A": 31,
   "B": 19,
   "OUT": 0
  }
 },
 "kgrid_16": {
  "label_or_J_changes": 0,
  "labels": {
   "A": 31,
   "B": 19,
   "OUT": 0
  }
 },
 "kgrid_200": {
  "label_or_J_changes": 0,
  "labels": {
   "A": 31,
   "B": 19,
   "OUT": 0
  }
 },
 "tol_1e-7": {
  "label_or_J_changes": 0,
  "labels": {
   "A": 31,
   "B": 19,
   "OUT": 0
  }
 },
 "tol_1e-11": {
  "label_or_J_changes": 0,
  "labels": {
   "A": 31,
   "B": 19,
   "OUT": 0
  }
 }
}

## Certificate arithmetic path

{
 "iv_arith": 1600,
 "iv_arith_fail": 0,
 "iv_trans": 561,
 "iv_trans_fail": 0,
 "arb_arith": 900,
 "arb_arith_fail": 0,
 "arb_trans": 576,
 "arb_trans_fail": 0,
 "arb_available": true
}

## Replay independence (static import audit)

{
 "stage12_full_arb_replay.py": {
  "imports": [
   "functools",
   "stage12_arb_backend"
  ],
  "imports_primary": [],
  "independent": true
 },
 "stage12_arb_backend.py": {
  "imports": [
   "csv",
   "decimal",
   "flint",
   "json",
   "math",
   "pathlib",
   "sys",
   "time"
  ],
  "imports_primary": [],
  "independent": true
 },
 "stage16_independent_arb.py": {
  "imports": [
   "__future__",
   "csv",
   "dataclasses",
   "flint",
   "fractions",
   "json",
   "math",
   "pathlib"
  ],
  "imports_primary": [],
  "independent": true
 },
 "bcert_arb.py (this run)": {
  "imports": [
   "__future__",
   "flint",
   "fractions",
   "math",
   "numpy",
   "pathlib",
   "sys"
  ],
  "imports_primary": [],
  "independent": true
 }
}

## Deterministic reproducibility

{
 "global_search_bitwise_identical": true,
 "n_minima": 77
}
