# Stage 0 — baseline reproduction

## Numerical reproduction (N=21, z=2, b=10, φ=π, rectangular)

| check | criterion (pre-registered) | result | pass |
|---|---|---|---|
| B0.1 float64 crossing | |Δε_c| ≤ 1e-10 | ε_c = 0.2481906301722716, |Δ| = 5.80e-15 | True |
| B0.2 50-digit joint solve | |Δ| ≤ 1e-30 vs 110-digit reference | ε_c = 0.248190630172277398203609851355192030880784691; max root/ε diff 5.1288e-51; J* diff 1.3364e-51 | True |
| B0.3 Hessians | both PD | A eig(u) = 0.054625, 5.2479; B eig(u) = 0.270032, 8.7393 | True |
| B0.4 reference global search | A,B lowest two; third margin 0.88296 ± 1e-5 | third minimum at (0.06713, 20.49489), margin 0.8829599621 | True |
| B0.5 ordering | A below / B above; also at certified interval ends | ε=0.247191: J_A=0.804648, J_B=0.812892; ε=0.249191: J_A=0.819003, J_B=0.810729; ε=0.213191: J_A=0.581623, J_B=0.853747; ε=0.283191: J_A=1.084376, J_B=0.777411 | True |

Overall numeric pass: **True** (evidence: GLOBAL_NUMERICAL reproduction of an interval-certified result). Branch A at the crossing: (-4.772447353919, 0.605170320277); B: (-0.057653561927, 12.463419356802); J* = 0.8118078719175; dΔJ/dε = 8.259051.

## Theorem-level certificate replays (sandbox copies; original scripts unmodified; cwd = sandbox)

| chain | commands (exit code, seconds) | outputs identical to frozen | substantive JSON differences | status |
|---|---|---|---|---|
| R1_arb_finite_z | `stage12_full_arb_replay.py 21` (0, 66s); `stage12_full_arb_replay.py 31` (0, 52s); `stage12_cover_audit.py` (0, 1s) | 6/8 (non-identical: stage12_replay/N21_summary.json, stage12_replay/N31_summary.json — timing fields or key order only) | 0 | CERTIFIED (replayed) |
| R2_generalN | `stage12_generalN.py 11` (0, 148s); `stage12_generalN.py 15` (0, 227s); `stage12_generalN.py 21` (0, 341s); `stage12_generalN.py 31` (0, 545s); `stage12_generalN.py 41` (0, 828s) | 20/20 | 0 | CERTIFIED (replayed) |
| R3_continuum | `stage14_continuum.py` (0, 5s); `stage14_continuum_interval.py` (0, 0s); `stage14_continuum_global_hybrid.py` (0, 98s); `stage14_replay_global.py` (0, 51s); `stage14_continuum_endpoint_interval.py` (0, 0s) | 7/8 (non-identical: continuum_global_certificate.json — timing fields or key order only) | 0 | CERTIFIED (replayed) |
| R4_independent_arb_N0 | `stage16_independent_arb.py` (0, 45s) | 1/2 (non-identical: independent_N0_replay.json — timing fields or key order only) | 0 | CERTIFIED (replayed) |
| R5_largeN_floor_10001 | `stage17_uniform_trial.py 10001` (0, 337s); `stage17_refine_uniform.py 10001` (0, 2s); `stage17_largeN_margin_summary.py` (0, 358s); `stage17_write_largeN.py` (0, 0s); `stage17_annotate_condition_thresholds.py` (0, 0s) | 5/5 | 0 | CERTIFIED (replayed) |
| R6_finite_width_eps_0p035 | `stage17_epsilon_refine.py 0.035 14287` (0, 162s); `stage17_epsilon_root_tiling.py 0.035 4096` (0, 1354s); `stage17_epsilon_inner_tiled.py _0p035 3000000` (0, 1824s); `stage17_epsilon_joints.py _0p035` (0, 1030s); `stage17_epsilon_coverage_audit.py 0.035` (0, 11s); `stage17_write_epsilon_certificate.py 0.035` (0, 8s); `stage17_epsilon_coeff_crosscheck.py` (0, 2s) | 9/9 | 0 | CERTIFIED (replayed) |
| R7_lower_transition | `stage13_build.py` (0, 10s) | 3/3 | 0 | CERTIFIED (replayed) |

Technical note: the README lists `stage17_epsilon_refine.py 0.035 40`, but the downstream coverage audit and certificate writer read the *complete* outer refinement `stage17_epsilon_outer_refinement_0p035_14287.json`; the replay therefore ran `stage17_epsilon_refine.py 0.035 14287`, which regenerates exactly that frozen file (byte-identical). `stage16_independent_arb.py` was run with the workspace-local python-flint 0.9.0 (`stage12_deps`) on `PYTHONPATH` instead of the WSL venv named in CODE_README; both are python-flint 0.9.0.

Integrity: `integrity_check.json` re-hashes all 766 frozen files after the replays — none changed.

**Stage 0 verdict: baseline reproduces; every theorem-level certificate chain replays to identical outputs. No NOT_READY trigger.**