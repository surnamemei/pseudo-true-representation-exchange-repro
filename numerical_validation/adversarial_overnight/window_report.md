# Stage 1 — window / metric invariance

Evidence level: GLOBAL_NUMERICAL (pre-registered crossing analysis CA; production search at 101 scan points, reference search at each crossing). Weights normalised to Σw = N; J_w = Σ w_t |x_t − (VB)_t|².

| window | role | verdict | ε_c | J* | A pair | B pair | third margin (rel.) | ε_c·slope/J* | coalescence |
|---|---|---|---|---|---|---|---|---|---|
| rect | primary | PERSISTS | 0.2481906302 | 0.811808 | (-4.7724, 0.6052) | (-0.0577, 12.4634) | 0.883 (1.088) | 2.53 | no |
| hann | primary | PERSISTS | 0.1471361491 | 0.130183 | (-7.1789, 0.0142) | (-0.2093, 15.0343) | 0.09449 (0.726) | 3.19 | no |
| hamming | primary | PERSISTS | 0.1835974128 | 0.337679 | (-5.3941, 0.2437) | (-0.1486, 13.3260) | 0.265 (0.785) | 2.71 | no |
| dpss | primary | PERSISTS | 0.1556046466 | 0.171524 | (-6.2951, 0.0774) | (-0.1968, 14.3552) | 0.1184 (0.690) | 2.95 | no |
| hann2 | secondary | PERSISTS | 0.1215600636 | 0.0433867 | (-8.9305, -0.1381) | (-0.2655, 16.9431) | 0.02671 (0.616) | 3.66 | no |
| homotopy:0.0 | secondary | PERSISTS | 0.2481906302 | 0.811808 | (-4.7724, 0.6052) | (-0.0577, 12.4634) | 0.883 (1.088) | 2.53 | no |
| homotopy:0.1 | secondary | PERSISTS | 0.2444495255 | 0.791898 | (-4.8452, 0.5697) | (-0.0577, 12.5236) | 0.8533 (1.078) | 2.54 | no |
| homotopy:0.2 | secondary | PERSISTS | 0.2401222164 | 0.763518 | (-4.9155, 0.5346) | (-0.0599, 12.5825) | 0.8137 (1.066) | 2.56 | no |
| homotopy:0.3 | secondary | PERSISTS | 0.2351081137 | 0.726572 | (-4.9824, 0.4994) | (-0.0643, 12.6417) | 0.7613 (1.048) | 2.58 | no |
| homotopy:0.4 | secondary | PERSISTS | 0.2292597767 | 0.680699 | (-5.0450, 0.4635) | (-0.0711, 12.7037) | 0.6966 (1.023) | 2.59 | no |
| homotopy:0.5 | secondary | PERSISTS | 0.2223595686 | 0.625213 | (-5.1031, 0.4260) | (-0.0809, 12.7728) | 0.6202 (0.992) | 2.61 | no |
| homotopy:0.6 | secondary | PERSISTS | 0.2140772399 | 0.559007 | (-5.1585, 0.3851) | (-0.0941, 12.8575) | 0.5317 (0.951) | 2.63 | no |
| homotopy:0.7 | secondary | PERSISTS | 0.2038866944 | 0.480361 | (-5.2189, 0.3381) | (-0.1118, 12.9760) | 0.4307 (0.897) | 2.65 | no |
| homotopy:0.8 | secondary | PERSISTS | 0.190882406 | 0.386578 | (-5.3135, 0.2786) | (-0.1353, 13.1754) | 0.3184 (0.824) | 2.68 | no |
| homotopy:0.9 | secondary | PERSISTS | 0.1732946486 | 0.273109 | (-5.5795, 0.1901) | (-0.1671, 13.6130) | 0.2012 (0.737) | 2.76 | no |
| homotopy:1.0 | secondary | PERSISTS | 0.1471361491 | 0.130183 | (-7.1789, 0.0142) | (-0.2093, 15.0343) | 0.09449 (0.726) | 3.19 | no |

Rectangular→Hann homotopy (secondary, pre-registered): local continuation in 101 steps of θ, complete=True; at every θ ∈ {0,0.1,…,1} the independently searched crossing coincides with the continued A/B pair: True.

Classification: every primary nonrectangular window (Hann, Hamming, DPSS) → **PERSISTS**. The crossing amplitude is window dependent (ε_c from 0.1471 (Hann) to 0.2482 (rectangular)), but the separated two-branch exchange with a large third-competitor margin is not a rectangular-window artifact.

Plots: `stage1_window/window_branches.png|pdf`, `stage1_window/window_homotopy.png`. Arrays: `stage1_window/window_minima.npz`. Table: `stage1_window/window_summary.csv`.