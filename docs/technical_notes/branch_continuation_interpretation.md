# Interpretation of amplitude continuation

N=21, z=2, b=10, phi=pi. Deterministic stationary continuation uses amplitude steps no larger than 0.0005, analytic gradients/Hessians, and 65-digit independent stationary refinement at the displayed landmarks. Complex amplitudes are optimized at each step; conjugate symmetry makes them real in this centered phase-pi family. The CSV still records real and imaginary parts explicitly.

## Observed branches

|Branch|epsilon|u1|u2|Loss|
|---|---:|---:|---:|---:|
|A|0.0000000000|-2.0000000000|2.0000000000|0.0000000000|
|A|0.0500000000|-2.6117933963|1.4331936394|0.0256680114|
|A|0.1000000000|-3.2315831902|1.0775903671|0.1120222739|
|A|0.2481906302|-4.7724473539|0.6051703203|0.8118078719|
|B|0.0000000000|-0.0200738086|14.6313005826|1.4420915912|
|B|0.0500000000|-0.0197411651|14.1099198513|1.2287423409|
|B|0.1000000000|-0.0249810774|13.6136562750|1.0709769042|
|B|0.2481906302|-0.0576535619|12.4634193568|0.8118078719|

A has 602 samples and B 601. Both were followed over the physical range 0<=epsilon<=approximately 0.30, across epsilon_c=0.248190630172277. Neither encountered a sampled fold, singularity, coalescence or loss of positive curvature. Minimum sampled Hessian eigenvalues: A=0.0489333571, B=0.0670904357. These trajectory checks are numerical; the crossing itself has the independent interval certificate.

A begins exactly at the generating strong pair (-2,2), with unit amplitudes and zero loss. Its increasingly asymmetric frequencies and amplitudes provide a continuous deformation of that representation. B already exists at zero omitted amplitude as a worse local approximation. It reaches the physical boundary epsilon=0 without termination, so there is no positive-amplitude birth event to identify on this backward segment.

At equality, B combines a nearly central component with a satellite at 12.4634193568, **not** at b=10. The satellite is displaced by 2.4634193568 normalized radians. Fixing the satellite to b and optimizing the other frequency produces u1=0.0366746458, loss=1.4247760600, versus 0.8118078719 for the unconstrained branches (penalty 0.6129681881). A 4096-interval full-period scan and local refinements found this as the best constrained minimum; this comparison is numerical, not an additional global certificate. The established tangent fixed-b loss 0.1036226553 also exceeds the equal tangent loss 0.0587626233.

## Signal-processing interpretation

Immediately below the certified crossing, the winning representation is a deformation of the exact strong-pair fit. Immediately above it, a central component represents the close pair collectively and a biased satellite accounts for the omitted component's influence. B should not be described as literally estimating the weak tone. Its movement depends on the finite-record projection geometry and the simultaneous fit of the central component.

The continuation figure plots A1/A2 and B1/B2 with the generating frequencies and crossing line. It establishes numerical branch identity from zero amplitude; it does not upgrade the entire displayed epsilon interval to an interval-certified global ordering. The separately certified below/above points and the narrow crossing bracket retain their original scope.

Outputs: `branch_continuation.csv`, `branch_continuation_figure.pdf`, `stage12_continuation_checks.json`.
