# Finite-width global epsilon certificate at N=21, z=2

**Certified closed interval:** \([0.213190630172278,\,0.283190630172277]\), width \(0.069999999999999\). It is contained strictly inside both the outer and root/inner working intervals, each approximately \(\epsilon_c\pm0.035\). The unique crossing lies in the outward-rounded bracket \([0.248190630172276,\,0.248190630172279]\); the exact directed endpoint-sign bracket is in the JSON certificate.

For every epsilon in the displayed interval, A and B are separated strict local minima with positive Gram determinants and finite fitted amplitudes. They are the only globally competitive two-tone fits. The optimized cost difference has positive derivative (uniform lower bound \(>6.844\)), so A is uniquely global below the single crossing, B above it, and both tie there. No third regular or coalescent fit is competitive.

## Exact quadratic structure

For a fixed fitted subspace \(P\), write \(x(\epsilon)=x_0+\epsilon x_1\). Then
\[
J(P;\epsilon)=\|(I-P)x(\epsilon)\|^2
=a(P)\epsilon^2+b(P)\epsilon+c(P),
\]
where \(a=\|(I-P)x_1\|^2\), \(b=2\operatorname{Re}\langle(I-P)x_0,(I-P)x_1\rangle\), and \(c=\|(I-P)x_0\|^2\). The code encloses these coefficients with directed intervals. For each outer spatial cell, the centered confluent basis gives a quadratic \(J_c(\epsilon)\), while fixed A/B fitted subspaces give feasible quadratic incumbents \(U_A,U_B\). The spatial projector displacement \(\theta\) is independent of epsilon. On each epsilon interval, the code bounds the minimum of the quadratic difference \(J_c-U_w\) and the maximum of \(2\theta\sqrt{J_c E}-\theta^2E\), after verifying \(\sqrt{J_c}-\theta\sqrt E>0\). Their strict inequality excludes the entire cell.

Inside the two accepted neighborhoods, Taylor lower bounds, gradients, and Krawczyk exclusion are also evaluated as coefficient-wise epsilon quadratics. Epsilon is subdivided only where a spatial cell needs it. Each root is enclosed by a parameterized Krawczyk map on 4,096 adjacent epsilon slabs; 8,190 endpoint checks certify that neighboring slabs share the same root. The Gram determinant stays above \(420.905\), and the slope stays above \(6.844\).

## Coverage accounting

The 14,287 archived outer exclusion parents become 19,054 terminal leaves after refinement of 3,738 parents. The 819 archived crossing inner parents generate 985,705 terminal three-dimensional leaves from 1,970,591 visited nodes. All leaves pass; none is unresolved. The exact-decimal audit checks containment and volume identity for every archived parent, and checks that all 4,096 epsilon slabs meet exactly. The original directed bracket has opposite certified endpoint signs. The per-parent and per-slab record is [epsilon_interval_obligations_0p035.csv](epsilon_interval_obligations_0p035.csv); the full inner leaves and root checks are in the Stage-17 machine logs.

## Limit

This certifies a finite epsilon interval for \(N=21,z=2\). It does not prove a connected \(z\)-continuation from the small-spacing theorem. The separately certified outer amplitude samples are distinct records; the displayed inward-rounded closed interval is the continuous claim. The endpoints are sufficient, not maximal.
