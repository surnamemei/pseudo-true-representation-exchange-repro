# Stage 18: connected \(p\)-bridge certificate status

**Verdict: TSP-READY-WITH-DISCONNECTED-RESULTS.** The N21 numerical A/B crossing continues from sampled \(p=10^{-6}\) to \(p=4\), but the global bridge from the existential small-\(p\) theorem to the finite \(z=2\) certificate remains unproved. No connected theorem is added to the manuscript.

## Model and charts

This ledger concerns \(N=21,b=10,\phi=\pi\), \(p=z^2\), and the amplitude-eliminated two-tone least-squares loss \(J(u_1,u_2;p,\epsilon)\). The finite chart solves

\[
F=(\partial_{u_{A1}}J_A,\partial_{u_{A2}}J_A,
\partial_{u_{B1}}J_B,\partial_{u_{B2}}J_B,J_A-J_B)=0.
\]

The sampled small-\(p\) chart is \((v_A,\kappa_A,\kappa_B,v_B,\lambda)\). On \(p>0\), its exact map into the finite chart is

\[
u_{A1}=v_A,\quad u_{A2}=p\kappa_A,\quad
u_{B1}=p\kappa_B,\quad u_{B2}=v_B,\quad \epsilon=p\lambda.
\]

The inverse divides the two central coordinates and \(\epsilon\) by \(p\), so the maps are algebraic inverses on \(p>0\). The original small-spacing theorem uses a desingularized residual/Gram chart at \(p=0\). Stage 18 did **not** produce a parameterized interval chart-overlap certificate at an explicit positive \(p\), so this algebraic map does not join the proof components.

## What was computed

`p_bridge_numerical.csv` has 44 high-precision crossing solves over \(10^{-6}\le p\le4\); all predictor-corrector calls converged without a tangent-seed fallback. At the sampled points, both ordinary-coordinate Hessian minimum eigenvalues, both Gram determinants, both within-fit separations, A/B pair distance, and the ordinary \(\partial_\epsilon(J_A-J_B)\) slope are positive. The smallest sampled Hessian eigenvalues are \(2.75127\times10^{-15}\) for A and \(1.97033\times10^{-14}\) for B at the near-tangent endpoint. The crossing moves from \(\lambda_c(10^{-6})=0.065980410109242\ldots\) to \(\lambda_c(4)=0.0620476575430693\ldots\), with \(\epsilon_c(4)=0.248190630172277398\ldots\), matching the archived finite-width baseline. These are numerical statements between validation intervals.

The archived Stage-4 parametric Krawczyk boxes certify the **local** transverse crossing over \(1.9999\le z\le2\). The archived full-domain transfer certifies **global** A/B competition only over \(1.9999997\le z\le2\), equivalently \(3.99999880000009\le p\le4\). It excludes 14,287 outer and rechecks 819 inner cells with zero unresolved cells. Its smallest outer margin is \(4.9186206118\times10^{-7}\). Isolated global checks at \(z=1,1.5\) do not connect these intervals. The small-\(p\) theorem gives some \(p_0>0\) but no explicit lower endpoint for a finite-\(p\) slab, so a maximal interval \([0,p_{\max}]\) cannot be numerically stated. The largest **explicit** connected global interval currently verified is the near-\(p=4\) interval above.

## Bounded subdivision repair and obstruction

The single bounded repair in `stage18_slab_pilot.py` shrank the ordinary-coordinate interval Krawczyk slabs. It first passed a local crossing slab of width \(10^{-9}\) around \(z=0.1\), \(10^{-7}\) around \(z=0.5\), \(10^{-6}\) around \(z=1\), and \(10^{-5}\) around \(z=1.5,1.9\). At \(z=0.1\), width \(10^{-5}\) still gave contraction bound 1554; width \(10^{-9}\) gave 0.168. Even at \(z=1\), width \(10^{-5}\) gave contraction below one but failed box inclusion, while width \(10^{-6}\) passed. The finite chart is therefore a poor route to a complete small-\(p\) cover. These isolated local slabs are **not** joined and contain no new full-domain global exclusion. A practical bridge would need an interval implementation of the desingularized chart and recentered parameterized global bounds; Stage 18 stops after the prescribed bounded repair.

The obstruction is interval dependency and associated computational cost, with coordinate ill-conditioning near \(p=0\). This is **not** an established mathematical obstruction: the sampled branches show no coalescence, vanishing Gram determinant, Hessian zero at positive \(p\), or zero crossing slope; archived numerical third-candidate margins are positive. A sampled absence of events does not rule out an event between samples. The Stage-2 outer cell recorded in `stage4_bridge_transfer_limit.json` has a baseline lower-bound gap of only \(5.89789476\times10^{-6}\) despite a Gram condition number about 4.47 and distance 3.06 from the coalescent strip; uniform parameter transfer loses that small margin. This is an enclosure bottleneck rather than evidence of a third branch.

`p_bridge_slabs.csv` labels each row as local validated, local failed, or archived global validated. `p_bridge_competitor_margins.csv` distinguishes numerical candidate margins from the rigorous full-domain exclusion margin. `p_bridge_certificate.json` gives machine-readable scope and source hashes. Neither CSV may be read as a continuous global certificate.
