# Proof outline: affine rigidity and Gaussian exactness

Track I: informal mathematical proof. This outline accompanies the manuscript inspected on 31 October 2026; the original snapshot, final wording, complete diff and SHA-256 manifests are in `proof_digest_snapshot/`. The final wording was checked and changes no theorem or proof estimate. This outline does not replace the separate academic review.

## Statement

The paper identifies exactly when Gaussian couplings solve the unrestricted squared-distance Gromov–Wasserstein problem. It also identifies which affine maps can be optimal for a much broader class of source measures and quantifies how much Gaussian restriction can cost.

The objective is

\[
Q(\pi)=\mathbb E[(\|X-X'\|^2-\|Y-Y'\|^2)^2],
\]

where the two pairs are independent with law \(\pi\). Inner distances are squared; the outer loss is squared; there is no prefactor. These conventions are part of the theorem.

**Affine statement.** Let \(\mu\) have finite fourth moment. Assume its support contains a nonempty relatively open set in its affine hull \(E\). If \(T:E\to\mathbb R^m\) is affine, its graph coupling to \(T_\#\mu\) is optimal exactly when its restriction to \(E\) is a similarity or constant. In orthonormal coordinates on \(E\), this is \(L^\top L=\lambda I\), \(\lambda\ge0\).

**Gaussian statement.** For Gaussian marginals of positive covariance rank, a jointly Gaussian unrestricted optimiser exists exactly when the ranks agree and the ordered positive covariance eigenvalues are proportional. In that case every unrestricted optimiser is an affine similarity coupling. If either covariance is zero, the unique coupling is Gaussian; this exceptional case is handled separately.

**Quantitative statement.** Let \(a,b\) be decreasing covariance eigenvalue lists padded with zeros. Write \(p=\|a\|\|b\|>0\), \(c=a\cdot b\), \(S=p^2-c^2\), and \(D=Q_G-Q_*\). Then

\[
\frac{S}{p}\le D\le8(p-c).
\]

The upper bound is established prior work. The proposed contribution is the constructive lower bound and its relation to the exact classification. The theorem does not compute the unrestricted optimum or optimiser for nonproportional spectra.

## Idea

The objective depends on a coupling through a mixed fourth moment and its cross-covariance. An optimal affine coupling must also solve a linearised transport problem. Local cycles inside the source support then force a certain vector field to have zero curl. Its quadratic coefficients force the affine map to scale every direction equally. Gaussian marginals have interior support after removing zero-variance coordinates, so this obstruction classifies Gaussian exactness. A second route rearranges a latent standard Gaussian through an exactly measure-preserving nonlinear flow. The flow improves the mixed fourth moment while the covariance term has zero first variation. Combining all coordinate-pair perturbations gives a dimension-free finite improvement.

## Outline

1. **Reduce to moments.** Expand the pairwise-distance product to obtain \(Q=C_{\mu,\nu}-4G\), where \(G=\mathbb E[\|X\|^2\|Y\|^2]+2\|\mathbb E[XY^\top]\|_F^2\). **Location:** Lemma 1, Section 2. **Tag:** standard; the manuscript attributes this reduction to Delon et al. (2022).

2. **Compute the Gaussian benchmark.** Block covariance positivity gives a contraction representation for the cross-covariance; squared contraction entries are doubly substochastic. Rearrangement bounds its squared Frobenius norm by \(a\cdot b\), and Gaussian fourth moments yield \(Q_G=4(\sum a-\sum b)^2+8\|a-b\|^2\). **Location:** Lemma 2. **Tag:** standard; singular covariance supports are handled explicitly.

3. **Turn affine optimality into a curl obstruction.** First-order variation through mixtures gives a linear transport cost. Equal-mass cyclic rerouting proves finite-cycle inequalities on the source support. Small loops force the field \(g(x)=2\|x\|^2Mx+4MAx\) to have zero curl; its polynomial degree-two part forces \(M=L^\top L\) to be scalar. **Location:** Theorem 3. **Tag:** the point; the linearisation is established background, while the claimed affine rigidity is the candidate result.

4. **Complete the affine equivalence and demonstrate its boundary.** A similarity attains equality in the Cauchy–Schwarz bound for the two squared distances. The four-corner square projected onto a sign has an optimal anisotropic affine map because its support has no relative interior. **Location:** Theorem 3 sufficiency and the example following Corollary 4. **Tag:** routine, with the finite-support example serving a necessary scope check.

5. **Classify Gaussian exactness.** Orient the known Gaussian benchmark as an affine map from the larger covariance support. Affine rigidity gives equal ranks and proportional spectra. In the proportional case, equality of pairwise distances extends from almost everywhere to the support of every optimal plan; the support is a similarity graph. **Location:** Theorem 5. **Tag:** the point for the classification; the final support-extension argument uses standard isometry facts.

6. **Construct a feasible descent.** The polynomial pair field \(V_i=z_i(z_j^2-1),\ V_j=z_j(1-z_i^2)\) has bounded radial growth and zero Gaussian-weighted divergence. Its complete flow preserves the Gaussian marginal, and \(Q'(0)=-16(b_i a_j-b_j a_i)\). **Location:** Lemma 6. **Tag:** the point; the underlying two-variable Lotka–Volterra system is not claimed as new.

7. **Find the strongest first-order direction in the stated norm.** A polynomial gradient decomposition projects the mixed-moment derivative onto the Gaussian-divergence-free subspace. The residual is \(W_i=z_i\sum_j(b_i a_j-a_i b_j)(z_j^2-1)\), with squared norm \(4S\). The maximal objective descent slope at unit latent \(L^2\) velocity is \(8\sqrt S\). **Location:** Proposition 7. **Tag:** the point. The norm is explicitly a latent velocity norm, not an intrinsic global GW gradient norm.

8. **Convert the slope into a finite certificate.** Gaussian invariance makes composition a unitary group. Generator identities and a tensor-product correlation estimate bound \(|G''|\) by \(32\|a\|\|b\|S\). At time \(1/(8\|a\|\|b\|)\), the resulting objective improvement is at least \(S/(\|a\|\|b\|)\). A single pair gives a sharper certificate for the plane-to-line example. **Location:** Theorem 8, Section 6.1, and Section 7. **Tag:** the point.

## New components and reuse

- **Affine rigidity on interior support:** a structural obstruction for any finite-fourth-moment source, including some singular measures with dense support. Its relative-affine-hull formulation makes it applicable to lower-dimensional Gaussian supports. It is not a theorem about arbitrary finite clouds.
- **Explicit Gaussian-preserving descent:** a feasible coupling path that certifies strict improvement. The construction can be reused to test Gaussian restrictions for objectives containing the same mixed radial moment and covariance term.
- **Projected first-variation field:** the formula identifies a maximiser over a precisely stated weak divergence-free tangent space. The paper separately verifies that its polynomial maximiser generates an actual flow.
- **Dimension-free finite gap:** the combined field and tensor curvature estimate quantify the size of failure. The proof supplies an explicit feasible competitor rather than relying on a numerical optimiser.

These are proposed contribution units. Their historical novelty remains bounded by the literature search. The moment reduction, Gaussian benchmark, proportional-spectrum sufficiency, general transport linearisation and upper gap bound are credited to prior work.

## Prior work

The inspected Salmona–Delon–Desolneux preprint gives the Gaussian-restricted value, the proportional exact case, and the upper-gap expression, and conjectures wider exactness. Peyré's author post of 23 October 2026 gives a non-Gaussian radial competitor in sufficiently large dimension. Dumont, Lacombe and Vialard provide the relevant squared-distance linearisation and analyse optimal map or two-map structure. These sources establish the problem boundary; they do not by themselves verify historical priority for the proposed necessity theorem.

The current book PDF linked by Peyré was subsequently acquired and its Gaussian section on printed pp. 302–304 independently read. Propositions 13.33–13.34 give the known Gaussian benchmark and the large-dimensional radial counterexample, followed by a discussion of proportional-spectrum sufficiency. No necessity classification, affine rigidity theorem or positive lower-gap bound appears in that inspected passage. The PDF is identified by SHA-256 `b47340a8f0a9ba8d9b838d9c412fbf67e2d6804c6c5b788cc444b4d675b8528a` (27,846,433 bytes). The whole book was not read page by page, and no absence-of-results search can certify priority.

## Verification status

No proof assistant was used. This reviewer independently reconstructed the general curl calculation, Gaussian first variation, pair generator norm, tensor curvature improvement, combined-flow norm identities, all-optimiser support argument and finite-square counterexample. The final exact checker reproduced all 441 integer pair cases and the \(76/7\) rational certificate. A separate exact symbolic checker reproduced 28 identities. Independent direct quadrature of generator polynomials checked 120 cases in dimensions two through four, with maximum scaled residual below \(5.6\times10^{-16}\). The finite-square calculation gave optimum eight, and an independently implemented two-ODE reduction matched the full flow in 60 cases, including zero coordinates and proportional spectra, with maximum scaled reconstruction discrepancy below \(4.3\times10^{-12}\).

The finite checks support implementation consistency. The general theorem rests on the analytic arguments. No unresolved mathematical gap was found within this inspection. Novelty, practical impact and any eventual journal judgement require separate assessment.

## References

Delon, J., Desolneux, A., & Salmona, A. (2022). Gromov–Wasserstein distances between Gaussian distributions. *Journal of Applied Probability, 59*(4), 1178–1198. https://doi.org/10.1017/jpr.2022.16. Inspected preprint: https://arxiv.org/pdf/2104.07970

Dumont, T., Lacombe, T., & Vialard, F.-X. (2024). *On the existence of Monge maps for the Gromov–Wasserstein problem* (Version 3) [Preprint]. https://arxiv.org/pdf/2210.11945

Peyré, G. (2026, October 23). *Gaussians are not closed for Gromov–Wasserstein*. https://www.gpeyre.com/blog/2026/10/23/gaussian-gromov-wasserstein-counterexample/

Peyré, G. (2026). *Optimal transport for machine learners* [Online book; cover dated 18 October 2026, file retrieved 31 October 2026]. https://www.gpeyre.com/ot4ml/OT4ML/OT4ML.pdf
