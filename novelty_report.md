# Novelty-gate report: Gaussian exactness and affine rigidity for squared-distance Gromov–Wasserstein transport

**Search date:** 31 October 2026.  
**Researcher/agent:** AI-assisted research session, with separate scouting and proof-audit passes. This is not external human peer review.  
**Corpus boundary:** Two live web-search indexes; openly indexed primary papers on arXiv, publisher and proceedings websites; author webpages; the downloaded current OT4ML book, public repository and chapter audit; an OEIS fingerprint check. Formula passages were opened when accessible. No subscription citation database or author correspondence was used.  
**Intended contribution units:** Classification theorem; general affine rigidity theorem; explicit perturbative proof; quantitative approximation-gap lower bounds; exact first-variation sensitivity; reproducible verification artefacts.

## 1. Frozen claims

### Exact statements being assessed

The objective throughout is

\[
 Q(\pi)=\iint (\|x-x'\|^2-\|y-y'\|^2)^2\,
 d\pi(x,y)d\pi(x',y'),
\]

with no prefactor, and independent pairs distributed according to the same coupling. Write \(Q_*\) for its infimum over all couplings and \(Q_G\) for its infimum over jointly Gaussian couplings.

1. **Gaussian exactness.** For finite-dimensional Gaussian marginals of positive covariance ranks, a jointly Gaussian unrestricted optimum exists if and only if their decreasing, zero-padded covariance spectra \(a,b\) are proportional. Equivalently, the positive ranks agree and their positive eigenvalues are proportional. A zero-rank marginal is a separate trivial exception.
2. **Affine rigidity.** For \(\mu\in\mathcal P_4(\mathbb R^d)\) whose support contains a nonempty open set, an affine graph coupling \((\mathrm{Id},Lx+c)_\#\mu\) is globally optimal between its marginals if and only if \(L^TL=\lambda I_d\), \(\lambda\ge0\). The statement can be made in the affine hull using relative interior and the restriction of \(L\) to its direction space. A positive \(C^1\) density on an open set gives the stronger local failure result against compactly supported \(C^1\) measure-preserving flows.
3. **Dimension-free gap.** Put \(A=\|a\|_2\), \(B=\|b\|_2\), \(c=a\cdot b\), and \(S=A^2B^2-c^2\). For positive ranks, the proposed new lower bound is

   \[
   \frac{S}{AB}\le Q_G-Q_*.
   \]

   The complementary upper bound \(Q_G-Q_*\le8(AB-c)\) is established background, identified precisely below. Consequently the lower estimate makes the actual Gaussian approximation gap comparable, with universal constants, to an already known spectral collinearity defect.
4. **Constructive pair certificates.** For \(\Delta_{ij}=b_i a_j-b_j a_i\), set

   \[
   P_{ij}(a)=24(a_i^2+a_j^2)-32a_i a_j,
   \quad K_{ij}=\sqrt{P_{ij}(a)P_{ij}(b)}+8(a_i b_i+a_j b_j).
   \]

   Whenever \(\Delta_{ij}\ne0\), the explicit feasible flow gives

   \[
   Q_G-Q_*\ge\frac{32\Delta_{ij}^2}{K_{ij}}>0.
   \]

5. **Specified latent sensitivity.** At the canonical Gaussian coupling, the steepest first-order decrease of \(Q\) over unit \(L^2(\gamma_d)\) Gaussian-divergence-free latent velocities is \(8\sqrt S\), with an explicit attaining complete field. This is a statement about the specified latent velocity norm, not an assertion about an intrinsic global GW gradient norm.

### Exclusions

- No claim to have computed the unrestricted optimal coupling or exact unrestricted value for nonproportional spectra.
- No novelty claim for the Gaussian restricted formula, moment expansion, proportional-spectrum sufficiency, old upper-gap bound, Gaussian-preserving-flow method, weighted Helmholtz decomposition, Lagrange's identity, or the Lotka–Volterra vector field.
- The ordinary-distance, inner-product, entropically regularised, fused, semi-relaxed and unbalanced problems have different objectives or constraints.
- A source with finitely many support points does not satisfy the interior-support affine theorem. The local flow theorem requires \(\mu|_U=\rho\,dx\), with no singular part in the active open set.
- The known large-dimensional non-Gaussian counterexample is credited to Peyré. Merely disproving universal Gaussian optimality is not claimed as new.

### Permitted language

“The classification, affine rigidity, lower-gap estimate and exact sensitivity were not found in the bounded search documented here.” “The new necessity argument complements known proportional-spectrum sufficiency.” “The construction strengthens the previously reported large-dimensional failure to every nonproportional pair.”

Do not write “first ever”, “previously unknown in all literature”, “the most important open problem in applied mathematics is solved”, or “the general Gaussian GW problem is completely solved”.

## 2. Normalisation and verified translations

| Source/object | Canonical notation | Transformation checked |
|---|---|---|
| Delon–Desolneux–Salmona's \(GGW_2^2\) and \(GW_2^2\), squared Euclidean ground dissimilarities | \(Q_G\) and \(Q_*\) | Their defining integral is the present integral, with the same absence of a prefactor. |
| Their diagonal covariance matrices \(D_0,D_1\), with the shorter list extended by zeros | Decreasing nonnegative \(a,b\) | Orthogonal changes of coordinates and translations preserve \(Q\); null covariance coordinates can be discarded. |
| Their \(\|\alpha\|\|\beta\|-\alpha^T\beta\), Eq. (5.4) | \(AB-c\) | Same ordered padded eigenvalue lists; this is a known difference between their Gaussian upper bound and lower bound. |
| General centred moment reduction | \(G(\pi)=E[\|X\|^2\|Y\|^2]+2\|E[XY^T]\|_F^2\) | \(Q\) equals a marginal constant minus \(4G\). Thus improvement means increasing \(G\). |
| Linearised GW transport cost in Dumont–Lacombe–Vialard | \(-\|x\|^2\|y\|^2-4x^TKy\) | Same cost after removing marginal-only terms and an irrelevant positive scale. |
| Inner-product GW | Objective depending only on \(\|E[XY^T]\|_F^2\), up to marginal terms | It lacks the mixed radius moment. Gaussian optimality for this objective does not establish Gaussian optimality here. |
| Pair Gaussian field | \(V_i=z_i(z_j^2-1),\ V_j=z_j(1-z_i^2)\) | On \(p=z_i^2,q=z_j^2\): \(p'=2p(q-1),q'=2q(1-p)\), the classical Lotka–Volterra centre. |
| Simultaneous field | \(W_i=z_i\sum_j C_{ij}(z_j^2-1),\ C=ba^T-ab^T\) | On squared coordinates: \(p_i'=2p_i\sum_j C_{ij}(p_j-1)\), a skew-symmetric Lotka–Volterra system. |
| Pair coefficients | \(\Delta_{ij}=C_{ij}\), \(i<j\) | Reversing the pair reverses the field and the sign of its derivative. |
| Exterior product | \(S=\sum_{i<j}\Delta_{ij}^2=A^2B^2-c^2\) | Lagrange's identity, not a new spectral invariant. |

## 3. Formula fingerprints

### Small exact instances

| \(a\) | \(b\) | \(Q_G\) | Improving ordered pair | \(Q'(0)\) | \(S\) |
|---|---|---:|---|---:|---:|
| \((1,1)\) | \((1,0)\) | 12 | \((1,2)\) | \(-16\) | 1 |
| \((2,1)\) | \((1,1)\) | 12 | \((2,1)\) | \(-16\) | 1 |
| \((2,1)\) | \((4,2)\) | 76 | none | 0 | 0 |

For standard Gaussian rank two versus rank one, \(K=8+8\sqrt6\), so the pair certificate gives \(Q_*\le12-4/(1+\sqrt6)\). Time \(1/7\) also certifies \(Q_*<76/7\). These are upper bounds on \(Q_*\), not its exact value.

The finite arrays in `target.yaml` are pairs of spectral lists, not a new family of matrices. Characteristic-polynomial or graph-isomorphism fingerprints would therefore be artificial for this target. Their relevant invariant is the two-by-two determinant.

### Sequence and generating function

For standard Gaussian ranks \(n\) and one, the known Gaussian value is \(q_n=4(n^2-1)\), \(n\ge1\):

`0, 12, 32, 60, 96, 140, 192, 252, 320, 396`.

Its division by four is OEIS A005563 with index shifted by one: `0, 3, 8, 15, 24, 35, 48, 63, 80, 99`. The first differences are `12, 20, 28, 36, 44, ...`; the second difference is eight. The canonical generating functions are

\[
\sum_{n\ge1}q_n z^n=\frac{4z^2(3-z)}{(1-z)^3},
\qquad
\sum_{n\ge1}q_n z^{n-1}=\frac{4z(3-z)}{(1-z)^3}.
\]

This elementary sequence is an object-identity check for the already known benchmark, not a proposed contribution.

### Differential and invariant fingerprints

For the pair field,

\[
\operatorname{div}V-z\cdot V=0,
\quad G'(0)=4\Delta_{ij},\quad Q'(0)=-16\Delta_{ij},
\quad p+q-\log p-\log q=\text{constant}
\]

in a positive quadrant. In Hermite notation its active components are \(H_1(z_i)H_2(z_j)\) and \(-H_2(z_i)H_1(z_j)\). For the simultaneous field, \(\sum_i(p_i-\log p_i)\) is conserved in the positive orthant.

For the sensitivity, with \(H_a=a\cdot z^2\), \(H_b=b\cdot z^2\),

\[
2H_aD_bz=\nabla\Phi+W,
\quad \Phi=\tfrac12[H_aH_b+(\sum a)H_b-(\sum b)H_a],
\quad \|W\|_{L^2(\gamma)}^2=4S.
\]

The dimension-free lower estimate uses

\[
G'(0)=4S,\quad
\|\mathcal L H_a\|_2^2
=8\|Ca\|_2^2+16\sum_{i<j}C_{ij}^2(a_i-a_j)^2
\le24A^2S,
\]

and the uniform curvature bound \(|G''(t)|\le32ABS\). The tested time is \(t=1/(8AB)\).

For affine rigidity, the exact polynomial fingerprint is

\[
g(x)=2\|x\|^2Mx+4MAx,
\quad
\partial_jg_i-\partial_ig_j
=4[x_j(Mx)_i-x_i(Mx)_j]+4(MA-AM)_{ij},
\]

where here \(A=\operatorname{Cov}\mu\) and \(M=L^TL\). This use of \(A\) is local to the affine theorem, not the scalar norm used above.

## 4. Alias map and search translations

| Analysis | Probability | Combinatorics/algebra | Mathematical physics | Optimisation/statistics/ML |
|---|---|---|---|---|
| Weighted divergence-free field; Gaussian Helmholtz projection; curl obstruction | Gaussian copula; mixed fourth moment; marginal-preserving rearrangement | Spectral minors; exterior product; Lagrange identity | Lotka–Volterra centre; skew-symmetric predator–prey system; invariant product measure | Squared-distance GW; Gaussian GW; linear Gromov–Monge; covariance alignment |
| Third Hermite chaos | Wick/Isserlis formula; Gaussian invariant flow | Rank-two skew matrix \(ba^T-ab^T\) | Hamiltonian/Poisson dynamics | Nonlinear improvement of a Gaussian or affine transport coupling |
| Similarity rigidity; closed differential form | Finite fourth moments; interior of support | Polynomial vanishing on an open set | Incompressible weighted transport | Necessary optimality conditions for affine encoders |

The combinatorial and representation-language translations do not generate an independent finite combinatorial object. Searches concentrated on the substantive aliases above instead of inventing one.

## 5. Search log

All entries below were searched on 31 October 2026. Quotation marks and Boolean terms are recorded so the query families can be repeated. Search engines may broaden or ignore operators; a returned item was not treated as evidence without checking its mathematical object.

| Corpus/index | Exact query or query family | Result/consequence |
|---|---|---|
| Both live web indexes; arXiv and publisher hits | `Gaussian Gromov Wasserstein optimal coupling conjecture nonGaussian 2025 2026` | Located the originating Gaussian paper and Peyré's 2026 counterexample. |
| Freshness check, independent index, 365-day recency window | `Gaussian Gromov Wasserstein exactness classification affine rigidity 2026`; `Gaussian Gromov Wasserstein proportional covariance eigenvalues non-Gaussian optimal 2026` | Again located Peyré and the original Gaussian benchmark; no matching classification or rigidity result found. Recency metadata did not replace checking publication dates. |
| Live web | `"Gaussian" "Gromov-Wasserstein" "if and only if"`; `"Gromov-Wasserstein" "Gaussian" "proportional"`; `"Gromov-Wasserstein" "Gaussian" "non-optimal"`; `"Gromov-Wasserstein" "Gaussian" "counterexamples"` | Known sufficiency and distinct inner-product results found; no matching necessity classification located. |
| Live web/full text | `"Gromov-Wasserstein" "Gaussian" "first variation"`; `"Gromov-Wasserstein" "Gaussian" "perturbation"`; `"Gromov-Wasserstein" "Gaussian" "measure-preserving"`; `"Gromov-Wasserstein" "Gaussian" "collinear"` | Originating moment formulas and bound comparison; no exact flow-derived necessity statement located. |
| Live web/full text | `"Gromov-Wasserstein" "affine" "optimal" "similarity"`; `"Gromov-Wasserstein" "linear" "optimal" "isometry"`; `"Gromov-Monge" "affine"`; `"Gromov-Wasserstein" "curl"` | General Monge-map structure and representation-learning papers, not the proposed affine classification. |
| Second independent index | `"Gromov-Wasserstein" "affine" ("similarity" OR "rigidity" OR "curl")` | No matching theorem found; contemporary application and solver papers returned. |
| Live web/full text | `"Lotka-Volterra" "Gromov"`; `"Lotka-Volterra" "Gaussian measure"`; `"Gaussian" "measure-preserving" "Lotka-Volterra"`; `"Lotka-Volterra" "Gaussian" "skew-symmetric"` | Classical dynamical-system literature and unrelated Gaussian-noise models. No application of this flow to the present GW classification located. |
| Live web/full text | `"Gaussian-preserving" "flow" cubic`; `"Gaussian preserving" "cubic" "divergence"`; `"Gaussian" "divergence-free" "Gromov"` | Established Gaussian-preserving rearrangement methods found, including Morel et al. |
| Both indexes | `"Gromov-Wasserstein" "Gaussian" "solenoidal"`; `"Gromov-Wasserstein" "Gaussian" "Helmholtz"`; `"Gromov-Wasserstein" "Gaussian" ("divergence-free" OR "solenoidal" OR "Helmholtz")`; `"Gromov-Wasserstein" "Hermite"` | No matching exact sensitivity formula found. Some hits used Helmholtz/Hermite in institutional names rather than mathematics. |
| Both indexes | `"Gromov-Wasserstein" "Gaussian" ("proportional spectra" OR "proportional eigenvalues" OR "affine rigidity")` | Existing proportional sufficiency verified. |
| Live web/full text | `"Gromov-Wasserstein" "Gaussian" ("24" "32" OR "sqrt(6)" OR "76/7")`; `"Gromov-Wasserstein" "Gaussian" "dimension-free" "gap"`; `"Gromov-Wasserstein" "Gaussian" "collinearity" "gap"`; `"Gromov-Wasserstein" "Gaussian" "norm" "8" "lower bound"` | Known upper/lower benchmark comparison found; no proposed positive lower-gap estimate or pair certificate located. |
| Live web and OEIS | `"Gromov-Wasserstein" "Gaussian" "0, 12, 32, 60"`; `site:oeis.org "0, 3, 8, 15, 24, 35"` | Scaled benchmark sequence matches A005563. No additional GW result follows from this elementary sequence identity. |
| Original source, cited references, proceedings | Inspection of 2021/2022 Gaussian paper, 2024/2025 Monge-map paper, 2022 entropic Gaussian paper, 2025 Gromov–Monge Gap paper | Normalisations and scope checked. No matching classification found in inspected passages. |
| Author webpage, downloaded book and repository | Peyré counterexample post; OT4ML HTML chapter, current PDF and chapter audit | Full posted counterexample inspected. After an initial size-limit failure, the current PDF was downloaded; printed pp. 302–304 were read and the full extracted text was searched. See S8 for its hash. |

## 6. Inspected primary sources and verified collisions

### S1. Originating Gaussian paper: known formulas and sufficiency

Delon, J., Desolneux, A., and Salmona, A. (2022), *Gromov–Wasserstein distances between Gaussian distributions*, **Journal of Applied Probability** 59(4), 1178–1198. DOI: <https://doi.org/10.1017/jpr.2022.16>. The openly inspected preprint has author order Salmona–Delon–Desolneux: <https://arxiv.org/pdf/2104.07970>.

Load-bearing locators: Theorem 2.1 (moment reduction), Theorem 4.1 (Gaussian optimum), Proposition 5.2 (proportional sufficiency), Section 5.3 (degenerate supports), Eq. (5.3)–(5.4) (spectral bound difference), and Section 6, PDF pp. 15–16 (numerical conjecture).

Verified collision: the benchmark

\[
Q_G=4(\sum a-\sum b)^2+8\sum_i(a_i-b_i)^2
\]

and the upper-gap bound \(Q_G-Q_*\le8(AB-c)\) are existing results after the translations in Section 2. Their special-case exactness is not a new claim. The present proposed contribution is necessity and a strictly positive lower-gap bound outside that case.

### S2. Recent non-Gaussian counterexample: must be credited

Peyré, G. (23 October 2026), *Gaussians are not closed for Gromov–Wasserstein*: <https://www.gpeyre.com/blog/2026/10/23/gaussian-gromov-wasserstein-counterexample/>.

The full author post was inspected. Its radial quantile construction compares one- and \(n\)-dimensional standard Gaussians and beats the Gaussian value for sufficiently large \(n\). It attributes the construction to Proposition 13.34 of *Optimal Transport for Machine Learners* and leaves the smallest failing dimension and unrestricted optimiser unresolved. Therefore a generic counterexample is already a collision; the proposed universal classification and rank-two example go further than the posted result.

### S3. General GW linearisation and Monge structure

Dumont, T., Lacombe, T., and Vialard, F.-X., *On the existence of Monge maps for the Gromov–Wasserstein problem*, arXiv:2210.11945v3, 29 July 2024, subsequently *Foundations of Computational Mathematics*: <https://arxiv.org/pdf/2210.11945>.

The full indexed PDF was available. Section 1.2.2 discusses the Gaussian question; Section 3.2 gives the squared-distance linearised cost; Theorem 3.6 treats map/2-map structure under its support and density hypotheses. The linearisation is established background. The affine curl criterion and similarity classification were not found in the inspected material. The current proof uses only the necessary linearisation and proves its needed finite-cycle condition directly.

### S4. Gaussian-preserving flows are established methods

Morel, G., Drumetz, L., Benaichouche, S., Courty, N., and Rousseau, F. (2023), *Turning Normalizing Flows into Monge Maps with Geodesic Gaussian Preserving Flows*, **Transactions on Machine Learning Research**: <https://arxiv.org/pdf/2209.10873>.

The 31-page indexed paper was inspected. It develops Gaussian-preserving rearrangements to improve maps for ordinary quadratic transport, using divergence-free constructions and polar factorisation. This prevents claiming novelty for the general rearrangement method. Its objective and conclusions differ from the squared-distance GW classification and spectral certificates here. Bibliographic check: the PDF names TMLR publication in March 2023; arXiv version 4 is dated 14 April 2023. The journal review record is <https://openreview.net/forum?id=2UQv8L1Cv9>.

### S5. Similar titles, different inner-product objective

Le, K., Le, D. Q., Nguyen, H., Do, D., Pham, T., and Ho, N. (2022), *Entropic Gromov–Wasserstein between Gaussian Distributions*, **PMLR** 162, 12164–12203: <https://proceedings.mlr.press/v162/le22a.html>, paper <https://proceedings.mlr.press/v162/le22a/le22a.pdf>.

The primary proceedings description and relevant objective were inspected. Its Gaussian optimality theorem concerns inner-product GW, including an entropy penalty. The missing mixed-radius term makes it a different problem.

Dandapanthula, S., Podkopaev, A., Kasiviswanathan, S. P., Ramdas, A., and Goldfeld, Z. (2025), *Optimal Transportation and Alignment Between Gaussian Measures*: <https://arxiv.org/abs/2512.03579>. The abstract was inspected; the direct author PDF was not successfully fetched. Its advertised ordinary-OT and inner-product-GW scope does not give a squared-distance classification.

### S6. A sphere result with a Gaussian subproblem

Arya, S., Auddy, A., Clark, R. A., Lim, S., Mémoli, F., and Packer, D., *The Gromov–Wasserstein Distance Between Spheres*, published online 2024; **Foundations of Computational Mathematics** 26 (2026), 75–130: <https://doi.org/10.1007/s10208-024-09678-3>.

Full primary HTML inspected, particularly Eq. (7)–(8), Section 2.3.3 and Proposition 2.16. Its Gaussian subproblem maximises \(J(\pi)=\|E[XY^T]\|_F^2\). For spheres the radius term is constant; for Gaussian measures it is not. Thus the Gaussian projection result is not a collision with the present non-optimality claim for the full squared-distance objective.

### S7. Adjacent Gromov–Monge Gap literature

Uscidda, T., Eyring, L., Roth, K., Theis, F., Akata, Z., and Cuturi, M. (2025), *Disentangled Representation Learning with the Gromov-Monge Gap*, **ICLR**: <https://proceedings.iclr.cc/paper_files/paper/2025/file/c4a581f8219f62e3e5c83d4d80be5c9d-Paper-Conference.pdf>.

Its objective, Proposition 3.2, and Appendix B.4 were inspected. No affine similarity classification was located. The present proof does not invoke its support-inheritance proposition.

### S8. Current OT4ML book: full PDF acquired and relevant section inspected

Peyré, G., *Optimal Transport for Machine Learners*: <https://www.gpeyre.com/ot4ml/OT4ML/OT4ML.pdf>, retrieved 31 October 2026. An initial web-tool size failure was subsequently resolved by downloading the PDF through the available research environment. The inspected file has **27,846,433 bytes** and SHA-256

`b47340a8f0a9ba8d9b838d9c412fbf67e2d6804c6c5b788cc444b4d675b8528a` (an earlier copy of this report truncated it to 58 characters; corrected 1 November 2026 from the release verification record).

Reproducibility locations at inspection: `/tmp/gw_literature_current/OT4ML.pdf` and extracted text `/tmp/gw_literature_current/OT4ML.txt`.

Printed **pp. 302–304**, extracted-text lines 20310–20485, were read directly. Proposition **13.33** records the Gaussian benchmark, including singular covariance handling; Proposition **13.34** is precisely the sufficiently-large-dimensional radial counterexample described in S2. The subsequent paragraph credits known proportional-spectrum sufficiency and says the general Gaussian-optimality conjecture is disproved. No necessity classification, affine rigidity theorem or positive lower-gap estimate appears in that passage.

The whole extracted book was also searched for `Gaussian`, `GW`, `proportional`, `Gaussian-constrained`, `Gaussian optimal`, `Gaussian-optimal`, `Gaussian closure` and `affine rigidity`. No matching proposed criterion was located. This targeted section inspection plus whole-text search does not amount to reading every page. The earlier HTML and repository audit are supporting records; the current book is no longer an inaccessible primary-source gap.

### S9. Sequence identity

OEIS A005563: <https://oeis.org/A005563>. Its formula \(n(n+2)\), with index starting at zero, agrees with the benchmark sequence divided by four. The equality was verified algebraically. It does not identify any additional GW theorem.

### S10. Foundational context cited in the manuscript

Mémoli, F. (2011), *Gromov–Wasserstein Distances and the Metric Approach to Object Matching*, **Foundations of Computational Mathematics** 11, 417–487: <https://doi.org/10.1007/s10208-011-9093-5>.

Publisher title, author, date, pagination and abstract were checked. The publisher's full text is subscription restricted; it was not used as a source for a new technical inference. Its abstract directly supports the manuscript's brief account of GW as comparing metric-measure objects through their internal geometry.

## 7. Negative evidence and missingness audit

The recorded searches located no exact prior statement of the proposed Gaussian necessity theorem, general affine rigidity, strictly positive dimension-free lower gap, pair certificates, or exact latent solenoidal sensitivity. This is a bounded negative result, not proof of priority.

Coverage gaps remain:

- The current OT4ML PDF was acquired and its relevant Gaussian section inspected; the remainder was searched by text rather than read page by page. Equivalent terminology outside the searched patterns could still evade retrieval.
- The Mémoli foundational paper's publisher page was accessible only as metadata and abstract; the 2025 Dandapanthula et al. paper was checked at abstract level. Neither is load-bearing for the new necessity or lower-gap proof.
- No exhaustive MathSciNet, zbMATH, Web of Science, Scopus or Google Scholar citation export was performed.
- Searches were predominantly in English, with author names and mathematical formulas; non-English literature was not systematically searched.
- Formula indexing is imperfect. Equivalent statements under unrelated notation can be missed even after the aliases above.
- No author enquiry was sent. Such correspondence would require separate authorisation to send messages, and is unnecessary for constructing a reviewable manuscript now.
- The newest quantitative simultaneous-flow estimate was checked after its derivation; it has fewer independent search passes than the earlier classification and affine theorem.
- “Not located” concerns written statements in the inspected corpus; it says nothing about unpublished knowledge or priority dates.

## 8. Outcome by contribution unit

| Unit | Gate outcome | Confidence and allowed wording | Forbidden wording |
|---|---|---|---|
| Gaussian restricted optimum, moment expansion and proportional sufficiency | **COLLISION / background** | Exact formula and source translation verified; credit S1. | Calling these formulas or the sufficiency direction new. |
| Universal claim that some Gaussian pairs need non-Gaussian couplings | **COLLISION / background** | S2 already establishes large-dimensional examples. | “We give the first counterexample.” |
| Necessity classification for all positive-rank Gaussian pairs | **CLEAR TO PROVE within recorded corpus; priority BOUNDED UNCERTAINTY** | No collision located; distinguish exactness classification from solving all optimisers. | “The unrestricted Gaussian GW problem is fully solved.” |
| General affine similarity rigidity, including interior-support theorem | **CLEAR TO PROVE within recorded corpus; priority BOUNDED UNCERTAINTY** | No matching statement located; known linearisation credited. | Extending the statement to arbitrary finite clouds without hypotheses. |
| Pair certificate and explicit rank-two improvement | **CLEAR TO PROVE within recorded corpus** | No matching bound located; full proof and exact normalisation required. | Presenting its bound as the exact unrestricted optimum. |
| Dimension-free positive lower-gap estimate | **CLEAR TO PROVE within recorded corpus; recent derivation** | No matching lower estimate located. The complementary upper estimate is old. | Claiming the entire two-sided comparison was absent from prior literature. |
| Exact solenoidal projection and specified maximal slope | **CLEAR TO PROVE within recorded corpus** | Exact formula not located; specify latent \(L^2(\gamma)\) norm. | Claiming a new general Helmholtz theorem or an intrinsic GW gradient formula. |
| Lotka–Volterra dynamics and Gaussian-preserving rearrangement method | **COLLISION / established mechanism** | Algebraic identification verified; S4 supplies relevant methodological precedent. | Claiming a new vector field, new invariant-measure principle, or first use of Gaussian-preserving flows in transport. |
| Complete constructive proof assembled here | **CLEAR within recorded boundary** | Can call it a complete independently checked derivation; proof review is separate from novelty. | “First complete proof” without resolving priority gaps. |
| Exact symbolic checks and numerical demonstration | **Verification artefact** | Reproducible support for identities and examples; not proof of historical novelty. | Treating successful numerical tests as proof of optimality or priority. |
| Similarity characterisation of all proportional-case optimisers | **Standard equality-case consequence; no separate priority claim** | Useful corollary with a complete support argument. | Inflating the Cauchy–Schwarz equality case into an independently new principle. |

## 9. Decision

**Proceed with the proof and a candid manuscript.** The contribution-specific search is clear enough to develop and review the proposed classification, affine rigidity and quantitative lower bounds. **Historical novelty remains BOUNDED UNCERTAINTY.** The manuscript should explicitly credit S1–S4, identify the distinct inner-product objective, and state that the unrestricted optimiser away from proportional spectra remains undetermined.

The strongest proposed advance is the combined structural and quantitative result: a necessary-and-sufficient Gaussian criterion, a general affine theorem under a natural support condition, and a dimension-free positive lower bound that complements the pre-existing upper bound. A further isolated counterexample would be substantially weaker because S2 has already established the failure phenomenon. The Gaussian-preserving mechanism, similarity sufficiency and elementary spectral identities are supporting tools. Claims of practical algorithmic superiority or improved downstream prediction are not supported by this theoretical and verification package.

The next useful priority check is human expert inspection of the exact statements, especially by researchers familiar with S1–S3. The current OT4ML Gaussian section has now been inspected and creates no identified collision. No message has been sent on the user's behalf.

## 10. Rerun triggers

Rerun this gate if the proof produces an exact unrestricted optimiser, sharp gap constant, a closed-form all-time flow value, a substantially weaker support hypothesis, a non-Gaussian universality theorem, or an identification with another named Hamiltonian or transport construction. Also rerun if a later version of an inspected paper or book appears.

The gate was already rerun after the following simplifications: the two-coordinate flow, its Lotka–Volterra identity, the affine curl polynomial, the Hermite/solenoidal decomposition, the tensor curvature improvement, and the simultaneous-flow dimension-free gap.


## 11. Addendum, 1 November 2026: second search pass during internal review

A separate domain review re-ran the search for each proposed contribution. It read in full or searched the following primary texts: arXiv:2104.07970 (only version v1, all sections); Salmona's doctoral thesis (HAL tel-04793664, 2023), Chapter 4; Dumont–Lacombe–Vialard (arXiv:2210.11945v3; *Found. Comput. Math.* 25(2), 2025, doi:10.1007/s10208-024-09643-0); Salmona–Delon–Desolneux, arXiv:2310.11256; Dandapanthula et al., arXiv:2512.03579v2; Le et al., arXiv:2108.10961; Zhang–Goldfeld–Mroueh–Sriperumbudur, *Ann. Statist.* 52(4), 2024, doi:10.1214/24-AOS2406; Rioux–Goldfeld–Kato, arXiv:2306.00182; Beinert–Heiss–Steidl, arXiv:2205.09006; Mémoli–Needham, arXiv:2212.14123; Arya et al., arXiv:2306.10586; Vayer's thesis, arXiv:2011.04447. It also listed all Semantic Scholar citers of Delon et al. and searched the 2025–2026 mathematical ones, including Tran–Soh (arXiv:2608.24184).

| Unit | Outcome of second pass |
|---|---|
| Necessity direction of the Gaussian classification | Not found. Delon et al. conjecture that the Gaussian-restricted plan solves the unrestricted problem (end of §6 and §7 of v1). Salmona's thesis (§4.2.5) describes the construction of counterexamples as open. Peyré covers one versus n dimensions, n large. Tran–Soh (2026, §6.1) state that the general Gaussian value is not known. |
| Affine rigidity (only-if direction) | Not found as a statement. Ingredients known: linearisation (Dumont et al., eq. (1.1); compare the variational form of Zhang et al., Corollary 1) and a sufficient radial condition for optimal maps (Vayer, Proposition 4.2.4), which for linear maps is the similarity condition. Now cited. |
| Positive lower gap | Not found. The upper bound is Delon et al. eqs. (5.3)–(5.4) with Propositions 4.2 and 5.3. |
| New in revision: factor-two corollary, its sharpness, comonotone refinement, plug-in inconsistency | The upper half of the factor-two statement is an immediate consequence of the known upper bound and is credited as such. The comonotone refinement uses the classical one-dimensional rearrangement inequality. The stability lemma uses the standard gluing argument for GW distances (Mémoli, 2011). These units were not separately searched beyond the queries above; they are presented as consequences, not as independent priority claims. |

The bounded-search caveats of Section 7 still apply.
