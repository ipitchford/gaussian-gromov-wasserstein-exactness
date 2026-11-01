---
title: "Affine rigidity and Gaussian exactness in squared-distance Gromov–Wasserstein transport"
subtitle: "When the closed-form Gaussian answer is exact, and how much it misses when it is not"
author: "Anonymous"
date: "1 November 2026 · Version 0.1.0-candidate · [doi:10.5281/zenodo.23850758](https://doi.org/10.5281/zenodo.23850758)"
lang: en-GB
---

# Abstract

Gromov–Wasserstein (GW) transport compares two datasets or shapes that may live in different spaces by matching their internal pairwise distances. It is expensive to compute, so a closed-form answer for Gaussian data, found by Delon, Desolneux and Salmona (2022) after restricting attention to Gaussian couplings, is an attractive shortcut. Those authors conjectured that the restriction loses nothing. For the standard squared-distance objective we show exactly when this is true: two Gaussians of positive rank admit an optimal coupling that is itself Gaussian **if and only if** their covariances have the same rank and proportional ordered eigenvalues. In every other case, including equal dimensions, an explicit nonlinear coupling that keeps both marginals exactly Gaussian does strictly better. With \(a,b\) the zero-padded eigenvalue lists and \(S=\|a\|^2\|b\|^2-(a\cdot b)^2\), the Gaussian excess \(D=Q_G-Q_*\) satisfies \(S/(\|a\|\|b\|)\le D\le8(\|a\|\|b\|-a\cdot b)\); the upper bound is due to Delon et al., and the lower bound is new. Consequently the Gaussian formula never overstates the optimum by more than a factor of two, and this factor is approached by nearly isotropic spectra in high dimension. For a standard Gaussian plane matched to a standard Gaussian line, the Gaussian answer is 12, an explicit coupling is certified below \(76/7\approx10.86\), and the true optimum lies between 8.68 and 10.85 (and is at least about 8.82 by a numerically evaluated refinement). The necessity direction comes from a general theorem: for a source whose support has interior points, an affine map is GW-optimal onto its image only if it is a similarity. We also show that the plug-in Gaussian formula is an inconsistent estimator of the population GW value. The unrestricted optimum for nonproportional spectra is not computed.

**Keywords:** Gromov–Wasserstein; optimal transport; Gaussian measures; Gaussian coupling; affine rigidity; linear Gromov–Monge; measure-preserving flow; approximation error.

# 1. Introduction

## 1.1 What GW transport does, and why a Gaussian shortcut matters

Ordinary optimal transport moves mass from one distribution to another along a shared ground cost, which requires both distributions to live in the same space. Gromov–Wasserstein (GW) transport (Mémoli, 2011) removes that requirement. It compares two distributions that may live in different spaces, for example a point cloud in the plane and one on a line, a three-dimensional shape and a two-dimensional scan, or two families of measurements of the same cells. Instead of moving mass along a ground cost, it seeks a coupling under which pairwise squared distances in one space match pairwise squared distances in the other as closely as possible. The matching is unaffected by rotations or translations of either dataset, which is why GW is used for shape matching, graph matching and the alignment of datasets or representations with no common coordinates.

The price is computational. The GW problem is a nonconvex quadratic programme over couplings, and practical solvers find local solutions. For Gaussian data, Delon, Desolneux and Salmona (2022) found an attractive closed form by restricting the search to Gaussian couplings. Anyone who summarises each dataset by its covariance and plugs the eigenvalues into that formula is implicitly assuming that the restriction is harmless. Delon et al. conjectured, on numerical evidence, that it is always harmless; Salmona (2023, §4.2.5) later described the construction of counterexamples as an open problem. Peyré (2026a, 2026b) has since shown that it fails for a standard Gaussian line against a standard Gaussian space of sufficiently high dimension. This paper determines exactly when the Gaussian answer is the true answer, and how large the error is when it is not.

**A concrete example.** Match a standard Gaussian in the plane, \(N(0,I_2)\), to a standard Gaussian on the line, \(N(0,1)\). The best Gaussian coupling simply projects onto one coordinate and has objective \(12\). Section 7 constructs a curved rearrangement of the plane, which keeps the target exactly standard Gaussian and has objective provably below \(76/7\approx10.86\), about 10.47 numerically. A better numerical coupling reaches about 10.03, while lower bounds show that the true optimum is at least 8.68, and about 8.82 with a numerically evaluated refinement. The Gaussian shortcut therefore exceeds the true optimum by at least 1.159, or 9.6% of the Gaussian value, and by at most about 3.18, or 26.5%.

## 1.2 Main results in plain words

> **Main results.**
>
> 1. *Exactness.* The closed-form Gaussian GW value is the true optimum if and only if one covariance is a rescaled copy of the other in its eigenvalues: same positive rank, eigenvalues in the same ratios (Theorem 5). In particular the conjecture of Delon et al. fails for every other pair, including pairs of the same dimension.
> 2. *Size of the error.* Otherwise the Gaussian value overstates the optimum by an amount between \(S/(\|a\|\|b\|)\) and \(8(\|a\|\|b\|-a\cdot b)\) (Theorem 8), and never by more than a factor of two (Corollary 9). Both bounds can be computed from the two eigenvalue lists; Section 6.3 gives a recipe.
> 3. *Why.* For a source with interior support, an affine map, such as a principal-component projection, is GW-optimal onto its image only if it multiplies all distances by one constant (Theorem 3). The mixed radial term in the objective is what rules out every other affine map.
> 4. *Samples.* Plugging sample covariances into the Gaussian formula gives an inconsistent estimate of the GW value between nonproportional Gaussian populations (Proposition 12).

## 1.3 What is known and what is new

**Known and reproduced to fix conventions.** The moment reduction (Lemma 1), the Gaussian-restricted value and its common-latent optimiser (Lemma 2), sufficiency of proportional spectra, and the upper inequality in Theorem 8 are due to Delon et al. (2022). We use the numbering of the openly available preprint, arXiv:2104.07970, version 1, which is its only version. The linearisation of the squared-distance problem used in Theorem 3 is established (Dumont et al., 2025, §1.1 and §3.2; compare the variational form of Zhang et al., 2024, Corollary 1). Vayer (2020, Proposition 4.2.4) gives a sufficient radial condition for optimal maps, which for linear maps reduces to the similarity condition. Gaussian-preserving rearrangements are a recognised method in ordinary optimal transport (Morel et al., 2023). The flow (10) is a classical Lotka–Volterra centre. Peyré's construction (2026a, Proposition 13.34; 2026b) already shows that Gaussian couplings are not always optimal, for one versus \(n\) dimensions with \(n\) sufficiently large.

**Proposed contributions.** The necessity direction of Theorem 5 and hence the exact classification; the affine rigidity theorem (Theorem 3) in its "only if" direction; the positive lower bound in Theorem 8 with the sharper certificates of Sections 6.1–6.2; the maximal Gaussian-preserving slope (Proposition 7); the factor-two corollary and its sharpness; the comonotone refinement of the known upper bound (Proposition 10); and the plug-in consequence (Proposition 12). A formula-level literature search, conducted on 31 October and 1 November 2026 and recorded in the accompanying novelty report, did not locate these statements. That record establishes a search boundary, not historical priority.

For Gaussians the necessity direction can also be read off from a single explicit computation (Lemma 6), and we say so where it is proved. We give the general affine theorem because it explains the mechanism and applies beyond Gaussian measures.

## 1.4 Roadmap and scope

Section 2 fixes notation and recalls the Gaussian benchmark. Section 3 proves affine rigidity. Section 4 deduces the Gaussian classification. Sections 5 and 6 construct the improving flow and prove the two-sided bound, the factor-two corollary and a practical recipe. Section 7 works the plane-to-line and an equal-dimension example. Section 8 treats samples, and Section 9 discusses interpretation and open problems.

The conclusions apply to objective (1). Inner-product GW has a different objective and admits Gaussian optimality much more broadly (Delon et al., 2022, Proposition 4.1; Le et al., 2022; Dandapanthula et al., 2026). We infer nothing about entropy-regularised GW or GW with ordinary rather than squared distances. Theorem 3 requires interior support and does not apply directly to finite point clouds; Section 3 shows why.

# 2. Moment reduction and the Gaussian benchmark

For probability measures \(\mu\) and \(\nu\) on Euclidean spaces, let \(\Pi(\mu,\nu)\) denote their couplings. We study

\[
Q_* = \min_{\pi\in\Pi(\mu,\nu)}Q(\pi),\qquad
Q(\pi)=\mathbb E\big[(\|X-X'\|^2-\|Y-Y'\|^2)^2\big],                 \tag{1}
\]

where \((X,Y),(X',Y')\) are independent with law \(\pi\). The inner dissimilarity is a **squared** Euclidean distance, the outer loss is also squared, and there is no prefactor. For empirical measures with uniform weights this is exactly the default objective of the POT library (Flamary et al., 2021), with squared Euclidean cost matrices and square loss. We report the objective itself, rather than a square root or fourth root. Finite fourth moments make every coupling admissible with finite objective. Tightness of the coupling set and lower semicontinuity of the nonnegative continuous integrand give existence of a minimiser. For Gaussian marginals, \(Q_G\) denotes the minimum of (1) over jointly Gaussian couplings.

Translations of either marginal do not change (1), so we centre both laws. Set

\[
A=\mathbb E[XX^\top],\quad B=\mathbb E[YY^\top],\quad
K_\pi=\mathbb E[XY^\top],\qquad
G(\pi)=\mathbb E[\|X\|^2\|Y\|^2]+2\|K_\pi\|_F^2.                \tag{2}
\]

**Lemma 1 (moment reduction).** For fixed centred marginals with finite fourth moments, minimising \(Q\) is equivalent to maximising \(G\), and

\[
Q(\pi)=C_{\mu,\nu}-4G(\pi),\qquad
C_{\mu,\nu}=\mathbb E\|X-X'\|^4+\mathbb E\|Y-Y'\|^4
             -4\operatorname{tr}A\operatorname{tr}B.              \tag{3}
\]

**Proof.** Expanding the product of squared distances and using centring and independence of the two sampled pairs gives

\[
\mathbb E[\|X-X'\|^2\|Y-Y'\|^2]
=2\mathbb E[\|X\|^2\|Y\|^2]
 +2\operatorname{tr}A\operatorname{tr}B+4\|K_\pi\|_F^2.
\]

Substitution in the expansion of (1) gives (3). This is the moment reduction of Delon et al. (2022, Theorem 2.1), written in our normalisation. \(\square\)

For Gaussian marginals, list the positive eigenvalues of \(A,B\) in nonincreasing order and pad with zeros to a common length \(d\), obtaining \(a,b\in[0,\infty)^d\). One may take \(d\) to be the larger covariance rank. Orthogonal coordinates identify the affine supports with the corresponding coordinate subspaces; deterministic zero coordinates can be removed or restored without changing costs. Throughout, \(\alpha=\sum_ia_i\), \(\beta=\sum_ib_i\), and \(\|a\|,\|b\|\) are Euclidean norms of the eigenvalue lists, not matrix norms.

**Lemma 2 (Gaussian benchmark and a covariance bound).** Among all couplings,
\(\|K_\pi\|_F^2\leq a\cdot b\). Among jointly Gaussian couplings the minimum is

\[
Q_G=4\left[(\alpha-\beta)^2+2\|a-b\|^2\right].                 \tag{4}
\]

It is attained by the common-latent coupling
\(X=D_a^{1/2}Z,\ Y=D_b^{1/2}Z\), where \(Z\sim N(0,I_d)\) and \(D_a=\operatorname{diag}(a)\).

**Proof.** On the positive supports, positive semidefiniteness of the joint covariance matrix gives
\(K=A^{1/2}RB^{1/2}\) with \(\|R\|_{\mathrm{op}}\leq1\). In eigenbases,
\(\|K\|_F^2=\sum_{ij}a_i b_j R_{ij}^2\). The matrix \(P=(R_{ij}^2)\), padded to a square if necessary, has row sums \(r_i\le1\) and column sums \(c_j\le1\). Let \(\delta=\sum_i(1-r_i)=\sum_j(1-c_j)\). If \(\delta>0\), the matrix
\[
P+\delta^{-1}(1-r)(1-c)^\top
\]
is doubly stochastic and entrywise at least \(P\); if \(\delta=0\), \(P\) is already doubly stochastic. Nonnegative coefficients allow this completion without decreasing the sum. Birkhoff's theorem and the rearrangement inequality, applied to a convex combination of permutation matrices, then bound the sum by \(a\cdot b\).

For a jointly Gaussian pair, the fourth-moment identity yields
\(\mathbb E\|X\|^2\|Y\|^2=\alpha\beta+2\|K\|_F^2\). Thus
\(G_G=\alpha\beta+4a\cdot b\), attained by the common-latent coupling. Finally, \(\mathbb E\|X-X'\|^4=4[(\operatorname{tr}A)^2+2\operatorname{tr}(A^2)]\); insertion into (3) gives (4). \(\square\)

Formula (4) and the corresponding affine construction are Theorem 4.1 of Delon et al. (2022). No novelty is claimed for them.

# 3. Affine rigidity beyond Gaussian measures

We first prove the general structural result. A similarity on a Euclidean space is a map that multiplies all distances by one fixed positive number; a constant map is treated separately.

**Theorem 3 (global affine rigidity).** Let \(\mu\) have finite fourth moment, and suppose its support has nonempty relative interior in its affine hull \(E\). Let \(T:E\to\mathbb R^m\) be affine and \(\nu=T_\#\mu\). Then the coupling \((\mathrm{id},T)_\#\mu\) minimises (1) if and only if \(T\) restricts to a similarity or a constant on \(E\). Equivalently, in orthonormal coordinates on the direction space of \(E\), the linear part \(L\) satisfies

\[
L^\top L=\lambda I\quad\text{for some }\lambda\geq0.              \tag{5}
\]

**Proof.** The zero-dimensional case is immediate. Identify \(E\) with \(\mathbb R^r\) and centre the variables, so \(T(x)=Lx\). Write \(A=\mathbb E[XX^\top]\) and \(M=L^\top L\).

*Sufficiency.* Under the affine coupling, \(\|Y-Y'\|^2=\lambda\|X-X'\|^2\). Under any competing coupling, the two squared distances have fixed marginal distributions with second moments \(m_2\) and \(\lambda^2m_2\), where \(m_2=\mathbb E\|X-X'\|^4\). Cauchy–Schwarz gives
\(\mathbb E[\|X-X'\|^2\|Y-Y'\|^2]\leq\lambda m_2\), attained by the affine coupling. Hence its objective \((1-\lambda)^2m_2\) is minimal. This argument does not require interior support; it is the elementary sufficiency direction, implicit in Delon et al. (2022, Proposition 5.2) and in Vayer's radial condition (Vayer, 2020, Proposition 4.2.4).

*Necessity: linearisation.* Let \(\pi_0=(\mathrm{id},L)_\#\mu\) be optimal. For any coupling \(\pi\), with \(K_0=AL^\top\) and \(\Delta K=K_\pi-K_0\), the exact expansion
\[
G((1-\epsilon)\pi_0+\epsilon\pi)-G(\pi_0)
=\epsilon\Big[\textstyle\int\|x\|^2\|y\|^2\,d(\pi-\pi_0)+4\langle K_0,\Delta K\rangle_F\Big]+2\epsilon^2\|\Delta K\|_F^2
\]
and optimality of \(\pi_0\) make the bracket nonpositive. Therefore \(\pi_0\) minimises the ordinary transport cost

\[
c_0(x,y)=-\|x\|^2\|y\|^2-4x^\top A L^\top y.                    \tag{6}
\]

All terms are integrable under every coupling: \(|c_0(x,y)|\leq(\|x\|^4+\|y\|^4)/2+2\|A\|\,\|L\|(\|x\|^2+\|y\|^2)\), and the right-hand side is integrable by the fourth-moment assumption. This is the linearised cost of Dumont et al. (2025, eq. (1.1)), and it appears in the variational form of Zhang et al. (2024, Corollary 1), up to marginal-only terms and a positive factor. We use only this necessary direction: optimality for the linearised cost alone does not control the positive quadratic remainder in the displayed expansion.

*Necessity: finite cycles.* For any distinct \(x_1,\ldots,x_N\) in \(\operatorname{supp}\mu\), with \(x_{N+1}=x_1\), optimality for the continuous cost (6) implies

\[
\sum_{k=1}^N c_0(x_k,Lx_k)
\leq\sum_{k=1}^N c_0(x_k,Lx_{k+1}).                              \tag{7}
\]

For completeness, suppose a strict violation. By continuity it persists on small bounded disjoint neighbourhoods \(U_k\) of the \(x_k\), each of positive \(\mu\)-mass. Let \(m>0\) be the smallest of these masses. From the graph plan restricted to each \(U_k\), remove the fraction \(m/\mu(U_k)\) of its mass, so that each removed piece has mass \(m\). Reconnect the normalised source part from \(U_k\) to the normalised target part from \(U_{k+1}\) by a product coupling of total mass \(m\). Both global marginals are unchanged, all modified costs are bounded on the bounded neighbourhoods, and continuity makes the total cost strictly smaller. This contradicts optimality and proves (7), including for measures with atoms, since fractions of atoms may be removed.

Fix coordinates \(i\ne j\) and choose a rectangle in the \((x_i,x_j)\)-plane inside a ball contained in \(\operatorname{supp}\mu\). Apply (7) to the vertices of increasingly fine partitions of its boundary, traversed in order. Taylor expansion in the second argument has a remainder bounded by a constant times \(\sum_k\|x_{k+1}-x_k\|^2\), which is at most the mesh times the perimeter and tends to zero. Since

\[
L^\top\nabla_y c_0(x,Lx)=-g(x),\qquad
g(x)=2\|x\|^2Mx+4MAx,                                          \tag{8}
\]

the left-endpoint sums converge to a line integral of the continuous field \(g\), which therefore has one sign. Reversing the curve gives the other sign, so \(\oint g\cdot dx=0\). Small rectangles therefore give \(\partial_jg_i-\partial_ig_j=0\) on the ball. Direct differentiation gives

\[
\partial_jg_i-\partial_ig_j
=4[x_j(Mx)_i-x_i(Mx)_j]+4(MA-AM)_{ij}.                           \tag{9}
\]

This polynomial vanishes on an open set, so all its coefficients vanish. Its degree-two part alone implies that \(Mx\) is parallel to \(x\) for every \(x\); the covariance \(A\) plays no role. Taking \(x=e_i\) makes \(M\) diagonal; taking \(x=e_i+e_j\) equates every diagonal entry. Thus \(M=\lambda I\), with \(\lambda\geq0\). For \(r=1\), every such \(M\) is already scalar. \(\square\)

The theorem concerns affine graph couplings. It does not say what the optimiser is when the affine coupling fails, and the cycle argument it uses is standard optimal-transport technique; what is new is the conclusion.

**Corollary 4 (local descent under a density assumption).** Suppose additionally that on some nonempty open set \(U\), the restricted measure is exactly \(\mu|_U=\rho(x)dx\), where \(\rho>0\) is continuously differentiable. If \(M\) is not scalar, compactly supported, continuously differentiable \(\mu\)-preserving flows produce strict improvements arbitrarily close to the affine coupling.

**Proof.** For \(\psi\in C_c^\infty(U)\), choose a two-coordinate vector field
\(V_i=\rho^{-1}\partial_j\psi,\ V_j=-\rho^{-1}\partial_i\psi\), with other components zero. It is continuously differentiable with compact support and \(\operatorname{div}(\rho V)=0\). Its complete flow preserves \(\mu\). Composing the affine map with this flow gives \(G'(0)=\int g\cdot V\,d\mu\), equal, up to the chosen orientation, to the integral of (9) against \(\psi\). A nonscalar \(M\) makes at least one curl polynomial nonzero on \(U\). Choose \(\psi\) supported where it has a fixed sign. Reversing the flow if necessary gives \(G'(0)>0\). \(\square\)

**Why the support condition matters.** Let \(X\) be uniform on the four points \((\pm1,\pm1)\), and let \(Y=X_1\). This dimension-reducing affine coupling is optimal. Indeed, for any coupling to a uniform sign \(Y\), write \(f(X)=\mathbb E[Y\mid X]\). Since \(X_1,X_2\) are orthonormal in \(L^2(\mu)\),
\((\mathbb E X_1Y)^2+(\mathbb E X_2Y)^2\leq\mathbb E f(X)^2\leq1\). Both marginal squared norms are constant, so (2) is at most \(2+2=4\), attained by \(Y=X_1\). Its objective is \(8\). The support has empty interior in its two-dimensional affine hull. Thus the theorem cannot be applied to arbitrary finite data clouds.

# 4. Exact Gaussian classification

**Theorem 5 (Gaussian exactness).** Let both Gaussian marginals have positive covariance rank. The following are equivalent:

1. Some unrestricted minimiser of (1) is jointly Gaussian.
2. \(Q_*=Q_G\).
3. Their positive covariance ranks agree and their ordered positive covariance eigenvalues are proportional.

If (3) holds, every unrestricted optimiser is induced by an affine similarity between the two Gaussian supports. If either covariance has rank zero, there is a unique coupling and it is Gaussian, without a rank-matching requirement.

**Proof.** The equivalence of (1) and (2) follows from Lemma 2 and its attained benchmark. Orient the common-latent Gaussian optimum from the larger positive support, of rank \(r\), to the smaller, of rank \(s\leq r\). It is the graph of the linear map with matrix \([\operatorname{diag}(\sqrt{b_1/a_1},\ldots,\sqrt{b_s/a_s})\ \ 0]\), and its source Gaussian has full support on that larger space. If \(Q_*=Q_G\), Theorem 3 forces that map to be a similarity or constant. Positive target rank excludes a constant. If \(r>s\), \(L^\top L\) has both positive and zero eigenvalues and cannot be scalar; if \(r=s\), being scalar is exactly proportionality of the ordered spectra. Conversely, proportional spectra permit a similarity coupling, which is optimal by the sufficiency proof of Theorem 3. This sufficiency is Proposition 5.2 of Delon et al. (2022), with degenerate supports treated in their Section 5.3.

To identify all optimisers in the proportional case, write \(b=c a\), \(c>0\). Equality in the distance Cauchy–Schwarz bound forces
\(\|y-y'\|^2=c\|x-x'\|^2\) for \(\pi\otimes\pi\)-almost every pair. The equality is continuous and therefore holds on \(\operatorname{supp}\pi\times\operatorname{supp}\pi\). Two support pairs with the same \(x\) must have the same \(y\). The projection of \(\operatorname{supp}\pi\) onto the source is \(\sigma\)-compact, hence Borel, and has full \(\mu\)-measure, so it is dense in the source support. The graph defines a scaled isometry on this dense set. It extends continuously to the entire support; closedness of the joint support puts each resulting limit pair back in that support. A distance-preserving map from a Euclidean space into another Euclidean space is affine: subtract its value at zero and use the polarisation identity on an orthonormal basis. Equal support dimensions and the target marginal then give an affine similarity onto the target support. \(\square\)

**The conjecture of Delon et al.** On numerical evidence, Delon et al. (2022, end of Section 6 and Section 7) conjectured that the optimal Gaussian-restricted plan also solves the unrestricted problem, for all Gaussian pairs. Theorem 5 shows that this holds exactly on the proportional set and fails everywhere else. Peyré's construction (2026a, Proposition 13.34; 2026b) covers a standard Gaussian line against a standard Gaussian space of sufficiently large, unspecified dimension. Theorem 5 adds every dimension \(n\geq2\) for that pair, so that \(Q_*<4(n^2-1)\) already for \(n=2\), and every nonproportional pair of equal dimension, a regime the radial construction does not reach.

For Gaussian marginals the necessity direction does not need Theorem 3. Lemma 6 below gives, by one explicit expectation, a Gaussian-preserving direction of strict descent whenever the spectra are nonproportional. Theorem 3 identifies the general mechanism; Lemma 6 makes it quantitative.

# 5. Constructive Gaussian-preserving descent

The following construction proves strict Gaussian failure directly and supplies quantitative bounds. Let \(\gamma_d=N(0,I_d)\). For \(i<j\), define

\[
V^{ij}_i(z)=z_i(z_j^2-1),\qquad
V^{ij}_j(z)=z_j(1-z_i^2),                                      \tag{10}
\]

with all other components zero.

**Lemma 6 (admissible pair flow).** The field (10) has a complete smooth flow \(T_t\), \(t\in\mathbb R\), preserving \(\gamma_d\). For the feasible coupling
\(X=D_a^{1/2}Z,\ Y_t=D_b^{1/2}T_t(Z)\),

\[
G'(0)=4\Delta_{ij},\qquad Q'(0)=-16\Delta_{ij},\qquad
\Delta_{ij}=b_i a_j-b_j a_i.                                  \tag{11}
\]

**Proof.** Direct calculation gives
\(\operatorname{div}V^{ij}=z\cdot V^{ij}=z_j^2-z_i^2\), hence
\(\operatorname{div}(\varphi_d V^{ij})=0\), where \(\varphi_d\) is the standard Gaussian density. Moreover,

\[
\frac{d}{dt}\|T_t(z)\|^2=2(T_{t,j}^2-T_{t,i}^2),\qquad
\|T_t(z)\|\leq e^{|t|}\|z\|.                                  \tag{12}
\]

Local existence and this bound give a complete flow in both time directions. The density–Jacobian identity then proves exact Gaussian preservation.

Each entry of \(\mathbb E[ZV^{ij}(Z)^\top]\) is zero, so the covariance-squared term has zero derivative. Put \(u=Z_i,v=Z_j\). The derivative of the mixed norm moment is

\[
2\mathbb E[(a_i u^2+a_j v^2)
 \{b_i u^2(v^2-1)+b_j v^2(1-u^2)\}]
=4(b_i a_j-b_j a_i).
\]

Terms involving other coordinates vanish by independence. Equation (12) and the polynomial field bound give Gaussian-integrable polynomial dominators for differentiation on compact time intervals. Equation (3) completes the proof. \(\square\)

Nonproportional nonzero vectors have a nonzero two-by-two minor \(\Delta_{ij}\). Choosing the sign of time in (11) gives strict improvement over the **whole** Gaussian class, because the starting coupling attains (4). No numerical optimisation is involved in this conclusion. Under \(p=u^2,q=v^2\), the two-coordinate dynamics are \(p'=2p(q-1),q'=2q(1-p)\), a classical Lotka–Volterra centre. We claim no novelty for that dynamical system.

## 5.1 Exact first-order sensitivity

Let

\[
C=ba^\top-ab^\top,\quad
S=\sum_{i<j}C_{ij}^2=\|a\|^2\|b\|^2-(a\cdot b)^2,\quad
W_i(z)=z_i\sum_j C_{ij}(z_j^2-1).                              \tag{13}
\]

**Proposition 7 (maximal Gaussian-preserving slope).** In latent standard-Gaussian coordinates, the maximal first variation of \(G\) over Gaussian-divergence-free vector fields of \(L^2(\gamma_d)\)-norm at most one is \(2\sqrt S\). The corresponding maximal descent slope of \(Q\) is \(8\sqrt S\). For \(S>0\), the field \(W/(2\sqrt S)\) attains these slopes and generates a complete Gaussian-preserving flow.

**Proof.** The fields \(V^{ij}\) are mutually orthogonal in \(L^2(\gamma_d)\), each with squared norm \(4\). Therefore \(W=\sum_{i<j}C_{ij}V^{ij}\) has \(\|W\|_2^2=4S\). Put \(H_a(z)=\sum a_i z_i^2\) and \(H_b(z)=\sum b_i z_i^2\). The derivative field for the mixed term is \(f=2H_aD_bz\), and direct differentiation gives

\[
f=\nabla\Phi+W,\qquad
\Phi=\tfrac12H_aH_b+\tfrac12(\alpha H_b-\beta H_a).               \tag{14}
\]

For a Gaussian-divergence-free field \(V\), defined weakly by
\(\int\nabla\psi\cdot V\,d\gamma_d=0\) for all smooth compactly supported \(\psi\), the gradient term is orthogonal to \(V\). Smooth cutoffs of the polynomial \(\Phi\) justify this identity in \(L^2(\gamma_d)\). The same weak condition, tested on cutoffs of \(\psi=z_iz_j\) and \(\psi=z_i^2/2\), gives \(\mathbb E[Z_iV_j(Z)+Z_jV_i(Z)]=0\) for all \(i,j\), so the first variation of the cross-covariance \(\mathbb E[ZT_t(Z)^\top]\) is antisymmetric, and the first variation of the covariance-squared term vanishes at the common-latent coupling, whose cross-covariance is diagonal. Hence \(G'[V]=\langle W,V\rangle\), and Cauchy–Schwarz proves the stated supremum. For a weak field that generates no flow, \(G'[V]\) denotes this formal first variation, \(\langle f,V\rangle\) plus the covariance term; it is a linear tangent functional; the maximising polynomial field supplies an actual admissible flow.

Finally, skew symmetry of \(C\) cancels all quartic terms in \(z\cdot W\), leaving
\(z\cdot W=\sum_i(\beta a_i-\alpha b_i)z_i^2\). A Gronwall bound proves completeness, and weighted divergence zero proves preservation as in Lemma 6. \(\square\)

# 6. Quantitative Gaussian approximation error

Write \(D=Q_G-Q_*\) for the Gaussian excess, and assume \(\|a\|,\|b\|>0\).

**Theorem 8 (dimension-free two-sided error).** For all Gaussian pairs of positive covariance rank,

\[
\boxed{\quad
\frac{S}{\|a\|\|b\|}\ \leq\ D\ \leq\ 8(\|a\|\|b\|-a\cdot b).
\quad}                                                        \tag{15}
\]

Both sides vanish exactly for proportional spectra. The upper inequality is the comparison of Delon et al. (2022, eqs. (5.3)–(5.4), with Propositions 4.2 and 5.3). The lower inequality is constructive: the flow of \(W\) in (13), run for time \(1/(8\|a\|\|b\|)\), achieves that improvement or more.

**Proof of the upper inequality.** The Gaussian marginal norm variances are \(2\|a\|^2\) and \(2\|b\|^2\). Cauchy–Schwarz and Lemma 2 give

\[
G(\pi)\leq\alpha\beta+2\|a\|\|b\|+2a\cdot b.
\]

Subtract the attained Gaussian score and use (3).

**Proof of the lower inequality.** The case \(S=0\) follows from Theorem 5. Suppose \(S>0\), use the complete flow \(T_t\) of \(W\), and let \(\mathcal Lh=W\cdot\nabla h\) be its generator. Composition \(U_th=h\circ T_t\) is a strongly continuous unitary group on \(L^2(\gamma_d)\). For every polynomial \(h\) and \(k=1,2\), \(\partial_t^k(h\circ T_t)=(\mathcal L^kh)\circ T_t\), because \(DT_t\,W=W\circ T_t\); and \(\mathcal L(h\circ T_t)=(\mathcal Lh)\circ T_t\). The radial bound for \(W\), analogous to (12), gives Gaussian-integrable polynomial dominators on compact time intervals. Hence time derivatives may be taken under the integral, and cutoff integration by parts gives \(\langle\mathcal Lf,g\rangle=-\langle f,\mathcal Lg\rangle\) whenever \(f,g\) are of the form \(h\circ T_t\) with \(h\) polynomial. This is the only form of skew symmetry used below.

The mixed moment \(m(t)=\langle H_a,U_tH_b\rangle\) has
\(m''(t)=\langle H_a,U_t\mathcal L^2H_b\rangle=-\langle \mathcal LH_a,U_t\mathcal LH_b\rangle\). Expanding in the independent centred variables \(\eta_i=z_i^2-1\) gives the exact identity

\[
P_C(a):=\|\mathcal LH_a\|_2^2
=8\|Ca\|^2+16\sum_{i<j}C_{ij}^2(a_i-a_j)^2.                    \tag{16}
\]

Indeed,
\(\mathcal LH_a=2\sum_{i<j}C_{ij}(a_i-a_j)\eta_i\eta_j
-2\sum_i(Ca)_i\eta_i\); these Hermite components are orthogonal, with \(\mathbb E\eta_i^2=2\) and \(\mathbb E(\eta_i\eta_j)^2=4\). Since \(\|Ca\|^2=\|a\|^2S\) and \((a_i-a_j)^2\leq\|a\|^2\) for nonnegative spectra,
\(P_C(a)\leq24\|a\|^2S\); similarly \(P_C(b)\leq24\|b\|^2S\). Thus \(|m''(t)|\leq24\|a\|\|b\|S\).

Every coordinate reflection \(R\) satisfies \(W(Rz)=RW(z)\), because \(W_k\) is odd in \(z_k\) and even in every other coordinate. Uniqueness of solutions gives \(T_t(Rz)=RT_t(z)\), so every off-diagonal correlation \(\mathbb E[Z_iT_{t,j}(Z)]\), \(i\ne j\), vanishes. Write \(r_i(t)=\langle z_i,U_tz_i\rangle\). To bound its square sharply, use \(f_i=z_i\otimes z_i\), the unitary group \(U_t\otimes U_t\) on \(L^2(\gamma_d\otimes\gamma_d)\), and its generator \(J=\mathcal L\otimes I+I\otimes\mathcal L\), to which the same justification applies. Then

\[
|(r_i^2)''(t)|
=|\langle Jf_i,(U_t\otimes U_t)Jf_i\rangle|
\leq\|Jf_i\|_2^2=4\sum_j C_{ij}^2.                             \tag{17}
\]

Here \(\|Jf_i\|_2^2=2\|\mathcal Lz_i\|_2^2+2\langle z_i,\mathcal Lz_i\rangle^2\), with \(\|\mathcal Lz_i\|_2^2=2\sum_jC_{ij}^2\) and \(\langle z_i,\mathcal Lz_i\rangle=0\). The covariance part of \(G\) is \(2\sum_i a_ib_ir_i(t)^2\), so its second derivative in absolute value is at most

\[
8\sum_{i<j}C_{ij}^2(a_ib_i+a_jb_j)\leq8\|a\|\|b\|S.
\]

Consequently \(|G''(t)|\leq32\|a\|\|b\|S\) for every real \(t\), while \(G'(0)=4S\). Taylor's integral remainder gives
\(G(t)-G(0)\geq4St-16\|a\|\|b\|St^2\). At \(t=1/(8\|a\|\|b\|)\), this is at least \(S/(4\|a\|\|b\|)\). Multiplication by four through (3) proves the lower inequality. \(\square\)

The proof also gives a sharper computable certificate. Define

\[
\kappa_C=\sqrt{P_C(a)P_C(b)}
       +8\sum_{i<j}C_{ij}^2(a_ib_i+a_jb_j).
\]

For \(S>0\), \(D\geq32S^2/\kappa_C\), attained as a certified improvement at flow time \(4S/\kappa_C\). The bound (15) follows by the simpler estimate \(\kappa_C\leq32\|a\|\|b\|S\). In our numerical checks along the whole flow path, \(\max_t|G''(t)|\) reached up to 92% of \(\kappa_C\), so this curvature constant is not loose. The certificate is a one-step guarantee: running the same flow much longer than \(4S/\kappa_C\) can increase the objective again.

If \(\cos\theta=(a\cdot b)/(\|a\|\|b\|)\), (15) reads

\[
\|a\|\|b\|\sin^2\theta\leq D\leq8\|a\|\|b\|(1-\cos\theta).              \tag{18}
\]

Thus the normalised Gaussian excess \(D/(\|a\|\|b\|)\) lies between \(\theta^2(1+o(1))\) and \(4\theta^2(1+o(1))\) as \(\theta\to0\). The two bounds bracket \(D\) within a factor \(8/(1+\cos\theta)\), at most eight and tending to four near proportionality. We do not claim that these constants are optimal.

**Corollary 9 (the Gaussian value is within a factor of two).** For all Gaussian pairs of positive covariance rank,
\[
\tfrac12Q_G\leq Q_*\leq Q_G.
\]
The constant \(\tfrac12\) cannot be improved: if \(u\ne v\) in \(\mathbb R^d\) are nonincreasing with zero sum, and \(a=\mathbf 1+\epsilon u\), \(b=\mathbf 1+\epsilon v\), then \(d/(2(d+2))-O(\epsilon)\leq D/Q_G\leq\tfrac12\) as \(\epsilon\to0\). Hence \(\sup D/Q_G=\tfrac12\) over all Gaussian pairs.

**Proof.** Since \(\|a\|^2+\|b\|^2\geq2\|a\|\|b\|\), formula (4) gives \(Q_G\geq8\|a-b\|^2\geq16(\|a\|\|b\|-a\cdot b)\), so the known upper inequality in (15) gives \(D\leq Q_G/2\). For the sharpness family, equal traces make \(Q_G=8\|a-b\|^2=8\epsilon^2\|w\|^2\), where \(w=v-u\ne0\) has zero sum, and \(C_{ij}=\epsilon(w_i-w_j)+O(\epsilon^2)\). Lagrange's identity and the zero sums give \(S=\epsilon^2d\|w\|^2(1+O(\epsilon))\). In \(\kappa_C\), \(P_C(a)=8\|a\|^2S(1+O(\epsilon^2))\), \(\|a\|\|b\|=d+O(\epsilon^2)\), and the covariance term is \(16S(1+O(\epsilon))\). Hence \(32S^2/\kappa_C=4S/(d+2)\,(1+O(\epsilon))\), and dividing by \(Q_G\) gives the claim. \(\square\)

The upper half of Corollary 9 is a direct consequence of the bound of Delon et al.; the sharpness statement uses the new lower certificate. For nearly isotropic spectra, common after whitening or normalising representations, the Gaussian formula can therefore overstate the GW optimum by almost a factor of two. In fifty dimensions the certified excess is already 48% of \(Q_G\) (Section 6.3).

**Proposition 10 (comonotone refinement of the upper bound).** Let \(F_U^{-1}\) denote the quantile function of a real random variable \(U\), and let
\[
M_{\mathrm{co}}=\int_0^1F^{-1}_{\|X\|^2}(u)\,F^{-1}_{\|Y\|^2}(u)\,du,
\]
where \(\|X\|^2=\sum_ia_iZ_i^2\) and \(\|Y\|^2=\sum_ib_iZ_i^2\) are weighted chi-square variables. Then \(D\leq4(M_{\mathrm{co}}-\alpha\beta-2a\cdot b)\leq8(\|a\|\|b\|-a\cdot b)\).

**Proof.** In one dimension the comonotone coupling maximises \(\mathbb E[UV]\) among couplings of fixed marginals; this is the classical rearrangement inequality, equivalently the optimality of the monotone coupling in one-dimensional optimal transport (Villani, 2009). Hence \(\mathbb E\|X\|^2\|Y\|^2\leq M_{\mathrm{co}}\) under every coupling, and \(G(\pi)\leq M_{\mathrm{co}}+2a\cdot b\) by Lemma 2. Subtract \(G_G=\alpha\beta+4a\cdot b\) and use (3). Cauchy–Schwarz gives \(M_{\mathrm{co}}\leq\alpha\beta+2\|a\|\|b\|\), which is the second inequality. \(\square\)

\(M_{\mathrm{co}}\) is a one-dimensional integral. For the plane-to-line example, \(F^{-1}_{\chi^2_2}(u)=-2\log(1-u)\) and \(F^{-1}_{\chi^2_1}(u)=\Phi^{-1}((1+u)/2)^2\); high-precision quadrature gives \(M_{\mathrm{co}}=4.795281\ldots\), against the Cauchy–Schwarz value \(2+2\sqrt2=4.828427\ldots\). This refinement is not sharp either, because the comonotone norm coupling and the covariance bound are generally not attained by one coupling.

## 6.1 A pair-specific certificate

For a single field (10), put
\(P(x,y)=24(x^2+y^2)-32xy\) and

\[
\kappa_{ij}=\sqrt{P(a_i,a_j)P(b_i,b_j)}+8(a_ib_i+a_jb_j).              \tag{19}
\]

The same proof yields \(|G''(t)|\leq\kappa_{ij}\). At \(t=4\Delta_{ij}/\kappa_{ij}\),

\[
D\geq\frac{32\Delta_{ij}^2}{\kappa_{ij}}\quad(\Delta_{ij}\ne0).       \tag{20}
\]

The formula for \(P\) follows from
\(\mathbb E[u^4(v^2-1)^2]=6\) and
\(\mathbb E[u^2v^2(v^2-1)(1-u^2)]=-4\). Formula (17) gives the coefficient eight in (19); bounding the two differentiated factors separately would give a weaker constant. In two dimensions (20) coincides with the \(\kappa_C\) certificate. In high dimensions it is much weaker, so \(32S^2/\kappa_C\) is the recommended default.

## 6.2 Evaluating the correction

The simultaneous field has the rank-two form

\[
W_i(z)=z_i\{b_i[H_a(z)-\alpha]-a_i[H_b(z)-\beta]\},
\]

so each field evaluation costs \(O(d)\). Its trajectory from a fixed \(z\) can also be reconstructed from two scalar ordinary differential equations:

\[
\begin{aligned}
u'&=\sum_i a_i z_i^2\exp(2b_i u-2a_i v)-\alpha,\\
v'&=\sum_i b_i z_i^2\exp(2b_i u-2a_i v)-\beta,
\qquad u(0)=v(0)=0,\\
T_t(z)_i&=z_i\exp(b_i u(t)-a_i v(t)).
\end{aligned}
\]

Differentiating the last expression recovers \(T'_t=W(T_t)\), and completeness of that flow supplies global existence of this representation. Zero initial coordinates remain zero. For a full-rank source the improved coupling is a deterministic map \(x\mapsto D_b^{1/2}T_t(D_a^{-1/2}x)\), which can be applied to new points at \(O(d)\) cost per right-hand-side evaluation. The theorem concerns the exact flow; a finite-step numerical integrator requires its own marginal-error diagnostics and is not an exactly feasible coupling merely because it approximates the differential equation.

## 6.3 A practical bracket

Given two sorted, zero-padded eigenvalue lists \(a,b\), the following quantities need \(O(d^2)\) arithmetic and no optimisation:

1. the Gaussian value \(Q_G\) from (4) and \(\cos\theta=a\cdot b/(\|a\|\|b\|)\);
2. the certified lower bound on the excess, \(\max\{32S^2/\kappa_C,\ \max_{i<j}32\Delta_{ij}^2/\kappa_{ij}\}\);
3. the upper bound on the excess, \(\min\{8(\|a\|\|b\|-a\cdot b),\ Q_G/2\}\), optionally refined by Proposition 10.

Then \(Q_*\) lies in \([Q_G-\text{upper},\ Q_G-\text{lower}]\). The accompanying `gw_gaussian_bracket.py` implements items 1–3 with the standard library. Table 1 applies it to representative spectra. Gaps are shown as fractions of \(Q_G\); except for the PCA row, both spectra are scaled to unit trace.

| Spectra | \(\cos\theta\) | \(Q_G\) | certified excess at least | lower bound in (15) | excess at most (Prop. 10) |
|------------------------------------------|--------:|--------:|----------:|----------:|----------:|
| power law, \(d=50\): \(i^{-1}\) vs \(i^{-1.5}\) | 0.959 | 0.368 | 4.4% | 2.9% | 10.1% |
| exponential, \(d=50\): \(e^{-0.1i}\) vs \(e^{-0.2i}\) | 0.943 | 0.130 | 16.9% | 6.0% | 23.6% |
| PCA truncation: \(i^{-1}\), \(d=200\), vs its top 10, unscaled | 0.972 | 35.5 | 0.44% | 0.25% | 0.98% |
| near isotropic, \(d=50\): \(1\) vs \(1+U(0,0.2)\) | 0.999 | \(4.0\times10^{-4}\) | 48.0% | 12.5% | 50.0% |
| embedding dimensions 768 vs 512, both \(i^{-1}\) | 0.9998 | \(9.9\times10^{-4}\) | 2.4% | 1.3% | 5.3% |

: Certified brackets for the Gaussian excess \(D\) as a fraction of \(Q_G\). Values computed by `extension_checks.py`; the last column evaluates \(M_{\mathrm{co}}\) numerically.

Three readings follow. For typical power-law spectra the Gaussian formula is a reasonable proxy, and the bracket says how reasonable. For nearly isotropic spectra it can be off by almost a factor of two. When one side is an unscaled truncation of the other, the trace mismatch dominates \(Q_G\) and the relative error is small; the objective is not scale invariant, so normalisation conventions matter.

# 7. Worked examples

**Plane to line.** Take \(\mu=N(0,I_2)\), \(\nu=N(0,1)\). Use \(a=(1,1)\), \(b=(1,0)\). Every jointly Gaussian coupling has \(Q\geq12\), with equality for projection onto either coordinate. Let \((U_t,V_t)\) solve

\[
\dot U_t=U_t(V_t^2-1),\qquad
\dot V_t=V_t(1-U_t^2),\qquad (U_0,V_0)=(u,v),                  \tag{21}
\]

and couple \(X=(u,v)\sim N(0,I_2)\) to \(Y=U_{1/7}\). Lemma 6 proves that \(Y\sim N(0,1)\) exactly. Here \(\Delta=1\) and \(\kappa=8+8\sqrt6<28\). Therefore

\[
G(1/7)-G(0)\geq\frac47-\frac{\kappa}{98}>\frac27,
\qquad
Q(\pi_{1/7})<12-\frac87=\frac{76}{7}.                          \tag{22}
\]

The proved reduction in objective is greater than \(2/21\), or \(9.52\%\), of the best Gaussian objective. Optimising the curvature certificate gives
\(Q_*\leq12-4/(1+\sqrt6)\approx10.8404\). Table 2 collects what is known about \(Q_*\) in this example.

| Quantity | Value | Status |
|--------------------------------------------|---------:|------------------------------|
| Gaussian optimum \(Q_G\) | 12 | exact |
| upper bound \(76/7\), coupling (21) at \(t=1/7\) | 10.8571 | proved |
| upper bound \(12-4/(1+\sqrt6)\) | 10.8404 | proved |
| coupling (21) at \(t=1/7\) | 10.4697 | numerical quadrature |
| coupling (21) at its best time \(t\approx0.211\) | 10.2995 | numerical quadrature |
| rearrangement coupling (described below) | 10.0321 | numerical |
| lower bound, Proposition 10 | 8.8189 | proved inequality, integral evaluated numerically |
| lower bound, Delon et al. (upper inequality in (15)) | 8.6863 | proved |

: Plane-to-line example. Numerical rows use deterministic quadrature; the rearrangement row has marginal moment residuals of order \(10^{-14}\). None of these values is claimed to equal \(Q_*\).

The rearrangement coupling is \(Y=\Phi^{-1}(F_h(h(X)))\), where \(h(x)=x_1+\operatorname{sign}(x_1)(\lambda x_2^2+c\,x_1^2)\), \(F_h\) is the distribution function of \(h(X)\) and \(\Phi\) the standard normal one; it makes \(Y\) exactly standard Gaussian, and \((\lambda,c)\approx(0.367,0.172)\) was found by numerical optimisation. So \(Q_*\) lies between 8.68 and 10.85 with proof, is at least about 8.82 by Proposition 10 with its one-dimensional integral evaluated numerically (two independent high-precision evaluations agree to better than \(10^{-20}\), but no interval enclosure is claimed), and is numerically below about 10.03. The certified flow bound captures a little over half of the improvement that a simple rearrangement achieves numerically.

**Equal dimensions.** Take \(a=(2,1)\), \(b=(1,1)\), a planar Gaussian with unequal variances matched to a standard planar Gaussian. Peyré's line-versus-space construction does not apply. Here \(Q_G=12\), \(S=1\), and the pair certificate (20), which equals \(32S^2/\kappa_C\) in two dimensions, gives \(D\geq4/(3+\sqrt{14})\approx0.593\). The upper bounds give \(D\leq8(\sqrt{10}-3)\approx1.2982\), refined to \(D\lesssim1.2730\) by Proposition 10 evaluated numerically. Hence \(10.701\leq Q_*\leq11.407\) with proof, and \(Q_*\gtrsim10.727\) numerically: the Gaussian answer overstates the optimum by between about 4.9% and 10.6% of its value.

The accompanying code performs exact polynomial arithmetic with standard-library integer and rational operations. It checks the Gaussian divergence, radial growth identity, moment derivatives, generator norms and rational counterexample bound. It also enumerates all 441 pairs of sorted two-coordinate integer spectra with entries between zero and five, including zeros and repeated eigenvalues. A separately written checker verifies the simultaneous-flow identities, spectral-minor identity and gradient decomposition. A further script checks the curvature bound \(|G''(t)|\leq\kappa_C\) along whole flow paths in seven nonproportional cases, the factor-two corollary and its sharpness family, and Tables 1–2. All checks fail loudly, including under Python's `-O` flag, and a negative-control script confirms that planted errors are detected. These are algebraic, finite-range and numerical checks, not a proof-assistant verification of the theorems.

![A nonlinear Gaussian-preserving coupling improves on the Gaussian benchmark in the plane-to-line example. The horizontal axis is the flow time \(t\) in (21); the vertical axis is the objective (1). The diamond marks the proved upper bound at \(t=1/7\); the curve is numerical and does not assert the unrestricted optimum.](results/nonlinear_improvement.png)

# 8. Consequence for samples

Practitioners apply the Gaussian formula to sample covariances, or run a GW solver on the samples themselves. The following stability argument connects the population results to both.

Write \(\mathcal D(\mu,\nu)=Q_*(\mu,\nu)^{1/2}\) for measures with finite fourth moments, and \(W_4\) for the Wasserstein distance of order four.

**Lemma 11 (triangle inequality and stability).** For finite-fourth-moment measures \(\mu,\mu',\rho,\nu\) on Euclidean spaces,
\[
\mathcal D(\mu,\nu)\leq\mathcal D(\mu,\rho)+\mathcal D(\rho,\nu),\qquad
\mathcal D(\mu,\mu')^2\leq32\sqrt2\,W_4(\mu,\mu')^2\big(m_4(\mu)+m_4(\mu')\big)^{1/2},
\]
where \(m_4(\mu)=\int\|x\|^4d\mu\) and \(\mu,\mu'\) live in the same space.

**Proof.** The first inequality is the standard gluing argument for GW distances (Mémoli, 2011). Take optimal couplings of \((\mu,\rho)\) and \((\rho,\nu)\) and glue them along \(\rho\) into a law of \((X,R,Y)\). With independent copies, Minkowski's inequality in \(L^2\) gives
\(\|d_X^2-d_Y^2\|_2\leq\|d_X^2-d_R^2\|_2+\|d_R^2-d_Y^2\|_2\), where \(d_X^2=\|X-X'\|^2\) and similarly for the others; the left side is at least \(\mathcal D(\mu,\nu)\). For the second inequality, let \((X,\tilde X)\) be a \(W_4\)-optimal coupling, \(u=X-X'\) and \(v=\tilde X-\tilde X'\). Then \(|\|u\|^2-\|v\|^2|\leq\|u-v\|\,\|u+v\|\), so Cauchy–Schwarz gives \(Q\leq(\mathbb E\|u-v\|^4)^{1/2}(\mathbb E\|u+v\|^4)^{1/2}\). Using \((s+t)^4\leq8(s^4+t^4)\), \(\mathbb E\|u-v\|^4\leq16W_4^4\) and \(\mathbb E\|u+v\|^4\leq128(m_4(\mu)+m_4(\mu'))\). \(\square\)

**Proposition 12 (plug-in inconsistency).** Let \(\mu,\nu\) be Gaussian with positive ranks and nonproportional ordered spectra, and let \(\mu_n,\nu_m\) be empirical measures of independent samples. Almost surely, as \(n,m\to\infty\), the empirical GW value \(Q_*(\mu_n,\nu_m)\) converges to \(Q_*(\mu,\nu)\), while the plug-in value \(\hat Q_G\), formula (4) evaluated at the sample covariance eigenvalues, converges to \(Q_G\). Hence
\[
\lim(\hat Q_G-Q_*(\mu_n,\nu_m))=D\geq S/(\|a\|\|b\|)>0 .
\]

**Proof.** Almost surely the empirical measures converge weakly and their fourth moments converge, by Varadarajan's theorem and the strong law of large numbers; hence \(W_4(\mu_n,\mu)\to0\) (Villani, 2009, Chapter 6), and likewise for \(\nu_m\). Lemma 11 gives \(|\mathcal D(\mu_n,\nu_m)-\mathcal D(\mu,\nu)|\leq\mathcal D(\mu_n,\mu)+\mathcal D(\nu_m,\nu)\to0\). Sample covariances converge almost surely, sorted eigenvalues depend continuously on symmetric matrices, and (4) is continuous in the zero-padded spectra. Theorem 8 completes the proof. \(\square\)

The plug-in Gaussian formula is therefore an inconsistent estimator of the population GW value whenever the spectra are not proportional, with asymptotic bias between \(S/(\|a\|\|b\|)\) and the upper bound in (15). Proposition 12 concerns the global minimum over couplings of the empirical measures. Practical solvers return local solutions, so it does not describe what a particular solver finds, and it gives no rate of convergence.

# 9. Interpretation and remaining questions

The classification supplies a direct diagnostic for the Gaussian restriction. After reducing to positive covariance supports and sorting the spectra, proportionality is exactly the condition for unrestricted Gaussian optimality. Because exact proportionality almost never holds for estimated spectra, the useful diagnostic is quantitative: \(\cos\theta\) together with the bracket of Section 6.3.

The affine theorem identifies the source of the obstruction: the radial fourth-moment term creates a polynomial curl that an anisotropic affine map cannot remove on interior support. Read for linear encoders, it says that a principal-component projection that discards directions, or keeps them with unequal weights, is never a GW-optimal coupling between a distribution with interior support and its own image. Using squared-distance GW to score linear embeddings, as in GW-based dimension reduction (Van Assel et al., 2025), therefore always leaves slack that is not evidence of a poor embedding. The four-point example shows that finite clouds can escape this conclusion.

The theorem also explains why the conclusion does not transfer to inner-product GW, which lacks the radial term; for spheres the radial term is constant, consistent with the Gaussian subproblem of Arya et al. (2026). In the Gaussian setting, replace \(G\) in (2) by \(w\mathbb E[\|X\|^2\|Y\|^2]+\eta\|K\|_F^2\), with \(w>0\) and \(\eta\geq0\). The existence of a jointly Gaussian optimiser is still equivalent to proportional positive spectra: at the common-latent coupling the covariance derivative remains zero, while (11) is multiplied by \(w\), and proportional spectra attain both separate upper bounds. We assert only this existence statement for the family. The all-optimiser statement of Theorem 5 is not claimed, and it fails when \(\eta=0\): for \(X\sim N(0,1)\), the coupling \(Y=X\) on \(\{|X|\leq1\}\) and \(Y=-X\) otherwise has \(Y^2=X^2\) and maximises the radial term without being affine or Gaussian. For \(w=0\), the covariance bound alone admits Gaussian optima for all spectra.

All quantitative statements assume Gaussian marginals. Applied to non-Gaussian data through their covariances alone, formula (4) need not equal or bound the GW value between the actual marginals, because the fourth moments in (3) differ from the Gaussian ones, and the bracket of Section 6.3 does not apply.

The paper does not provide the unrestricted nonproportional optimiser, a closed formula for its value, or an algorithm guaranteed to find it. Open questions include: the exact plane-to-line optimiser and whether it is a Monge map (cf. Dumont et al., 2025); the optimal constants in (15), in particular the limit of \(D/S\) as \(S\to0\); whether iterating the flow, or using it to initialise a solver, approaches \(Q_*\); analogues for entropy-regularised GW and for GW with ordinary distances; and a version of Theorem 3 for finite clouds in general position. Moment–sum-of-squares relaxations (Tran & Soh, 2026) offer a route to numerical lower bounds, and the certified brackets above provide a test for such computations.

# References

Arya, S., Auddy, A., Clark, R. A., Lim, S., Mémoli, F., & Packer, D. (2026). The Gromov–Wasserstein distance between spheres. *Foundations of Computational Mathematics, 26*(1), 75–130. <https://doi.org/10.1007/s10208-024-09678-3>

Dandapanthula, S., Podkopaev, A., Kasiviswanathan, S. P., Ramdas, A., & Goldfeld, Z. (2026). *Optimal transportation and alignment between Gaussian measures* (Version 2) [Preprint]. arXiv:2512.03579. <https://arxiv.org/abs/2512.03579>

Delon, J., Desolneux, A., & Salmona, A. (2022). Gromov–Wasserstein distances between Gaussian distributions. *Journal of Applied Probability, 59*(4), 1178–1198. <https://doi.org/10.1017/jpr.2022.16>. Theorem and equation numbers refer to the preprint arXiv:2104.07970v1 (author order Salmona, Delon, Desolneux).

Dumont, T., Lacombe, T., & Vialard, F.-X. (2025). On the existence of Monge maps for the Gromov–Wasserstein problem. *Foundations of Computational Mathematics, 25*(2), 463–510. <https://doi.org/10.1007/s10208-024-09643-0>. Section numbers refer to arXiv:2210.11945v3.

Flamary, R., Courty, N., Gramfort, A., et al. (2021). POT: Python Optimal Transport. *Journal of Machine Learning Research, 22*(78), 1–8. <https://jmlr.org/papers/v22/20-451.html>

Le, K., Le, D. Q., Nguyen, H., Do, D., Pham, T., & Ho, N. (2022). Entropic Gromov–Wasserstein between Gaussian distributions. In *Proceedings of the 39th International Conference on Machine Learning* (PMLR Vol. 162, pp. 12164–12203). <https://proceedings.mlr.press/v162/le22a.html>

Mémoli, F. (2011). Gromov–Wasserstein distances and the metric approach to object matching. *Foundations of Computational Mathematics, 11*(4), 417–487. <https://doi.org/10.1007/s10208-011-9093-5>

Morel, G., Drumetz, L., Benaïchouche, S., Courty, N., & Rousseau, F. (2023). Turning normalizing flows into Monge maps with geodesic Gaussian preserving flows. *Transactions on Machine Learning Research*. <https://openreview.net/forum?id=2UQv8L1Cv9>; arXiv:2209.10873v4.

Peyré, G. (2026a). *Optimal transport for machine learners* [Book, online version retrieved 31 October 2026]. <https://www.gpeyre.com/ot4ml/OT4ML/OT4ML.pdf>. Propositions 13.33–13.34.

Peyré, G. (2026b, 23 October). *Gaussians are not closed for Gromov–Wasserstein* [Blog post]. <https://www.gpeyre.com/blog/2026/10/23/gaussian-gromov-wasserstein-counterexample/>

Salmona, A. (2023). *Transport of probability distributions across different Euclidean spaces* [Doctoral thesis, Université Paris-Saclay]. <https://doi.org/10.70675/a8e13508z4193z4651za1f4z4f4199ba110c>; HAL tel-04793664.

Tran, H. A., & Soh, Y. S. (2026). *Quadratic optimization over probability measures with coupling constraints* [Preprint]. arXiv:2608.24184. <https://arxiv.org/abs/2608.24184>

Van Assel, H., Vincent-Cuaz, C., Courty, N., Flamary, R., Frossard, P., & Vayer, T. (2025). *Distributional reduction: Unifying dimensionality reduction and clustering with Gromov–Wasserstein* (Version 3) [Preprint]. arXiv:2402.02239. <https://arxiv.org/abs/2402.02239>

Vayer, T. (2020). *A contribution to optimal transport on incomparable spaces* [Doctoral thesis; preprint]. arXiv:2011.04447. <https://arxiv.org/abs/2011.04447>

Villani, C. (2009). *Optimal transport: Old and new*. Grundlehren der mathematischen Wissenschaften, Vol. 338. Springer. <https://doi.org/10.1007/978-3-540-71050-9>

Zhang, Z., Goldfeld, Z., Mroueh, Y., & Sriperumbudur, B. K. (2024). Gromov–Wasserstein distances: Entropic regularization, duality and sample complexity. *The Annals of Statistics, 52*(4), 1616–1645. <https://doi.org/10.1214/24-AOS2406>
