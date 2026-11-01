# Verification record

Release date: 1 November 2026. Target: *Affine rigidity and Gaussian exactness in squared-distance Gromov–Wasserstein transport*.

## Scope

The theorems are proved analytically in the paper. The scripts below check identities, constants, finite cases and numerical illustrations. They are not a proof-assistant formalisation, and numerical quadrature is not interval arithmetic: values labelled "numerical" in the paper are not certified enclosures.

## Environment

Python 3.12 (CPython), NumPy 2.3.5, SciPy 1.17.0, Matplotlib 3.10.8, mpmath 1.4.1 (`requirements.txt`). PDF: Pandoc and XeLaTeX with xdvipdfmx. The exact checkers need only the standard library.

## Checks

| Script | What it checks | Result |
|---|---|---|
| `verify_gw.py --exact-only --out DIR` | 24 exact rational checks: Gaussian divergence of the pair field, radial identity, moment derivatives, generator norms, the rational certificate \(Q<76/7\), and all 441 sorted two-coordinate integer spectrum pairs with entries 0–5 | pass |
| `verify_gw.py --out results` | The above plus Gauss–Hermite quadrature (orders 64–384) and DOP853 integration of the plane-to-line flow; figure and `results/flow_quadrature.csv` | pass; objective about 10.4697 at \(t=1/7\) |
| `independent_symbolic_audit.py` | 28 identities in a separately written rational-polynomial implementation, including symbolic spectral coefficients in dimension four and the gradient decomposition (14) | pass |
| `verify_combined_flow.py --out DIR` | Simultaneous-field identities on 120 integer spectra in dimensions 2–4, by exact-degree tensor quadrature | pass |
| `proof_digest_checks.py --out DIR` | Six balanced sign couplings for the four-point example, 85 conditional-mean cases, 60 full-versus-reduced ODE comparisons | pass |
| `extension_checks.py --out DIR` | (a) factor-two corollary on 1,156 integer pairs and 20,000 random pairs; (b) comonotone integral \(M_{\mathrm{co}}\) by two independent mpmath evaluations; (c) plane-to-line bracket (Table 2) with deterministic quadrature and marginal-moment residuals; (d) equal-dimension example; (e) Table 1; (f) \(|G''(t)|\le\kappa_C\) and the Taylor bound along whole flow paths in seven nonproportional cases (worst ratio 0.921); (g) sharpness family of Corollary 9 for \(d\in\{2,5,50,500\}\) | pass (about 6 minutes, peak 1.3 GB) |
| `gw_gaussian_bracket.py` | Demonstration of the Section 6.3 recipe | runs |
| `negative_controls.py --out DIR` | Plants one false value in each of six checkers and confirms each run fails, under `python3 -I` and `python3 -I -O` | 12/12 detected |

Every checker raises `SystemExit` on failure; none relies on `assert`, so the checks remain active under `python -O`. Status fields are written only after the checks run.

## Numerical values quoted in the paper

| Quantity | Value | Status |
|---|---:|---|
| \(M_{\mathrm{co}}\), plane to line | 4.79528098495884 | mpmath, two methods agree to better than \(10^{-20}\); not an interval enclosure |
| Lower bound on \(Q_*\) from Proposition 10 | 8.81888 | as above |
| Coupling (21) at \(t=1/7\) | 10.46966 | numerical quadrature |
| Coupling (21) at best time \(t\approx0.2110\) | 10.29950 | composite Gauss–Legendre, about \(10^{-9}\) accuracy |
| Rearrangement coupling, \((\lambda,c)\approx(0.3668,0.1722)\) | 10.03206 | deterministic quadrature; \(\mathbb E Y^2-1\), \(\mathbb E Y^4-3\) about \(10^{-14}\) |

Proved values: \(76/7\), \(12-4/(1+\sqrt6)\approx10.8404\), \(12-8(\sqrt2-1)\approx8.6863\).

## Replay

The original handoff bundle was replayed from a fresh extraction before revision; outputs matched the recorded files except run-time fields. The revised package is replayed again from its frozen archive during publication, under `python -O -B`.
