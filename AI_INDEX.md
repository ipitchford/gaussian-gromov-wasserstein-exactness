# Agent navigation

## Identity

- Release slug: `gaussian-gromov-wasserstein-exactness`. Version `0.1.0-candidate`, tag `v0.1.0-candidate`. DOI 10.5281/zenodo.23850758; concept DOI 10.5281/zenodo.23850757.
- Repository: https://github.com/ipitchford/gaussian-gromov-wasserstein-exactness. Release page: https://evidencepress.org/releases/gaussian-gromov-wasserstein-exactness/
- Scholarly creator: Anonymous. Status: unrefereed candidate ([STATUS.md](STATUS.md)).
- Manuscript: [PDF](affine_rigidity_gaussian_gw.pdf) (16 pages), Markdown source [manuscript.md](manuscript.md), generated LaTeX [affine_rigidity_gaussian_gw.tex](affine_rigidity_gaussian_gw.tex), builder [build_pdf.sh](build_pdf.sh) with [paper_header.tex](paper_header.tex); bibliography [refs.bib](refs.bib) with registry records in [crossref/](crossref/).
- Machine-readable claims: [CLAIMS.json](CLAIMS.json). Frozen target and claim fingerprints: [target.yaml](target.yaml). Overview: [README.md](README.md). Proof architecture: [proof_outline.md](proof_outline.md).
- Assurance limits: [ASSURANCE.md](ASSURANCE.md). Review record (rating-free): [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md). Provenance: [PROVENANCE.md](PROVENANCE.md).
- Prior art: [RESEARCH_GATES.json](RESEARCH_GATES.json), [novelty_report.md](novelty_report.md).
- Citation: [CITATION.cff](CITATION.cff), [.zenodo.json](.zenodo.json). Licences: [LICENSES.md](LICENSES.md), [LICENSE](LICENSE), [LICENSE-MIT.txt](LICENSE-MIT.txt), [PUBLIC_DOMAIN.md](PUBLIC_DOMAIN.md).

## Scientific map

Object: the squared-distance Gromov–Wasserstein objective Q(pi) = E[(|X-X'|^2 - |Y-Y'|^2)^2] between Gaussian measures on Euclidean spaces of possibly different dimensions; Q_G is the minimum over jointly Gaussian couplings (closed form of Delon, Desolneux and Salmona, 2022), Q_* the unrestricted minimum. Question: when is Q_G = Q_*, and how large is Q_G - Q_* otherwise?

- Section 1: GW explained, the Gaussian shortcut, the plane-to-line example, results in plain words, known versus new.
- Section 2: Lemma 1 (moment reduction), Lemma 2 (Gaussian benchmark and covariance bound); both known.
- Section 3: Theorem 3 (affine rigidity: affine optimal onto its image iff similarity or constant, interior support); Corollary 4 (local descent under a density assumption); the four-point example showing why interior support is needed.
- Section 4: Theorem 5 (exactness iff equal positive rank and proportional ordered eigenvalues); necessity is new.
- Section 5: Lemma 6 (admissible pair flow), Proposition 7 (maximal Gaussian-preserving slope).
- Section 6: Theorem 8 (two-sided excess bound), Corollary 9 (factor two, sharp), Proposition 10 (comonotone refinement), Sections 6.1–6.3 (certificates, evaluation, practical bracket recipe), Table 1.
- Section 7: plane-to-line coupling below 76/7, Table 2 bracket [8.68, 10.85]; equal-dimension example.
- Section 8: Lemma 11, Proposition 12 (plug-in inconsistency). Section 9: interpretation and open questions.

## Evidence and replay

- Script guide and what each check does not establish: [VERIFICATION.md](VERIFICATION.md). One-command build: [build.sh](build.sh). Comparer: [tools/replay_compare.py](tools/replay_compare.py) (temporary copy; every checker with -I -O -B; regenerated outputs compared within a relative tolerance of 1e-7; `--quick` omits the two slow numerical parts). Environment: [requirements.txt](requirements.txt), [ENVIRONMENT.txt](ENVIRONMENT.txt). CI: [.github/workflows/replay.yml](.github/workflows/replay.yml). Record: [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md). Integrity: [SHA256SUMS](SHA256SUMS).
- Scripts and outputs: [verify_gw.py](verify_gw.py) -> [results/verification_results.json](results/verification_results.json), [results/run_log.json](results/run_log.json), [results/flow_quadrature.csv](results/flow_quadrature.csv), [results/nonlinear_improvement.png](results/nonlinear_improvement.png); [independent_symbolic_audit.py](independent_symbolic_audit.py) -> [independent_symbolic_audit.json](independent_symbolic_audit.json); [verify_combined_flow.py](verify_combined_flow.py) -> [combined_flow_verification.json](combined_flow_verification.json); [proof_digest_checks.py](proof_digest_checks.py) -> [proof_digest_checks.json](proof_digest_checks.json); [extension_checks.py](extension_checks.py) -> [extension_checks.json](extension_checks.json); [gw_gaussian_bracket.py](gw_gaussian_bracket.py) (standard-library bracket recipe) -> [gw_gaussian_bracket_demo.log](gw_gaussian_bracket_demo.log); [negative_controls.py](negative_controls.py) -> [negative_controls.json](negative_controls.json) (12/12 planted errors detected with and without -O).

## Trust boundary

All review was internal and AI-assisted. There was no peer review, no independent reproduction and no formal verification. Exact checks verify identities and finite cases; numerical values (for example 10.47, 10.03 and about 8.82 in the plane-to-line example) are quadrature or high-precision evaluations, not interval enclosures. Sufficiency, the Gaussian-restricted value and the upper excess bound are due to Delon et al. (2022); Peyré (2026) gave the first counterexample. The unrestricted optimum for nonproportional spectra is not computed. Priority holds only within the documented bounded searches.
