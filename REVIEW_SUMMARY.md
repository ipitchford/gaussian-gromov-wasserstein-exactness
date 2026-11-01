# Internal review summary (AI-assisted; not external peer review)

Dates: 31 October – 1 November 2026. Routes: five reviewer roles (editor-in-chief, methodology/proof, domain, perspective, devil's advocate) run as independent Claude Code subagents following the academic-paper-reviewer protocol; an adversarial review and, after revision, an audit of the new proofs by a separate model (OpenAI Codex "Astra", via a local companion tool); a separate proof audit of the revision; a clean replay of the handoff bundle from fresh extraction; and a re-review of the revision against every round-one objection. This is internal checking, not independent reproduction or human peer review.

**Outcome.** No reviewer or audit found a critical or major error in the mathematics. Several reviewers recomputed all load-bearing identities and constants independently, symbolically and numerically. The substantive objections concerned the verification code, positioning against the literature, readability and practical meaning. All were handled in one revision round.

| Objection | Response |
|---|---|
| The checking scripts used bare `assert`, so they passed vacuously under `python -O`; one 441-case grid recorded "passed" unconditionally | All asserts replaced by explicit checks; status written only after checks; outputs to `--out`; a negative-control script shows planted errors are caught with and without `-O` |
| The conjecture of Delon et al. was understated as "wider exactness" | Now stated as their conjecture that the Gaussian-restricted plan solves the unrestricted problem; Theorem 5 shows it holds exactly on proportional spectra and fails everywhere else, including equal dimensions |
| Necessity is essentially a routine first-variation computation | Agreed and stated: for Gaussians, Lemma 6 alone gives necessity; Theorem 3 is presented as the general mechanism |
| Missing literature: linearisation and duality (Zhang et al. 2024), Vayer's radial condition, Salmona's thesis, the OT4ML book, recent numerical work | Added from Crossref, arXiv and HAL records; Theorem 3 positioned against Vayer's sufficient condition |
| No plain-language explanation of GW, no roadmap, example too late | New introduction explaining GW and the Gaussian shortcut, a plain-language results box, roadmap, and the plane-to-line example in Section 1 |
| The size of the error relative to the Gaussian value was never stated | New Corollary 9: \(Q_G/2\le Q_*\le Q_G\), with sharpness; Table 1 of certified relative brackets for realistic spectra; Section 6.3 recipe and `gw_gaussian_bracket.py` |
| How far below 12 is the plane-to-line optimum? | Table 2: proved bracket [8.68, 10.85], about 8.82 by a numerically evaluated comonotone refinement (Proposition 10), numerical couplings down to about 10.03 |
| Practitioners use sample covariances | New Lemma 11 and Proposition 12: the plug-in Gaussian formula is inconsistent for nonproportional spectra |
| An equal-dimension example would show the regime not covered by earlier counterexamples | Added: \(a=(2,1)\), \(b=(1,1)\), \(10.70\le Q_*\le11.41\) |
| Curvature bound only checked at time zero | `extension_checks.py` checks it along whole flow paths |
| Generator-domain justification in Theorem 8 compressed; density step in Theorem 5; equal-mass detail in Theorem 3 | Written out |
| The weighted-objective remark could be read as covering all optimisers | Restricted to existence, with a counterexample to the stronger reading |
| Decimal bounds from numerically evaluated integrals were called proved; some rounded inward | Relabelled as numerical; all proved endpoints now rounded outward |
| Notation clashes (\(A\) and \(K\) each with two meanings) | Curvature constants renamed \(\kappa\); eigenvalue norms written \(\|a\|,\|b\|\) |

Not done (recorded as open): the unrestricted optimiser for nonproportional spectra; optimal constants; finite-sample rates; entropic and ordinary-distance GW; a finite-cloud version of Theorem 3; a solver experiment on samples.
