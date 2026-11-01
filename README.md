# Affine rigidity and Gaussian exactness in squared-distance Gromov–Wasserstein transport

Evidence Press research release, 1 November 2026. Version 0.1.0-candidate. DOI: [10.5281/zenodo.23850758](https://doi.org/10.5281/zenodo.23850758) (concept DOI, all versions: [10.5281/zenodo.23850757](https://doi.org/10.5281/zenodo.23850757)). Release page: <https://evidencepress.org/releases/gaussian-gromov-wasserstein-exactness/>. Author: Anonymous. UK English. Text, data and diagrams CC0-1.0; code MIT.

## In plain words

Gromov–Wasserstein (GW) transport compares two datasets or shapes that may live in different spaces, such as a cloud of points in the plane and one on a line, by matching their internal pairwise distances. It is costly to compute, so practitioners like closed-form answers. For Gaussian data, Delon, Desolneux and Salmona (2022) found one by restricting attention to Gaussian couplings, and conjectured that nothing is lost by the restriction.

This paper shows exactly when nothing is lost: the Gaussian answer is the true optimum **if and only if** the two covariances have the same rank and proportional ordered eigenvalues. In every other case, including pairs of the same dimension, an explicit nonlinear coupling that keeps both distributions exactly Gaussian does strictly better. The size of the error is bracketed from both sides, and the Gaussian answer is never more than twice the true optimum.

A concrete example: matching a standard Gaussian plane to a standard Gaussian line, the Gaussian answer is 12; the paper constructs a coupling with objective provably below 76/7 ≈ 10.86 (about 10.47 numerically), and shows the true optimum lies between 8.68 and 10.85.

## Main results

For the unregularised objective \(Q(\pi)=\mathbb E[(\|X-X'\|^2-\|Y-Y'\|^2)^2]\):

- **Theorem 5 (exactness).** Two positive-rank Gaussians admit a jointly Gaussian unrestricted optimiser iff their positive ranks agree and ordered positive eigenvalues are proportional. Sufficiency is due to Delon et al. (2022); necessity is the proposed contribution.
- **Theorem 3 (affine rigidity).** If the source has finite fourth moment and support with nonempty relative interior in its affine hull, an affine map is optimal onto its image iff it is a similarity or a constant.
- **Theorem 8 (two-sided error).** With zero-padded eigenvalue lists \(a,b\) and \(S=\|a\|^2\|b\|^2-(a\cdot b)^2\), \(S/(\|a\|\|b\|)\le Q_G-Q_*\le8(\|a\|\|b\|-a\cdot b)\). The upper bound is due to Delon et al.; the lower bound is new.
- **Corollary 9.** \(Q_G/2\le Q_*\le Q_G\), and the constant 1/2 is approached by nearly isotropic spectra in high dimension.
- **Proposition 10, Proposition 12.** A comonotone refinement of the known upper bound; the plug-in Gaussian formula is an inconsistent estimator of the population GW value for nonproportional spectra.

Not done: the unrestricted optimiser and value for nonproportional spectra; entropic or ordinary-distance GW; finite point clouds (Theorem 3 needs interior support, and a four-point example shows why).

## Files

| File | Contents |
|---|---|
| `affine_rigidity_gaussian_gw.pdf` | The paper (16 pages), complete proofs |
| `manuscript.md` | Source of the paper; `build_pdf.sh` and `paper_header.tex` rebuild the PDF |
| `proof_outline.md` | Proof architecture and dependencies |
| `target.yaml` | Frozen objective, hypotheses and claim fingerprints |
| `novelty_report.md` | Formula-first prior-art search record and its limits |
| `REVIEW_SUMMARY.md` | Internal AI-assisted review: objections and how each was handled |
| `VERIFICATION.md` | What each script checks, how to run it, and what it does not establish |
| `verify_gw.py`, `independent_symbolic_audit.py`, `verify_combined_flow.py`, `proof_digest_checks.py` | Exact and finite checks of the identities and constants |
| `extension_checks.py` | Whole-path curvature checks, factor-two corollary and sharpness family, Tables 1–2 |
| `gw_gaussian_bracket.py` | Standard-library recipe: Gaussian value and certified bracket for \(Q_*\) from two spectra |
| `negative_controls.py` | Plants a false identity in each checker and confirms it is caught, with and without `python -O` |
| `results/`, `*.json` | Recorded outputs and the figure |

## Reproduce

Python 3.12 with the pinned packages in `requirements.txt`:

```bash
uv run --python 3.12 --with-requirements requirements.txt bash build.sh
```

`build.sh` runs every checker (about seven minutes, mostly `extension_checks.py`) and then `build_pdf.sh`, which needs Pandoc and XeLaTeX. The exact checkers also run alone with the standard library, for example `python3 -O verify_gw.py --exact-only --out OUTDIR`. All checks raise on failure, including under `-O`.

## Status

The proofs, code and reviews were produced in an AI-assisted research process and checked by internal AI-assisted review (five reviewer roles, an adversarial review and a proof audit by a separate model, a clean replay and a re-review). This is not external human peer review or a formal proof-assistant verification. Historical novelty rests on the bounded search recorded in `novelty_report.md`.

## Release records

[STATUS.md](STATUS.md), [ASSURANCE.md](ASSURANCE.md), [PROVENANCE.md](PROVENANCE.md), [CLAIMS.json](CLAIMS.json), [RESEARCH_GATES.json](RESEARCH_GATES.json), [AI_INDEX.md](AI_INDEX.md), [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md), [CITATION.cff](CITATION.cff), [LICENSES.md](LICENSES.md), `refs.bib` with registry records in `crossref/`, the replay comparer `tools/replay_compare.py` (`python3 -I -B tools/replay_compare.py . [--quick]`) and `SHA256SUMS`.
