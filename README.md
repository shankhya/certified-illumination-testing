# Certified illumination testing

Associated data and computational certificates for **Certified illumination-test counts: limits imposed by spectral regularity and measurement tolerance**, by Shankhya Debnath.

This private repository accompanies an author-review draft. It has not been submitted, peer reviewed, or publicly released. The nine-page draft uses the author-supplied Optica universal single-column template, including references. Section 7 is a Discussion; there is no standalone reproducibility section. Its AI-generated narrative is not submission-compliant under Optica's current AI-text policy; see `AUTHOR_REVIEW.md`. No claim of journal acceptance or exhaustive literature priority is made.

## Main result and exact scope

For the test constituents **A, D65, FL2, FL7, FL11**, deployment **FL3, FL8, FL12**, CIE 1931 two-degree observer, 380--780 nm, individual reflectance slope at most **0.0355/nm**, and normalized XYZ infinity mismatch limit **0.02**, exactly **five tests** are necessary and sufficient for pair acceptance tolerance from zero through **0.000181293988**. Five exact witness pairs exclude all **31** smaller selections, even with exact nominal tests.

The continuous full-library precision--regularity boundary is enclosed, with maximum unresolved slope width **0.001104210242/nm** on tolerance `[0,0.002]`. The declared `0.001/nm` refinement goal remains unmet. The robust count regions and calibrated affine prediction corollary are also verified.

These are universal results for the declared continuous mathematical class. They are not dataset accuracy scores, perceptual color differences, arbitrary-mixture test-exposure counts, or hardware specifications. See the manuscript for all hypotheses.

## Verify the mathematical claims

Python **3.10 or later**, standard library only:

```console
python -I -O -B verify_all.py
```

This checks frozen-file hashes, reconstructs the CIE kernels, verifies all 90 continuous dual certificates and 32 exact witnesses, checks the continuous frontier and all-smaller-protocol exclusions, verifies both robust regions and the calibrated predictor, and rejects 16 deliberately invalid extension claims. It writes reports to `verification/`. No optimizer, network, MATLAB, historical project checkout, or saved PASS flag is needed. Typical runtime is several minutes; exact rational arithmetic is intentionally used.

The four directories in `certificates/` are frozen portable bundles. Their original manifests and historical README wording are retained byte for byte. In particular, the core's older frontier-width number is superseded by the sibling frontier extension. The top-level driver runs the strengthened combination. Internal certificate identifiers record provenance; they are not manuscript section numbers.

## Independent MATLAB checks and figures

The scripts use paths relative to their own location. They do not change global startup files or call an optimizer. From the repository root:

```matlab
addpath('matlab');
if ~isfolder('verification'), mkdir('verification'); end
check_170;
check_171;
make_figures;
```

MATLAB R2025b was used for the checked release. The first two scripts independently rebuild kernels from the original CIE CSV files and integrate absolute polynomials using roots. These floating crosschecks support but do not replace the exact verifier. `make_figures` regenerates the three vector PDFs and PNG previews from the stored numerical records. It runs no scientific search. AI-assisted plotting code is documented in `AUTHOR_REVIEW.md`.

`matlab/affine_predict.m` evaluates a supplied 9-by-15 matrix. Exact coefficients and an exact standard-library CLI are in `certificates/illumination-testing-prediction/`. Its README declares the observation and white-calibration information required. An arbitrary input file cannot itself establish the physical calibration assumptions.

## Contents

| Path | Purpose |
|---|---|
| `paper/` | LaTeX source, official template files, figures, compiled author-review PDF |
| `certificates/illumination-testing/` | Original CIE tables/metadata, 81 duals, 32 witnesses, core verifier |
| `certificates/illumination-testing-frontier/` | Nine saved-dual blends and strengthened continuous frontier |
| `certificates/illumination-testing-stability/` | Residual-aware kernel/model and shared-SPD count regions |
| `certificates/illumination-testing-prediction/` | Exact affine coefficients, calibration bounds, prediction CLI |
| `matlab/` | Independent checks, figure code, input fixtures |
| `figures/` | Vector figures, previews, generation report |
| `verification/` | Recomputed machine-readable verification results |
| `CLAIMS.md` | Manuscript-to-evidence map and important limitations |
| `MANIFEST.sha256` | Snapshot integrity, excluding itself and generated logs |

## Build the paper

Use a standard LaTeX installation or Overleaf with `paper/main.tex`. Only one manuscript `.tex` file is used. The official `optica-article` class, required style files, and vector figures are supplied. With a local installation:

```console
cd paper
latexmk -pdf main.tex
```

Alternatively, Tectonic 0.17.0 can compile the source in untrusted mode. The local research build uses a separate support cache; that executable/cache is not redistributed. No font-size or margin reduction was used to meet the requested ten-page cap.

## Data, rights, and provenance

Original CIE source tables and metadata are unchanged, with their source DOI, checksums, and CC BY-SA 4.0 notices. See `certificates/illumination-testing/DATA_NOTICES.md` and `RIGHTS.md`. No third-party full-text article, external optimizer, credential, machine configuration, or unrelated manuscript is included.

The research evidence was frozen before this publication package. The private parent research commit is recorded in `PROVENANCE.json`; this repository does not contain the parent's private working history. The verification package checks the retained numerical evidence independently of that history. Source-generation history is distinct from certificate verification: the original optimization campaign is not reproduced by `verify_all.py`.

Prepared 29 October 2026. Before any journal submission, the author must prepare the narrative in accordance with journal policy, review the scientific claims and outstanding literature gaps, and arrange appropriate reviewer access or public release of the data. No automatic release or submission is configured.
