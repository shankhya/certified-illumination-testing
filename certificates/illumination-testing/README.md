# Continuous illumination-testing certificate bundle

This internal research bundle verifies the strengthened five-light results
without an optimizer, MATLAB, network access, or the historical project.
It requires Python 3.10 or later and only the standard library. From this
folder run:

    python -I verify.py --output verification.json
    python -I check_faults.py --output fault-checks.json

The output files are created only after the corresponding checks pass.
Assertions are explicit exceptions and remain enabled with Python -O.
The exporter used to assemble proof.json lives in the parent research project;
it is not required to verify the bundle. No saved PASS flag is trusted.

## Precisely what is verified

* All 81 continuous dual certificates underlying the final tolerance–slope
  envelope, with exact rational source reconstruction and residual costs.
* 32 distinct piecewise-linear spectral differences, including 45 strict
  frontier points and five exact four-light counterexamples.
* Safe and strict-failure envelopes for every pairwise tolerance in [0,0.002],
  using exact affine intersections, not a wavelength or parameter sample claim.
* At individual slope g=0.0355/nm, exactly five constituent tests are necessary
  and sufficient for mismatch at most 0.02, for all pairwise tolerances from
  zero through 0.000181293988. All 31 smaller constituent protocols are covered.
* The maximum remaining enclosing gap 0.001345896087/nm (rounded upward).
  The predeclared resolution goal 0.001/nm was missed.

The independent MATLAB crosschecks, original optimization generation,
three/four-test transition claims, passive-class exact-metamer result,
dataset audit and perturbation-budget calculations are outside this bundle's
coverage. They remain in the parent research records. This checker is a
standalone implementation, not an independent research team or proof assistant.

## Mathematical contract

The domain is [380,780] nm. Linearly interpolate each tabulated illuminant and
each CIE 1931 two-degree colour-matching function separately, then multiply.
Normalize each light by the exact integral against the interpolated y-bar.
The response kernels are thus continuous and piecewise quadratic. XYZ weights
are all one. Tests are A, D65, FL2, FL7 and FL11; targets are FL3, FL8 and FL12.
These assumptions fix a nominal continuous model, not unknown physical
between-sample behaviour. The error metric is normalized XYZ infinity norm,
not a perceptual colour-difference formula.

For two passive reflectances with individual Lipschitz constant at most g,
their difference d satisfies |d|<=1 and |d'|<=2g almost everywhere. Conversely
every such d is realized by r=(1+d)/2 and s=(1-d)/2. Tests constrain |Md|<=delta
coordinatewise. For target row q, any real coefficient vector v and absolutely
continuous p with zero endpoints imply

    |integral q d| <= delta * sum|v| + integral|q-vM+p'| + 2g * integral|p|.

This follows by integration by parts and triangle inequalities. No dual
attainment or numerical optimality assumption is needed. The supplied flux
is a piecewise cubic with matching endpoints. Its derivative cancels the
nonconstant residual within each cell. The checker computes the constant
residual exactly and encloses integral|p| using nonnegative Bernstein basis
functions and exact dyadic subdivision. It uses a tighter integration target
than the original research verifier before accepting the stored upper bound.

For each target, maximize the safe slope across valid duals; take the minimum
across targets. For failure bounds, each supplied piecewise-linear d is
integrated exactly, with amplitude and actual slope checked. Sign-oriented
strict witnesses in the same target coordinate may be convexly combined.
Their test tolerance and slope are bounded by the corresponding convex
combinations, and their strict target failure persists. Witnesses from
different target coordinates cannot be interpolated by this argument.

The minimum-count proof uses one exact zero-tolerance witness for each of the
five four-light subsets. A witness also defeats each contained smaller subset.
These witnesses belong to the declared slope class. Full-library sufficiency
at the stated positive tolerance follows from nine checked dual bounds.
Monotonicity extends both conclusions to the entire stated tolerance interval.

For independently admissible additive measurement errors bounded by epsilon
in each test coordinate, a common observation exists exactly when pairwise
responses differ by at most 2*epsilon. Unrestricted infinity-norm optimal
prediction error is half the pairwise mismatch bound. This classical
optimal-recovery interpretation is supporting theory, not a novelty claim.

## Limitations and integrity

Counts refer to selections of the five named constituent lights. They do not
exclude fewer arbitrary mixtures, other illuminants, other response sensors,
or data-dependent protocols. The class is defined mathematically; membership
of an empirical dataset does not imply it models every real surface.
The bundle establishes certificate validity and stated bounds. It does not
prove the exact optimal frontier inside its unresolved bands or literature
priority. Worthey (1988) remains a full-text comparison gap.

Source CSVs and original metadata are unchanged. See DATA_NOTICES.md.
proof.json records exact proof objects and hashes of historical candidate
files for provenance. The source hashes are checked during verification.
MANIFEST.sha256 lists the complete deliverable's bytes, excluding generated
reports and the manifest itself. The package has not been publicly released
or submitted to a journal.
