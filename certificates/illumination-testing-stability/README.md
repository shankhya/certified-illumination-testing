# Continuous illumination-test count: strengthened stability verification

Place this folder beside the unchanged `illumination-testing` core folder.
Only Python 3.10+ and its standard library are required:

    python -I -O -B verify_stability.py --output verification.json
    python -I -O -B check_faults.py --output fault-check.json

No optimizer, MATLAB, network or historical project file is used. The earlier
`illumination-testing-robust` and `illumination-testing-spd` extensions are not
needed. Original CIE source attribution, licenses and checksums are supplied
by the core's DATA_NOTICES.md and data directory; no source data are modified.
This is an internal research artifact; no public release is authorized.

The verifier checks nine continuous duals, five exact spectral witnesses and
all 31 constituent subsets with fewer than five tests. It reconstructs the
CIE kernels and computes absolute polynomial integrals with exact rational
enclosures. It never accepts a saved research PASS flag as mathematical input.

At g=0.0355/nm, pair acceptance tolerance 0.0001 and mismatch limit 0.02,
five constituent tests are necessary and sufficient throughout either region:

    mu/0.000073935629 + rho/0.005633164098 <= 1;
    kappa/0.000103842397 + rho/0.005633164098 <= 1.

All parameters are nonnegative. In the first region, every actual normalized
test and target kernel is within mu in L1 of its nominal kernel. In the second,
every actual illuminant is a positive gain times E+DeltaE, the fixed observer
is used, |DeltaE|<=kappa*E continuously, and actual white-Y normalization is
used. In both, actual reflectances are passive and each is within rho in
uniform norm of the nominal slope class. The exact rational regions in
proof.json are slightly larger than the conservative printed versions.

The test library is A, D65, FL2, FL7, FL11 and deployment is FL3, FL8, FL12,
on 380–780nm with CIE1931 2-degree CMFs and white Y=1. Nominal SPD and CMF
tables are separately linearly interpolated and their products integrated.
The count is for constituent selection, not arbitrary mixture exposures.
These are sufficient budgets, not optimal robustness or physical calibration
requirements. No perceptual-error equivalence or empirical prediction accuracy
is claimed. Known-operator prediction remains a separate interpretation.

For one saved dual, define u=q-sum_i v_i*m_i, U=integral|u|, N=sum_i|v_i|,
R=integral|u+p'| and P>=integral|p|, with zero endpoint flux. Put
base=R+2*g*P+delta*N. The first region uses

    base+(1+N)*mu+2*rho*U <= tau.

For each test light, combine its three signed multiplier kernels into s_l.
Let S=integral(q)+sum_l integral|s_l|. The SPD region uses

    base+kappa*(delta*N+S)+2*rho*U <= tau*(1-kappa).

The actual target bound is the left side divided by 1-kappa. This retains
common-channel SPD errors and white-denominator changes. The model term uses
U, not R. For a nominal exact witness d, its selected actual test response is
at most kappa*integral(m_i*|d|)/(1-kappa); a target of nominal magnitude t
remains at least (t-kappa*integral(q*|d|))/(1+kappa). All such integrals are
recomputed. The same witness for a given subset works for every allowed
perturbation in the stated region.

check_faults.py checks five analytic controls and rejects seven intentionally
invalid claims, including understated residual norms, excessive allowances and
an incorrect normalization assumption. Successful checking establishes these
mathematical statements under the declared model, not literature priority or
journal acceptance.
