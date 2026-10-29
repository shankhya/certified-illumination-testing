# Manuscript-to-evidence map

The manuscript uses the fixed observer, illuminants, spectral interval, interpolation rule, white-Y normalization, coordinate weights, and mismatch threshold stated in Section 2. Changing these changes the problem.

| Manuscript result | Checkable evidence | What it does not prove |
|---|---|---|
| Theorem 2: five constituents necessary and sufficient | Core `minimumCountClaim`, nine selected duals, five exact four-light witnesses, `verify_count`; frontier driver rechecks | Minimum among arbitrary mixture exposures or arbitrary lamps |
| Table 1: all four-light exclusions | Core witness objects referenced by `excludedFourLightWitnessIds`; exact continuous integration | Typical behavior of measured materials |
| Critical-slope endpoints and complete frontier enclosure | Core and frontier `safePieces`, `failurePieces`, exact intersections, 45 strict points | Exact boundary within the unresolved band |
| Kernel/model uncertainty region | Stability proof and verifier; independent MATLAB `check_170` | Physically optimal robustness radius |
| Shared-SPD/model uncertainty region | Stability proof and verifier; independent MATLAB `check_171` | Observer uncertainty or leakage at nominal SPD zeros |
| Explicit calibrated affine predictor | Prediction `proof.json`, `verify_prediction.py`, exact `predict.py` | Exact minimax optimality of the supplied matrix, or physical-spectrum reconstruction |
| Classical prediction/ambiguity identity | Self-contained proof in the paper; classical optimal-recovery references | Euclidean/perceptual identities or arbitrary noise geometry |
| Nonnegative deployment mixtures | White-normalized convex-combination identity in the paper | Lower counts for mixtures used as tests |
| Deterministic adaptive prediction lower count | Central-observation proof using the same subset witnesses | Perturbed noiseless adaptive pair testing or randomized expected-error performance |

The certificate package establishes correctness of these mathematical bounds under its assumptions. It does not establish literature priority or SOTA performance. It does not treat an inability to certify a point as proof that the physical requirement fails.

The figure source fixtures are floating exports for independent crosschecks and plotting. The rational certificate files, rather than those plotting fixtures, determine the rigorous claims. Reported upper bounds are rounded outward upward; lower bounds are rounded downward; uncertainty-axis denominators are rounded inward.

The uniform frontier target `0.001/nm` was missed: the verified remaining width is at most `0.001104210242/nm`. The fixed-slope tolerance bracket also missed its tighter resolution goal. Both limitations remain explicit. No further optimization was used for this manuscript package.
