# Refined continuous tolerance–regularity frontier

Keep this folder beside the unchanged `illumination-testing` core. Run:

    python -I -O -B verify_refinement.py --check-faults --output verification.json

Python 3.10+ and its standard library suffice. The checker uses original CIE
data in the core; it requires no optimizer, MATLAB, network or historical
research files. Data attribution and source licenses remain in the core's
DATA_NOTICES.md. This is an internal artifact, not an authorized public release.

All 81 original continuous duals, nine new blended duals, 32 exact witnesses
and 45 strict frontier points are checked. Each new dual is reconstructed from
two previously checked parents and its exact rational weight. The same nominal
continuous model, five-test library, three deployment illuminants, coordinate
normalization and tau=0.02 apply. No new reflectance or illumination assumption
is introduced by blending.

The certified maximum unresolved critical-slope gap over pair acceptance
tolerance [0,0.002] is at most 0.001104210242/nm. The original goal of
0.001/nm remains unmet. This extension preserves that shortfall explicitly.
The lower failure envelope is unchanged. The original sharp five-test count
at g=0.0355 and positive tolerance is also rechecked.

Four deliberate faults must be rejected: altered blend weight with retained
bound, mixed target coordinates, understated whole-interval gap, and a false
precision-success flag. The exact comparison checks every affine-envelope
segment, not just plotted sample points.

The stability extension is separate and is not required here. This package
does not establish an optimal frontier, arbitrary-mixture test count, measured
perceptual accuracy or historical novelty of dual-certificate combinations.
