# Plan: Robust SFG Coating Design

## Problem

Trial 6 (N=16, 260306_021817) hits all nominal targets but is fragile to
manufacturing perturbations:
- T@1064: median 6510 ppm vs target 5000 (30% over)
- R@700: 95th percentile 11,500 ppm vs target 200 (57x over)
- R@2050: 95th percentile 850 ppm vs target 200 (4x over)

Root cause: first derivatives dT/dx approx 0 but curvature d2T/dx2 is huge
(ratios 155-753x). The design sits at an extremum that collapses under
0.5% manufacturing perturbations.

## What has been done

### Curvature penalties implemented (costs.py)

Added d2TdnH, d2TdnL, d2Tdd cost terms. Second-derivative penalties
computed from the same finite-difference stencil as the first-derivative
terms (zero extra transfer matrix calls beyond one shared T0 evaluation).
Opt-in via YAML weights. All tests pass.

### Landscape analysis (this session)

1. Cost landscape mapping: 500 random samples at each perturbation radius.
   At 0.5% perturbation cost jumps from 1.08 to median 837 (100% feasible).
   At 5% median 4.3M (40% feasible). At 10% median 45B (3% feasible).
   The feasible valley is extremely narrow in 33D.

2. Basin search: 200 random starts + L-BFGS-B local optimization.
   60/200 found feasible solutions (30% success rate).
   Trial 6 basin is still the best (cost 373 vs next-best 1163).
   All feasible basins have similar curvature (d2TdnL around 0.8-1.2).
   Conclusion: the curvature is structural, not basin-specific.

3. Npairs scan with curvature (quick runs, maxiter=2000):
   N=16: T1=3489 ppm, MC 95th=3788 (tight!) but R@700=16000 ppm (bad)
   N=18: T1=4332, MC 95th=4768, R@700=2787 (still bad)
   N=20: T1=5619, MC 95th=6237
   N=22: T1=5519, MC 95th=6584, d2Tdd drops to 1.87

4. N=16 with curvature, 3 seeds, maxiter=4000:
   Curvature penalties dramatically tighten MC spread (7% vs 60%).
   But optimizer trades R@700 for robustness: R@700 goes to 2000-16000 ppm.
   Trans2 weight (250) loses to curvature weights in multiplicative cost.

### Key insight

At N=16, the 3-wavelength constraint (1064 + 700 + 2050) plus robustness
is over-constrained. The 700 nm target requires precise interference
cancellation far from the design wavelength, and that cancellation is
inherently fragile. More layers give more degrees of freedom.

## Next steps

### Step 1: Broader Npairs sweep with curvature (N=16..24)

Run optimalbragg sweep with curvature penalties enabled, wider range.
Use maxiter=4000 per point. Key question: at what N does R@700 become
simultaneously achievable AND robust?

Save as new trial params: _trial10_params.yml with curvature costs.

### Step 2: Cost weight tuning

The multiplicative cost C = prod(1 + w_i * c_i) means weights interact
nonlinearly. Current weights may be poorly balanced with curvature added.
Specific experiments:
- Increase Trans2 weight from 250 to 1000-5000 to prevent R@700 collapse
- Decrease d2Tdd weight (currently most dominant curvature at around 2.7)
- Try ratio: Trans2_weight / d2TdnL_weight > 100

### Step 3: Evaluate optimizer adequacy

Current evidence suggests dual_annealing is finding the right basins.
The issue is the cost landscape topology, not the search algorithm.
However, if the curvature-augmented landscape has many more local minima
(plausible since we added 3 competing objectives), consider:

- CMA-ES (covariance matrix adaptation): good for 30-50D problems with
  correlated variables and multimodal landscapes. pip install cma.
  Natural fit because it adapts its search covariance to the local landscape
  shape, exactly what we need for a narrow curved valley.
- Bayesian optimization (scikit-optimize): too expensive at 33D. Skip.
- Multi-start L-BFGS-B with Latin hypercube: cheap, embarrassingly
  parallel. Basin search showed 30% feasible rate. 50 starts would
  give around 15 feasible minima to pick from. Good complement to DA.

Only pursue CMA-ES if Step 1-2 fail to produce acceptable robustness.

### Step 4: Validation

For the best candidate:
1. Full MC (5000 samples) with emcee
2. All 3 wavelength 95th percentiles must be within 2x of target
3. Generate report + PDF for sharing

## Files

- OptimalBragg/costs.py: curvature penalties (DONE)
- projects/SFG/_trial10_params.yml: new trial with curvature weights (TODO)
- projects/SFG/PLAN_robust_sfg.md: this file
