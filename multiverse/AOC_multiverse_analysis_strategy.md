## AOC Multiverse Analysis Strategy

- Primary inference uses trial-level mixed-effects models.
- Subject-level models are retained as aggregation robustness checks.
- Trial-level models should include subject random slopes for `Condition` where identifiable, with fallback to random intercepts if convergence or singularity requires simplification.
- Subject-level models should not add a `trial` random effect because trial variance is already collapsed during aggregation.
- FOOOF fit quality is evaluated with threshold sensitivity analyses (`R2` grid), not a single hard cutoff only.

## Minimum Reporting Set

- Trial-level estimate and CI for each key hypothesis.
- Subject-level estimate and CI for the same hypothesis.
- Threshold-sensitivity trajectory across at least `R2 >= 0.50, 0.60, 0.70, 0.80` (optionally `0.90` strict check).
- Retention profile at each threshold: retained rows, subjects, trials, and universes.
