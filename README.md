# Can the standard designs see harm?

Simulation code and supporting documents for a methods paper on the evidence about prenatal antidepressant exposure and autism. The paper asks two questions:

1. What do the influential studies adjust for, and which of those choices are consistent with a causal diagram of the problem? (adjustment-set audit)
2. If prenatal antidepressant exposure did raise autism risk by a plausible amount, how often would the standard designs report reassurance? (simulation, calibrated to the published crude associations)

Author: Vera Wilde, Heilbronn Data Science Center, Technical University of Munich (ORCID 0000-0002-6470-4109).

## Contents

- `R/` — simulation code, versioned. The latest data-generating model is `10-27-26_ad-pregnancy-sim_v0.4.R`; earlier versions are kept so pilot results can be reproduced.
  - v0.1: base model, estimators M1–M8, two-target calibration
  - v0.2: maternal–paternal use link and the paternal control adjusted for maternal use (M8b)
  - v0.3: persistent within-family non-adherence
  - v0.4: autism subtypes (with / without intellectual disability) with separate effects, dose heterogeneity, age at diagnosis with administrative censoring and Cox estimators, severity-dependent exposure misclassification, and causal truths computed from potential-outcome risks in the target-trial population
- `R/*run-grid*` — scenario grids; `R/*summarize-grid*` — summaries and figure.
- `docs/` — simulation protocol (ADEMP structure) and the adjustment-set audit.
- `results/` — pilot summaries (50 replicates per scenario) and the final results.
  - `11-12-26_final-results_v1.0.csv`: preregistered final run, 74 scenarios with 1,000 replicates each (one row per scenario and method).
  - `final-replicates/`: the stored replicate-level estimates from the final run, one file per scenario.
- `R/explore/` and `results/explore/` — exploratory analyses added after the final run, outside the preregistration (see below).

## Running

From the `R/` folder:

```
Rscript 10-27-26_run-grid_v0.4.R    # v0.4 pilot grid, 50 replicates per scenario
Rscript 10-27-26_run-final_v1.0.R   # final run: 74 scenarios x 1,000 replicates
```

Requires R (≥ 4.2) with the `survival` package. The final run uses all cores but one and saves each scenario as it finishes; if interrupted, rerun the same command and finished scenarios are skipped.

## Status

The final run (74 scenarios, 1,000 replicates each) was preregistered on OSF before it began: https://osf.io/9wfq8/. The code used for it is frozen at commit de8ae14.

## Exploratory analyses (not preregistered)

These were added after the final run and are reported as exploratory.

- `R/explore/11-13-26_posthoc-pooled-sibling-check.R`: how often a single simulated sibling comparison returns an estimate as low as the pooled sibling estimate in Chan et al. (2026), 0.98, or an interval ending at or below 1.12. Output: `results/11-13-26_posthoc-pooled-sibling-check.csv`.
- `R/explore/11-13-26_simulated-meta-analysis.R`: random-effects meta-analyses of seven simulated sibling studies built from the stored replicates. Output: `results/explore/11-13-26_simulated-meta-analysis_results.csv`.
- `R/explore/11-12-26_ad-pregnancy-sim_v0.4-explore.R`: a copy of the v0.4 model with one added parameter (`assort`), through which fathers' depression can share the mother's family component of severity. The frozen v0.4 file is unchanged.
- `R/explore/11-12-26_run-explore-nc.R`: twelve uncalibrated negative-control scenarios, 150 replicates each.
- `R/explore/11-12-26_run-explore-dose.R`: the dose scenario under stronger pregnancy-specific confounding or a lower between-pregnancy correlation of severity, plus a no-effect scenario, 300 replicates each.
- `R/explore/11-14-26_summarize-explore.R`: summaries of the two sets of runs above, in `results/explore/`.

The two run scripts are run from `R/explore/`; the others from the repository root.

## License

Code: MIT. Documents in `docs/`: CC BY 4.0.
