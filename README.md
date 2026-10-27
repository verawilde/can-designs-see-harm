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
- `results/` — pilot summaries (50 replicates per scenario; Monte Carlo SEs are reported with every performance measure). `10-27-26_grid-results_v0.4.csv` holds the v0.4 pilot.

## Running

From the `R/` folder:

```
Rscript 10-27-26_run-grid_v0.4.R    # v0.4 pilot grid, 50 replicates per scenario
Rscript 10-27-26_run-final_v1.0.R   # final run: 74 scenarios x 1,000 replicates
```

Requires R (≥ 4.2) with the `survival` package. The final run uses all cores but one and saves each scenario as it finishes; if interrupted, rerun the same command and finished scenarios are skipped.

## Status

Pilot results only. The final run (1,000 replicates per scenario) is preregistered on OSF: https://osf.io/9wfq8/

## License

Code: MIT. Documents in `docs/`: CC BY 4.0.
