# EXPLORATORY, not preregistered. Dose effect under stronger pregnancy-specific
# confounding or a low between-pregnancy correlation of severity; calibrated as
# in the final run.
source("11-12-26_ad-pregnancy-sim_v0.4-explore.R")
methods <- methods[c("M1", "M2", "M2high", "M7", "M7high", "M8")]
nrep <- as.integer(Sys.getenv("N_REP", "300"))
grid <- list(
  dose_nonshared0.25 = make_scenario(rr_drug = 1.05, rr_high_extra = 1.7, beta_specific = 0.25),
  dose_rho0.3        = make_scenario(rr_drug = 1.05, rr_high_extra = 1.7, rho_severity = 0.3),
  null_nonshared0.25 = make_scenario(rr_drug = 1.0, beta_specific = 0.25))
out <- list()
for (k in names(grid)) {
  r <- run_scenario(grid[[k]], nrep, seed0 = 7000)
  r$perf <- NULL; out[[k]] <- r
  cat(k, format(Sys.time()), r$cal$compatible, r$cal$note, "\n"); saveRDS(out, "explore_dose.rds")
}
