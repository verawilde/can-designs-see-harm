# Final run: all 75 preregistered scenarios, 1,000 replicates each (N_REP to override).
# Run from the folder that holds this file:  Rscript 10-27-26_run-final_v1.0.R
# Uses all cores but one. Each scenario is saved as it finishes, so the run can be
# stopped and restarted: finished scenarios are skipped.
here <- tryCatch(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))),
                 error = function(e) getwd())
source(file.path(here, "10-27-26_ad-pregnancy-sim_v0.4.R"))
suppressPackageStartupMessages(library(parallel))
n_rep <- as.integer(Sys.getenv("N_REP", "1000"))
cores <- max(1, detectCores() - 1)
out_dir <- file.path(here, "final"); dir.create(out_dir, showWarnings = FALSE)
reg <- list(fu_min = 3, fu_max = 15)
one_factor <- list(
  base = list(), rho_0.9 = list(rho_severity = 0.9), rho_0.3 = list(rho_severity = 0.3),
  capture_1.0 = list(capture = 1.0), capture_0.70 = list(capture = 0.70),
  halfsib_0 = list(half_sib = 0), halfsib_0.30 = list(half_sib = 0.30),
  nonshared_0 = list(beta_specific = 0), nonshared_0.25 = list(beta_specific = 0.25),
  paternal_1.2 = list(rr_paternal = 1.2), surveil_1.15 = list(surveillance = 1.15),
  brown_size = list(n_fam = 7300),
  couple = list(couple_logor = log(3.42)),
  paternal_1.2_couple = list(rr_paternal = 1.2, couple_logor = log(3.42)),
  nonadh_0.20 = list(p_nonadh = 0.20), nonadh_0.35 = list(p_nonadh = 0.35),
  nonadh_0.35_rho0.3 = list(p_nonadh = 0.35, rho_severity = 0.3),
  censor = reg, miscsev = list(misc_sev = 0.5),
  target_crude_1.5 = list(target_crude = 1.5), target_crude_1.9 = list(target_crude = 1.9),
  target_pat_1.16 = list(target_paternal = 1.16), target_pat_1.40 = list(target_paternal = 1.40))
grid <- list()
for (nm in names(one_factor)) for (rr in c(1.0, 1.2, 1.5))
  grid[[sprintf("%s__rr%.1f", nm, rr)]] <- c(one_factor[[nm]], list(rr_drug = rr))
grid$subtype_mod       <- list(rr_noid = 1.3, rr_id = 1.0)
grid$subtype_strong    <- list(rr_noid = 1.6, rr_id = 1.0)
grid$dose              <- list(rr_drug = 1.05, rr_high_extra = 1.7)
grid$combined__null    <- c(reg, list(rr_drug = 1.0, misc_sev = 0.5))
grid$combined__subtype <- c(reg, list(rr_noid = 1.3, rr_id = 1.0, misc_sev = 0.5))
stopifnot(length(grid) == 74)
run_one <- function(key) {
  f <- file.path(out_dir, paste0(key, ".rds")); if (file.exists(f)) return(key)
  sc <- do.call(make_scenario, grid[[key]])
  r <- run_scenario(sc, n_rep = n_rep, seed0 = 50000 + 1000 * match(key, names(grid)))
  r$key <- key; r$session <- sessionInfo(); saveRDS(r, f)
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "done", key, "\n", file = file.path(out_dir, "progress.log"), append = TRUE)
  key
}
cat("Running", length(grid), "scenarios x", n_rep, "replicates on", cores, "cores\n")
invisible(mclapply(names(grid), run_one, mc.cores = cores, mc.preschedule = FALSE))
cat("FINAL RUN COMPLETE\n", file = file.path(out_dir, "progress.log"), append = TRUE)
