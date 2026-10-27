# v0.4 pilot grid: heterogeneity by subtype and dose, censoring, severity-dependent
# misclassification, and a combined scenario. 50 replicates each (pilot).
source("10-27-26_ad-pregnancy-sim_v0.4.R")
suppressPackageStartupMessages(library(parallel))
n_rep <- as.integer(Sys.getenv("N_REP", "50"))
out_dir <- "grid_v04"; dir.create(out_dir, showWarnings = FALSE)
reg <- list(fu_min = 3, fu_max = 15)
grid <- list(
  subtype_mod__null   = list(rr_drug = 1.0),
  subtype_mod         = list(rr_noid = 1.3, rr_id = 1.0),
  subtype_strong      = list(rr_noid = 1.6, rr_id = 1.0),
  dose                = list(rr_drug = 1.05, rr_high_extra = 1.7),
  censor__rr1.0       = c(reg, list(rr_drug = 1.0)),
  censor__rr1.2       = c(reg, list(rr_drug = 1.2)),
  censor__rr1.5       = c(reg, list(rr_drug = 1.5)),
  miscsev__rr1.0      = list(rr_drug = 1.0, misc_sev = 0.5),
  miscsev__rr1.2      = list(rr_drug = 1.2, misc_sev = 0.5),
  miscsev__rr1.5      = list(rr_drug = 1.5, misc_sev = 0.5),
  combined__null      = c(reg, list(rr_drug = 1.0, misc_sev = 0.5)),
  combined__subtype   = c(reg, list(rr_noid = 1.3, rr_id = 1.0, misc_sev = 0.5))
)
run_one <- function(key) {
  f <- file.path(out_dir, paste0(key, ".rds")); if (file.exists(f)) return(key)
  sc <- do.call(make_scenario, grid[[key]])
  r <- run_scenario(sc, n_rep = n_rep, seed0 = 50000 + 1000 * match(key, names(grid)))
  r$key <- key; saveRDS(r, f)
  cat(format(Sys.time(), "%H:%M:%S"), "done", key, "\n", file = file.path(out_dir, "progress.log"), append = TRUE)
  key
}
invisible(mclapply(names(grid), run_one, mc.cores = 2, mc.preschedule = FALSE))
cat("GRID COMPLETE\n", file = file.path(out_dir, "progress.log"), append = TRUE)
