# Scenario grid, one factor at a time around the base case, at true RR 1.0, 1.2, 1.5.
# Pilot grid: 50 replicates per scenario (protocol target is 1,000 for the final run).
source("10-26-26_ad-pregnancy-sim_v0.1.R")
suppressPackageStartupMessages(library(parallel))

n_rep <- as.integer(Sys.getenv("N_REP", "50"))
out_dir <- "grid"; dir.create(out_dir, showWarnings = FALSE)

settings <- list(
  base            = list(),
  rho_0.9         = list(rho_severity = 0.9),
  rho_0.3         = list(rho_severity = 0.3),
  capture_1.0     = list(capture = 1.0),
  capture_0.70    = list(capture = 0.70),
  halfsib_0       = list(half_sib = 0),
  halfsib_0.30    = list(half_sib = 0.30),
  nonshared_0     = list(beta_specific = 0),
  nonshared_0.25  = list(beta_specific = 0.25),
  paternal_1.2    = list(rr_paternal = 1.2),
  surveil_1.15    = list(surveillance = 1.15),
  brown_size      = list(n_fam = 7300)   # ~620 discordant pairs, as in Brown 2017
)
grid <- list()
for (nm in names(settings)) for (rr in c(1.0, 1.2, 1.5))
  grid[[sprintf("%s__rr%.1f", nm, rr)]] <- c(settings[[nm]], list(rr_drug = rr))

run_one <- function(key) {
  f <- file.path(out_dir, paste0(key, ".rds"))
  if (file.exists(f)) return(key)
  sc <- do.call(make_scenario, grid[[key]])
  seed0 <- 10000 + 1000 * match(key, names(grid))
  r <- run_scenario(sc, n_rep = n_rep, seed0 = seed0)
  r$key <- key
  saveRDS(r, f)
  cat(format(Sys.time(), "%H:%M:%S"), "done", key, "\n", file = file.path(out_dir, "progress.log"), append = TRUE)
  key
}
invisible(mclapply(names(grid), run_one, mc.cores = 2, mc.preschedule = FALSE))
cat("GRID COMPLETE\n", file = file.path(out_dir, "progress.log"), append = TRUE)
