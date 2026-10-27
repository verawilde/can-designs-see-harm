source("10-27-26_ad-pregnancy-sim_v0.3.R")
suppressPackageStartupMessages(library(parallel))
out <- "grid_v03"; dir.create(out, showWarnings = FALSE)
settings <- list(
  nonadh_0.20        = list(p_nonadh = 0.20),
  nonadh_0.35        = list(p_nonadh = 0.35),
  nonadh_0.35_rho0.3 = list(p_nonadh = 0.35, rho_severity = 0.3),
  paternal_1.2_couple = list(rr_paternal = 1.2, couple_logor = log(3.42))
)
grid <- list(); for (nm in names(settings)) for (rr in c(1.0, 1.2, 1.5)) grid[[sprintf("%s__rr%.1f", nm, rr)]] <- c(settings[[nm]], list(rr_drug = rr))
invisible(mclapply(names(grid), function(k) {
  f <- file.path(out, paste0(k, ".rds")); if (file.exists(f)) return(k)
  r <- run_scenario(do.call(make_scenario, grid[[k]]), n_rep = 50, seed0 = 50000 + 1000 * match(k, names(grid)))
  r$key <- k; saveRDS(r, f)
  cat(format(Sys.time(), "%H:%M:%S"), "done", k, "\n", file = file.path(out, "progress.log"), append = TRUE); k
}, mc.cores = 2, mc.preschedule = FALSE))
cat("GRID COMPLETE\n", file = file.path(out, "progress.log"), append = TRUE)
