source("10-27-26_ad-pregnancy-sim_v0.2.R")
suppressPackageStartupMessages(library(parallel))
dir.create("grid_couple", showWarnings = FALSE)
keys <- c(1.0, 1.2, 1.5)
invisible(mclapply(keys, function(rr) {
  sc <- make_scenario(rr_drug = rr, couple_logor = log(3.42))
  r <- run_scenario(sc, n_rep = 50, seed0 = 90000 + rr * 1000)
  saveRDS(r, sprintf("grid_couple/couple__rr%.1f.rds", rr)); rr
}, mc.cores = 2, mc.preschedule = FALSE))
for (rr in keys) { r <- readRDS(sprintf("grid_couple/couple__rr%.1f.rds", rr))
  cat("\nRR", rr, "compatible", r$cal$compatible, "crude", round(mean(r$raw$crude_rr),2), "pat", round(mean(r$raw$crude_rr_pat),2), "\n")
  print(format(r$perf$by_method[, c("method","mean_est_RR","coverage","p_reassuring")], digits=2), row.names=FALSE) }
