# EXPLORATORY, not preregistered. Uncalibrated negative-control comparison:
# confounding held fixed (no recalibration) while the maternal effect, the
# maternal-side confounding and the link between the parents' use vary.
source("11-12-26_ad-pregnancy-sim_v0.4-explore.R")
methods <- methods[c("M1", "M7", "M8", "M8b")]
nrep <- as.integer(Sys.getenv("N_REP", "200"))
cb <- 1.03125   # maternal-side confounding calibrated in the base case (true RR 1.0)
grid <- expand.grid(rr = c(1.0, 1.5), conf = c(0, cb), link = c("none", "couple", "assort"),
                    stringsAsFactors = FALSE)
out <- list()
for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]
  sc <- make_scenario(rr_drug = g$rr,
                      couple_logor = if (g$link == "couple") log(3.42) else 0,
                      assort = if (g$link == "assort") 0.4 else 0)
  cal <- list(conf_scale = g$conf, conf_pat = 0)
  res <- do.call(rbind, lapply(seq_len(nrep), function(r) { x <- one_rep(sc, cal, 5000 + r); x$rep <- r; x }))
  key <- sprintf("nc_rr%.1f_conf%s_%s", g$rr, ifelse(g$conf > 0, "on", "off"), g$link)
  out[[key]] <- res
  cat(key, format(Sys.time()), "\n"); saveRDS(out, "explore_nc.rds")
}
