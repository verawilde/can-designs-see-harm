# Summaries of the exploratory (not preregistered) runs.
# Run from the repository root: Rscript R/explore/11-14-26_summarize-explore.R
summ <- function(d, truth = NULL) {
  do.call(rbind, lapply(split(d, d$method), function(s) {
    s <- s[is.finite(s$est) & is.finite(s$se), ]
    lo <- s$est - 1.96 * s$se; hi <- s$est + 1.96 * s$se
    data.frame(method = s$method[1], n = nrow(s), mean_RR = exp(mean(s$est)),
               p_lower_above1 = mean(lo > 0), p_includes1 = mean(lo <= 0 & hi >= 0))
  }))
}
nc <- readRDS("results/explore/explore_nc.rds")
o1 <- do.call(rbind, lapply(names(nc), function(k) cbind(scenario = k, summ(nc[[k]]))))
write.csv(o1, "results/explore/11-14-26_explore-negative-control_summary.csv", row.names = FALSE)
dz <- readRDS("results/explore/explore_dose.rds")
o2 <- do.call(rbind, lapply(names(dz), function(k) {
  r <- dz[[k]]; s <- summ(r$raw); t <- r$cal$truths
  s$truth_RR <- ifelse(grepl("high", s$method), exp(t["high"]), exp(t["all"]))
  cbind(scenario = k, s)
}))
write.csv(o2, "results/explore/11-14-26_explore-dose_summary.csv", row.names = FALSE)
