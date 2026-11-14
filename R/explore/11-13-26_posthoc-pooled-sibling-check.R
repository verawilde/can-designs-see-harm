# Exploratory (not preregistered): how often a single simulated sibling comparison (M7)
# returns an estimate as low as the pooled sibling estimate in Chan et al. 2026 (0.98),
# and how often its 95% interval ends at or below that estimate's upper bound (1.12).
# Run from the repository root: Rscript R/explore/11-13-26_posthoc-pooled-sibling-check.R
dir <- "results/final-replicates"
files <- list.files(dir, pattern = "\\.rds$")
out <- do.call(rbind, lapply(files, function(f) {
  x <- readRDS(file.path(dir, f))$raw
  x <- x[x$method == "M7" & is.finite(x$est) & is.finite(x$se), ]
  data.frame(s = sub("\\.rds$", "", f), n = nrow(x),
             est098 = mean(x$est <= log(0.98)),
             up112  = mean(x$est + 1.96 * x$se <= log(1.12)))
}))
write.csv(out, "results/11-13-26_posthoc-pooled-sibling-check.csv", row.names = FALSE)
