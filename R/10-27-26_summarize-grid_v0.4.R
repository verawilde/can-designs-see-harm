# Summaries for a grid directory (default: v0.4 pilot). Adds the preregistered measures:
# success share, interval entirely above 1, compatible with both 1 and 1.5, excludes 1.2,
# excludes 1.5, published-pattern probability, calibrated confounding parameters.
args <- commandArgs(TRUE)
dir <- if (length(args)) args[1] else "grid_v04"
out <- if (length(args) > 1) args[2] else "10-27-26_grid-results_v0.4.csv"
files <- list.files(dir, pattern = "\\.rds$", full.names = TRUE)
rows <- do.call(rbind, lapply(files, function(f) {
  r <- readRDS(f); b <- r$perf$by_method; raw <- r$raw
  n_tot <- max(raw$rep)
  raw$lo <- raw$est - 1.96 * raw$se; raw$hi <- raw$est + 1.96 * raw$se
  ok <- is.finite(raw$est) & is.finite(raw$se)
  m <- function(x) tapply(x[ok], raw$method[ok], mean)
  data.frame(setting = r$key, method = b$method, truth_RR = b$truth_RR, mean_RR = b$mean_est_RR,
    coverage = b$coverage, p_includes1 = b$p_reassuring,
    p_above1 = as.numeric(m(raw$lo > 0)[b$method]),
    p_both_1_and_1.5 = as.numeric(m(raw$lo <= 0 & raw$hi >= log(1.5))[b$method]),
    p_excl_1.2 = as.numeric(m(raw$hi < log(1.2) | raw$lo > log(1.2))[b$method]),
    p_excl_1.5 = as.numeric(m(raw$hi < log(1.5) | raw$lo > log(1.5))[b$method]),
    p_success = as.numeric(tapply(ok, raw$method, mean)[b$method]),
    p_pattern = r$perf$p_pattern, compatible = r$cal$compatible,
    conf_scale = r$cal$conf_scale, conf_pat = r$cal$conf_pat,
    crude = mean(raw$crude_rr[raw$method == "M1"]), n_rep = n_tot)
}))
write.csv(rows, out, row.names = FALSE)
keep <- c("M1","M2","M4","M6","M7","M7c","M8","M8b","M7noid","M7id","M2high","M7high")
for (s in unique(rows$setting)) {
  x <- rows[rows$setting == s & rows$method %in% keep, ]
  cat(sprintf("\n== %s  compatible=%s crude=%.2f pattern=%.2f conf_scale=%.2f conf_pat=%.2f\n", s,
      x$compatible[1], x$crude[1], x$p_pattern[1], x$conf_scale[1], x$conf_pat[1]))
  print(format(x[, c("method","truth_RR","mean_RR","coverage","p_includes1","p_above1","p_both_1_and_1.5","p_success")], digits = 2), row.names = FALSE)
}
