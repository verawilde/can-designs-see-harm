# Summarize the scenario grid: per scenario x method, mean estimated RR, coverage,
# P(reassuring reading), and the published-pattern likelihood ratio vs the null.
files <- list.files("grid", pattern = "\\.rds$", full.names = TRUE)
res <- lapply(files, readRDS)

rows <- do.call(rbind, lapply(res, function(r) {
  b <- r$perf$by_method
  setting <- sub("__rr.*", "", r$key)
  data.frame(setting = setting, rr = r$scenario$rr_drug, method = b$method,
             mean_RR = b$mean_est_RR, coverage = b$coverage, p_reassuring = b$p_reassuring,
             mcse_reas = b$mcse_reas, n_rep = b$n_rep, compatible = r$cal$compatible,
             p_pattern = r$perf$p_pattern,
             disc_pairs = mean(r$raw$n_disc_pairs), crude = mean(r$raw$crude_rr),
             crude_pat = mean(r$raw$crude_rr_pat))
}))
write.csv(rows, "10-27-26_grid-results_v0.1.csv", row.names = FALSE)

# Key table: sibling (M7) and HDPS-like (M4) reassurance probability under true harm
key <- subset(rows, method %in% c("M2", "M4", "M7", "M8"))
wide <- reshape(key[, c("setting", "rr", "method", "p_reassuring")],
                idvar = c("setting", "rr"), timevar = "method", direction = "wide")
wide <- wide[order(wide$setting, wide$rr), ]
print(format(wide, digits = 2), row.names = FALSE)

# Pattern likelihood ratio per setting
pat <- unique(rows[, c("setting", "rr", "p_pattern", "compatible", "disc_pairs", "crude", "crude_pat")])
lr <- do.call(rbind, lapply(split(pat, pat$setting), function(p) {
  p0 <- p$p_pattern[p$rr == 1.0]
  data.frame(setting = p$setting[1],
             p_null = p0, p_rr1.2 = p$p_pattern[p$rr == 1.2], p_rr1.5 = p$p_pattern[p$rr == 1.5],
             compatible_all = all(p$compatible),
             disc_pairs = round(mean(p$disc_pairs)))
}))
print(format(lr, digits = 2), row.names = FALSE)
write.csv(lr, "10-27-26_grid-pattern_v0.1.csv", row.names = FALSE)

# Figure: P(sibling interval includes 1) when the true RR is 1.2 and 1.5, by setting. Black and white.
pdf("10-27-26_grid-figure_v0.1.pdf", width = 9, height = 6.5, pointsize = 14)
s7 <- subset(rows, method == "M7" & rr > 1)
lev <- unique(s7$setting[order(s7$p_reassuring[s7$rr == 1.2])])
s7$setting <- factor(s7$setting, levels = lev)
par(mar = c(4.5, 10, 3, 1))
y12 <- as.numeric(s7$setting[s7$rr == 1.2]); y15 <- as.numeric(s7$setting[s7$rr == 1.5])
plot(NA, xlim = c(0, 1), ylim = c(0.5, length(lev) + 0.5), yaxt = "n",
     xlab = "Share of replicates where the sibling interval includes RR = 1",
     ylab = "", main = "Sibling comparison under true harm")
axis(2, at = seq_along(lev), labels = lev, las = 1, cex.axis = 0.85)
abline(h = seq_along(lev), col = "grey85", lty = 3)
points(s7$p_reassuring[s7$rr == 1.2], y12, pch = 1, cex = 1.4)
points(s7$p_reassuring[s7$rr == 1.5], y15, pch = 19, cex = 1.2)
legend("bottomright", c("true RR 1.2", "true RR 1.5"), pch = c(1, 19), bty = "n")
dev.off()
cat("summary written\n")
