# Exploratory (not preregistered): random-effects meta-analysis of simulated sibling studies,
# built from the stored M7 replicates of the final run.
suppressMessages(library(stats))
# Run from the repository root: Rscript R/explore/11-13-26_simulated-meta-analysis.R
res <- read.csv("results/11-12-26_final-results_v1.0.csv")
m7  <- res[res$method == "M7", ]
dir <- "results/final-replicates"
load_m7 <- function(key) { x <- readRDS(file.path(dir, paste0(key, ".rds")))$raw
  x <- x[x$method == "M7" & is.finite(x$est) & is.finite(x$se) & x$se < 2 & abs(x$est) < 5, c("est","se")]; x }
dl <- function(y, s) { w <- 1/s^2; mu <- sum(w*y)/sum(w); Q <- sum(w*(y-mu)^2); k <- length(y)
  tau2 <- max(0, (Q-(k-1))/(sum(w)-sum(w^2)/sum(w))); ws <- 1/(s^2+tau2)
  m <- sum(ws*y)/sum(ws); se <- sqrt(1/sum(ws)); c(m, se, tau2) }
set.seed(20261113)
NMA <- 20000; K <- 7
out <- list()
for (rr in c("1.0","1.2","1.5")) {
  keys <- m7$setting[grepl(paste0("__rr", rr, "$"), m7$setting) & m7$compatible %in% c(TRUE,"TRUE")]
  large <- setdiff(keys, paste0("brown_size__rr", rr))
  truth <- log(as.numeric(rr))
  mm <- m7[match(large, m7$setting), ]
  nullward <- mm$setting[order(mm$mean_RR)][1:5]
  data <- lapply(setNames(c(large, paste0("brown_size__rr", rr)), c(large, paste0("brown_size__rr", rr))), load_m7)
  designs <- list(
    base_only      = list(pool = paste0("base__rr", rr), nlarge = 2),
    mixed          = list(pool = large, nlarge = 2),
    nullward_only  = list(pool = nullward, nlarge = 2))
  for (d in names(designs)) {
    pool <- designs[[d]]$pool; nl <- designs[[d]]$nlarge
    r <- t(replicate(NMA, {
      sc <- c(sample(pool, nl, replace = TRUE),
              rep(paste0("brown_size__rr", rr), K - nl))
      ys <- ss <- numeric(K)
      for (j in seq_len(K)) {
        src <- if (grepl("brown_size", sc[j])) {
          # small study: draw a scenario's bias from the pool, add Brown-size sampling error
          b <- data[[sample(pool, 1)]]; bb <- data[[sc[j]]]
          i <- sample(nrow(bb), 1); c(mean(b$est) + (bb$est[i] - mean(bb$est)), bb$se[i])
        } else { x <- data[[sc[j]]]; i <- sample(nrow(x), 1); c(x$est[i], x$se[i]) }
        ys[j] <- src[1]; ss[j] <- src[2] }
      dl(ys, ss) }))
    lo <- r[,1] - 1.96*r[,2]; hi <- r[,1] + 1.96*r[,2]
    out[[length(out)+1]] <- data.frame(truth = rr, design = d, n_scen = length(pool),
      median_pooled_se = median(r[,2]), mean_RR = exp(mean(r[,1])),
      p_est_le_098 = mean(r[,1] <= log(0.98)), p_up_le_112 = mean(hi <= log(1.12)),
      p_incl1 = mean(lo <= 0 & hi >= 0))
  }
}
o <- do.call(rbind, out); print(o, digits = 3)
write.csv(o, "results/explore/11-13-26_simulated-meta-analysis_results.csv", row.names = FALSE)

