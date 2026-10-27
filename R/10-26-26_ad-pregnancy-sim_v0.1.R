# =============================================================================
# Can the standard designs see harm?
# Simulation of prenatal antidepressant exposure and child autism diagnosis
# Protocol: 10-26-26_Antidepressants-pregnancy_simulation-spec_ADEMP_v0.1.pdf
# Vera Wilde, TUM Campus Heilbronn. Skeleton v0.1, 26 Oct 2026.
#
# Structure
#   1. Scenario parameters
#   2. Data-generating mechanism (one row per pregnancy, two per family)
#   3. Calibration: tune confounding strength so the crude RR is ~1.7
#   4. Estimators M1-M8
#   5. One replicate -> all estimates
#   6. Performance measures (with Monte Carlo SEs)
#   7. Pilot run
#
# Decisions locked 26 Oct 2026: autism primary outcome; surveillance-bias
# factor included; binary outcome at age 10; target journal Epidemiology.
#
# Model notes
#   - Outcome risk is log-linear with no interactions, so the causal risk
#     ratio is collapsible: the target log RR equals log(rr_drug) for the
#     total effect. M3 (adjusting for preterm birth) therefore shows bias by
#     construction when preterm birth mediates part of the drug effect.
#   - Risks are rare (~1-2%), so the conditional odds ratio from the
#     sibling analysis (M7) approximates the risk ratio.
#   - Population-level methods use cluster-robust SEs by family.
# =============================================================================

suppressPackageStartupMessages(library(survival))

# -----------------------------------------------------------------------------
# 1. Scenario parameters
# -----------------------------------------------------------------------------
base_scenario <- list(
  rr_drug        = 1.0,   # true causal RR of maternal use on autism (1.0, 1.2, 1.5)
  rho_severity   = 0.6,   # between-pregnancy correlation of depression severity (0.9, 0.6, 0.3)
  capture        = 0.85,  # share of true use captured by dispensing records (1.0, 0.85, 0.70)
  half_sib       = 0.15,  # share of second pregnancies with a different father (0, 0.15, 0.30)
  rr_paternal    = 1.0,   # paternal preconception effect RR (1.0, 1.2)
  surveillance   = 1.0,   # diagnosis-rate multiplier for recorded-exposed children (1.0, 1.15)
  beta_specific  = 0.10,  # log-RR per SD of pregnancy-specific severity (non-shared confounding)
  rr_preterm_drug = 1.3,  # drug effect on preterm birth (makes preterm birth a mediator)
  rr_preterm_y   = 1.5,   # preterm birth effect on autism
  target_crude   = 1.7,   # calibration target for crude RR (Chan 2026 pooled crude ASD RR 1.69)
  target_paternal = 1.28, # calibration target for crude paternal RR (Chan 2026 paternal ASD RR 1.28)
  base_risk      = 0.015, # autism risk by age 10 in an unexposed child at average covariates
  n_fam          = 100000 # families per replicate (two pregnancies each); to be set so
                          # discordant-pair counts match the source studies
)

make_scenario <- function(...) modifyList(base_scenario, list(...))

# -----------------------------------------------------------------------------
# 2. Data-generating mechanism
#    conf_scale multiplies every confounding path into the outcome; it is
#    set by calibrate() so that crude RR hits target_crude in each scenario.
# -----------------------------------------------------------------------------
simulate_cohort <- function(sc, conf_scale = 1, seed = NULL, conf_pat = 1) {
  if (!is.null(seed)) set.seed(seed)
  n <- sc$n_fam
  fam <- rep(seq_len(n), each = 2)
  k   <- rep(1:2, times = n)

  # Family-level (shared) variables
  G_M <- rnorm(n)                       # maternal genetic liability
  S   <- rnorm(n)                       # socioeconomic / other confounders
  U   <- 0.5 * G_M + sqrt(0.75) * rnorm(n)  # shared component of severity

  # Fathers: second pregnancy has a different father with prob half_sib
  G_F1 <- rnorm(n)
  new_father <- rbinom(n, 1, sc$half_sib)
  G_F2 <- ifelse(new_father == 1, rnorm(n), G_F1)
  PD1  <- 0.5 * G_F1 + sqrt(0.75) * rnorm(n)        # paternal depression
  PD2  <- ifelse(new_father == 1, 0.5 * G_F2 + sqrt(0.75) * rnorm(n), PD1)

  # Expand to pregnancy level
  G_Mp <- G_M[fam]; Sp <- S[fam]; Up <- U[fam]
  G_F  <- ifelse(k == 1, G_F1[fam], G_F2[fam])
  PD   <- ifelse(k == 1, PD1[fam], PD2[fam])
  N    <- 2 * n

  # Severity in each pregnancy: shared + pregnancy-specific
  D <- sqrt(sc$rho_severity) * Up + sqrt(1 - sc$rho_severity) * rnorm(N)

  # Recorded prior depression history (noisy, family-level proxy of U)
  H <- as.integer((U + rnorm(n, sd = 0.8)) > 1.0)[fam]
  # HDPS-style continuous proxy of current severity (noisy)
  Dproxy <- D + rnorm(N, sd = 1.0)

  # Use in the 90 days before LMP (eligibility for the target trial)
  P0 <- rbinom(N, 1, plogis(-3.2 + 1.2 * D + 0.3 * Sp))
  # Continuation after recognition among prior users; rare new initiation
  A <- ifelse(P0 == 1,
              rbinom(N, 1, plogis(0.3 + 0.8 * D)),
              rbinom(N, 1, plogis(-6.0 + 1.2 * D)))

  # Dispensing records: capture true use with prob `capture`;
  # some discontinuers have a fill they do not take
  A_rec <- ifelse(A == 1, rbinom(N, 1, sc$capture),
                  ifelse(P0 == 1, rbinom(N, 1, 0.5 * (1 - sc$capture)), 0L))

  # Paternal use around conception
  F_use <- rbinom(N, 1, plogis(-3.0 + 0.9 * PD))
  F_rec <- F_use   # assume paternal dispensing recorded without error (v0.1)

  # Preterm birth: affected by drug and severity (mediator of drug effect)
  PT <- rbinom(N, 1, pmin(1, 0.07 * sc$rr_preterm_drug^A * exp(0.2 * D)))

  # Outcome: autism diagnosis by age 10
  # Direct drug effect is set so that the TOTAL effect (direct + via PT)
  # equals rr_drug on the risk-ratio scale (approximately; PT path is small).
  # Confounding paths into the outcome. Shared paths (family severity
  # component, both genetic liabilities, SES) are scaled by conf_scale and
  # calibrated. The pregnancy-specific severity path (beta_specific) is a
  # fixed scenario parameter: it is what sibling designs cannot remove.
  D_specific <- D - sqrt(sc$rho_severity) * Up
  lp <- log(sc$base_risk) +
        conf_scale * (0.35 * Up + 0.30 * G_Mp - 0.10 * Sp) +
        conf_pat * 0.30 * (G_F + PD) +
        sc$beta_specific * D_specific +
        log(sc$rr_paternal) * F_use +
        log(sc$rr_preterm_y) * (PT - mean(PT))
  rr_direct <- sc$rr_drug / ((1 + 0.07 * (sc$rr_preterm_drug - 1) *
                              (sc$rr_preterm_y - 1)) )
  lp <- lp + log(rr_direct) * A
  risk <- pmin(1, exp(lp) * ifelse(A_rec == 1, sc$surveillance, 1))
  Y <- rbinom(N, 1, risk)

  data.frame(fam, k, parity2 = as.integer(k == 2), S = Sp, H, Dproxy, D,
             P0, A, A_rec, F_rec, PT, Y, new_father = new_father[fam])
}

# -----------------------------------------------------------------------------
# 3. Calibration: bisection on conf_scale so crude RR (recorded exposure)
#    hits target_crude. Uses one large draw per scenario.
# -----------------------------------------------------------------------------
crude_rr <- function(d) mean(d$Y[d$A_rec == 1]) / mean(d$Y[d$A_rec == 0])

crude_rr_pat <- function(d) mean(d$Y[d$F_rec == 1]) / mean(d$Y[d$F_rec == 0])

bisect <- function(f, lo, hi, tol, label) {
  flo <- f(lo); fhi <- f(hi)
  if (flo > 0) return(list(value = lo, ok = FALSE,
                           note = paste(label, "target exceeded with this path switched off")))
  if (fhi < 0) return(list(value = hi, ok = FALSE,
                           note = paste(label, "target not reached at upper bound")))
  for (i in 1:25) {
    mid <- (lo + hi) / 2; fm <- f(mid)
    if (abs(fm) < tol) break
    if (fm > 0) hi <- mid else lo <- mid
  }
  list(value = mid, ok = TRUE, note = "")
}

# Two-step calibration to two observed quantities:
#   step 1: paternal-side confounding (conf_pat) so crude paternal RR = target_paternal
#   step 2: maternal-side confounding (conf_scale) so crude maternal RR = target_crude
# A scenario where step 2 fails because the crude maternal RR already exceeds
# the target with maternal confounding switched off is INCOMPATIBLE with the
# observed data: that true effect size, plus the familial confounding implied
# by the paternal association, produces a crude RR larger than observed.
calibrate <- function(sc, tol = 0.02, seed = 1) {
  sc_big <- modifyList(sc, list(n_fam = max(sc$n_fam, 400000)))
  p <- bisect(function(cp) crude_rr_pat(simulate_cohort(sc_big, 0, seed, cp)) - sc$target_paternal,
              0, 4, tol, "paternal")
  m <- bisect(function(cs) crude_rr(simulate_cohort(sc_big, cs, seed, p$value)) - sc$target_crude,
              0, 4, tol, "maternal")
  list(conf_scale = m$value, conf_pat = p$value,
       compatible = p$ok && m$ok, note = trimws(paste(p$note, m$note)))
}

# -----------------------------------------------------------------------------
# 4. Estimators. Each returns c(est = log RR, se = SE of log RR)
# -----------------------------------------------------------------------------
# Poisson regression with cluster-robust (family) sandwich SE
pois_rr <- function(formula, d, term, weights = NULL) {
  d$.w <- if (is.null(weights)) 1 else weights
  fit <- suppressWarnings(glm(formula, family = poisson(), data = d, weights = .w))
  X  <- model.matrix(fit); mu <- fitted(fit); y <- fit$y; w <- d$.w
  bread <- solve(crossprod(X, X * (w * mu)))
  sc <- X * (w * (y - mu))
  sc_cl <- rowsum(sc, d$fam)
  V <- bread %*% crossprod(sc_cl) %*% bread
  j <- which(colnames(X) == term)
  c(est = unname(coef(fit)[term]), se = sqrt(V[j, j]))
}

m1_crude  <- function(d) pois_rr(Y ~ A_rec, d, "A_rec")
m2_adj    <- function(d) pois_rr(Y ~ A_rec + S + H + parity2, d, "A_rec")
m3_badctl <- function(d) pois_rr(Y ~ A_rec + S + H + parity2 + PT, d, "A_rec")
m4_iptw   <- function(d) {
  ps <- fitted(glm(A_rec ~ S + H + parity2 + Dproxy, family = binomial(), data = d))
  p  <- mean(d$A_rec)
  w  <- ifelse(d$A_rec == 1, p / ps, (1 - p) / (1 - ps))   # stabilized ATE weights
  w  <- pmin(w, quantile(w, 0.999))                          # trim extremes
  pois_rr(Y ~ A_rec, d, "A_rec", weights = w)
}
m5_restrict <- function(d) pois_rr(Y ~ A_rec + S + parity2, d[d$H == 1, ], "A_rec")
m6_discont  <- function(d) pois_rr(Y ~ A_rec + S + H + parity2, d[d$P0 == 1, ], "A_rec")
m7_sibling  <- function(d) {
  disc <- ave(d$A_rec, d$fam, FUN = function(a) length(unique(a)) > 1) == 1
  ds <- d[disc, ]
  if (nrow(ds) < 20 || sum(ds$Y) < 5) return(c(est = NA, se = NA))
  fit <- tryCatch(clogit(Y ~ A_rec + parity2 + strata(fam), data = ds,
                         method = "approximate"), error = function(e) NULL)
  if (is.null(fit)) return(c(est = NA, se = NA))
  c(est = unname(coef(fit)["A_rec"]), se = sqrt(vcov(fit)["A_rec", "A_rec"]))
}
m8_paternal <- function(d) pois_rr(Y ~ F_rec + S + parity2, d, "F_rec")

methods <- list(M1 = m1_crude, M2 = m2_adj, M3 = m3_badctl, M4 = m4_iptw,
                M5 = m5_restrict, M6 = m6_discont, M7 = m7_sibling, M8 = m8_paternal)

# -----------------------------------------------------------------------------
# 5. One replicate
# -----------------------------------------------------------------------------
one_rep <- function(sc, cal, seed) {
  d <- simulate_cohort(sc, cal$conf_scale, seed = seed, conf_pat = cal$conf_pat)
  out <- t(sapply(methods, function(m) m(d)))
  data.frame(method = rownames(out), est = out[, "est"], se = out[, "se"],
             n_disc_pairs = sum(tapply(d$A_rec, d$fam, function(a) length(unique(a)) > 1)),
             crude_rr = crude_rr(d), crude_rr_pat = crude_rr_pat(d),
             prev_exposed = mean(d$A_rec),
             risk_unexposed = mean(d$Y[d$A_rec == 0]), row.names = NULL)
}

# -----------------------------------------------------------------------------
# 6. Performance measures
#    truth: log(rr_drug) for M1-M7; for M8 the "truth" for a negative control
#    is log(rr_paternal) (zero when the paternal pathway is switched off).
# -----------------------------------------------------------------------------
performance <- function(res, sc) {
  res$truth <- ifelse(res$method == "M8", log(sc$rr_paternal), log(sc$rr_drug))
  res$lo <- res$est - 1.96 * res$se; res$hi <- res$est + 1.96 * res$se
  res$covers     <- res$lo <= res$truth & res$hi >= res$truth
  res$reassuring <- res$lo <= 0 & res$hi >= 0     # interval includes RR = 1
  agg <- do.call(rbind, lapply(split(res, res$method), function(r) {
    r <- r[!is.na(r$est), ]; n <- nrow(r)
    b <- mean(r$est - r$truth); esd <- sd(r$est)
    cov <- mean(r$covers); rea <- mean(r$reassuring)
    data.frame(method = r$method[1], n_rep = n,
               bias_logRR = b, mcse_bias = esd / sqrt(n),
               emp_se = esd, mean_est_RR = exp(mean(r$est)),
               coverage = cov, mcse_cov = sqrt(cov * (1 - cov) / n),
               p_reassuring = rea, mcse_reas = sqrt(rea * (1 - rea) / n),
               mean_ci_width = mean(r$hi - r$lo))
  }))
  # Published evidence pattern: crude above 1, sibling includes 1, paternal above 1
  w <- reshape(res[, c("rep", "method", "lo", "hi")], idvar = "rep",
               timevar = "method", direction = "wide")
  pattern <- with(w, lo.M1 > 0 & lo.M7 <= 0 & hi.M7 >= 0 & lo.M8 > 0)
  list(by_method = agg, p_pattern = mean(pattern, na.rm = TRUE),
       n_pattern = sum(!is.na(pattern)))
}

run_scenario <- function(sc, n_rep, seed0 = 1000) {
  cal <- calibrate(sc)
  res <- do.call(rbind, lapply(seq_len(n_rep), function(r) {
    x <- one_rep(sc, cal, seed = seed0 + r); x$rep <- r; x
  }))
  list(scenario = sc, cal = cal, raw = res, perf = performance(res, sc))
}

# -----------------------------------------------------------------------------
# 7. Pilot: checks code and calibration only (not results; see protocol s.9)
# -----------------------------------------------------------------------------
if (sys.nframe() == 0) {
  n_rep_pilot <- as.integer(Sys.getenv("N_REP", "10"))
  pilot <- list(
    null   = make_scenario(rr_drug = 1.0),
    harm12 = make_scenario(rr_drug = 1.2),
    harm15 = make_scenario(rr_drug = 1.5)
  )
  results <- lapply(pilot, run_scenario, n_rep = n_rep_pilot)
  for (nm in names(results)) {
    r <- results[[nm]]
    cat("\n=== Pilot scenario:", nm, " (true RR =", r$scenario$rr_drug, ") ===\n")
    cat(sprintf("compatible with observed data: %s %s\n", r$cal$compatible, r$cal$note))
    cat(sprintf("conf_scale = %.3f; conf_pat = %.3f; crude RR = %.2f; crude paternal RR = %.2f\n",
                r$cal$conf_scale, r$cal$conf_pat, mean(r$raw$crude_rr), mean(r$raw$crude_rr_pat)))
    cat(sprintf("exposed = %.1f%%; unexposed risk = %.2f%%; discordant pairs = %.0f\n",
                100 * mean(r$raw$prev_exposed), 100 * mean(r$raw$risk_unexposed),
                mean(r$raw$n_disc_pairs)))
    print(format(r$perf$by_method, digits = 3), row.names = FALSE)
    cat(sprintf("P(published pattern) = %.2f over %d reps\n",
                r$perf$p_pattern, r$perf$n_pattern))
  }
  p_0 <- results$null$perf$p_pattern
  for (h in c("harm12", "harm15")) {
    p_h <- results[[h]]$perf$p_pattern
    cat(sprintf("Pattern ratio, P(pattern | %s) / P(pattern | null) = %s\n", h,
                ifelse(p_0 > 0, sprintf("%.2f", p_h / p_0), "undefined (null count 0)")))
  }
  saveRDS(results, "pilot_results.rds")
}
