# ---------------------------------------------------------------------------
# 06_inference.R
#
# Formal tests and uncertainty intervals added in response to supervisor
# comments on the 2026-08-26 draft:
#
#   C4   Do the observed rates STATISTICALLY exceed the corrected null?
#        -> Monte-Carlo exceedance tests against the corrected null (global
#           and study-specific), plus a study-bootstrap exceedance test that
#           respects the clustering of comparisons within studies.
#   C6   Error bars for every bar in Figs 1-2, and the material for the
#        slope-with-band overlay in Fig. 1a.
#   C10  Supporting statistics for the secondary-archive paragraph.
#   C16  One reporting template used for all three analyses (GPT-4 x crowd,
#        GPT-4 x other LLMs, GPT-4 x experts): joint-miss rate, share of
#        joint misses in the same direction, corrected null, exceedance P.
#   C17  The error-on-error slope reported for BOTH archives, with
#        cluster-robust confidence intervals.
#
# INPUT   data/derived/comparisons_primary.csv
#         data/derived/comparisons_all_models.csv
#         data/derived/comparisons_archive2.csv
# OUTPUT  output/tables/exceedance_tests.csv
#         output/tables/consistent_reporting.csv
#         data/derived/inference_objects.rds   (consumed by 05_figures.R)
#
# Run AFTER 03_main_results.R and BEFORE 05_figures.R (see run_all.R).
#
# Wilhelmsen, Esfandiari & Gollwitzer
# ---------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(purrr)
library(here)

source(here("R", "00_functions.R"))

derived_dir <- here("data", "derived")
table_dir   <- here("output", "tables")
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

wide <- read.csv(file.path(derived_dir, "comparisons_primary.csv"), check.names = FALSE)

miss_g <- abs(wide$zg) > 1.96
miss_h <- abs(wide$zh) > 1.96
joint  <- miss_g & miss_h
same   <- joint & (sign(wide$zg) == sign(wide$zh))
n_cmp  <- nrow(wide)

obs <- c(joint     = mean(joint),
         same_side = mean(same),
         direction = sum(same) / sum(joint))

nulls <- null_exact(mean(miss_g), mean(miss_h))


# --- 1. Exceedance test against the corrected null (C4) ---------------------
# The corrected null gives the EXPECTED rates for independent predictors
# scored against the same noisy estimate; it does not by itself say whether
# the observed rates sit outside what sampling variation under that null
# could produce. Simulating full archives of n = 1,678 comparisons from the
# null supplies that finite-sample distribution, and the one-sided P is the
# share of null archives whose rate reaches the observed one.
#
# Under this null, comparisons are independent, whereas real comparisons
# cluster within studies; the companion tests in sections 2 and 3 address
# that. All three lead to the same conclusion.

set.seed(41)
R_NULL <- 10000

sd_for <- function(target) {
  # total spread s of z = s*b - u reproducing a marginal miss rate, matching
  # the calibration used in 05_figures.R and 04_robustness.R
  uniroot(function(s) {
    2 * (1 - pnorm(1.96 / sqrt(s^2 + 1))) - target
  }, c(0.01, 40))$root
}
sg <- sd_for(mean(miss_g))
sh <- sd_for(mean(miss_h))

null_draw <- function(n) {
  u  <- rnorm(n)
  zg <- sg * rnorm(n) - u
  zh <- sh * rnorm(n) - u
  j  <- abs(zg) > 1.96 & abs(zh) > 1.96
  s  <- j & sign(zg) == sign(zh)
  c(joint = mean(j), same_side = mean(s),
    direction = if (sum(j)) sum(s) / sum(j) else NA_real_)
}

null_dist <- t(replicate(R_NULL, null_draw(n_cmp)))

p_mc <- vapply(names(obs), function(k) {
  (1 + sum(null_dist[, k] >= obs[[k]], na.rm = TRUE)) / (1 + R_NULL)
}, numeric(1))

null_ci <- apply(null_dist, 2, quantile, c(.025, .975), na.rm = TRUE)

message("Global corrected null, finite-sample Monte Carlo (", R_NULL, " archives):")
for (k in names(obs)) {
  message(sprintf("  %-9s observed %5.1f%%  null mean %5.1f%% [%4.1f, %4.1f]  P %s",
                  k, 100 * obs[[k]], 100 * mean(null_dist[, k], na.rm = TRUE),
                  100 * null_ci[1, k], 100 * null_ci[2, k],
                  format.pval(p_mc[[k]], digits = 2)))
}

report_value("Main text, para 9 (C4)", "P, same-side rate vs corrected null", p_mc[["same_side"]], 6)
report_value("Main text, para 9 (C4)", "P, direction share vs corrected null", p_mc[["direction"]], 6)
report_value("Main text, para 9 (C4)", "P, joint rate vs corrected null", p_mc[["joint"]], 6)


# --- 2. Exceedance test against the study-specific null (C4) ----------------
# The stricter benchmark of SI 1.5 lets every study set its own difficulty.
# Re-run per replicate archive so that a null DISTRIBUTION, not just a null
# mean, is available. Spreads are calibrated once per study (as in
# 04_robustness.R) and reused across replicates.

set.seed(42)
R_SS <- 2000

POOL   <- 1e6
pool_b <- rnorm(POOL); pool_u <- rnorm(POOL)
rate_for <- function(s) mean(abs(s * pool_b - pool_u) > 1.96)
calib_s  <- function(target) {
  if (target <= rate_for(0))  return(0)
  if (target >= rate_for(40)) return(40)
  uniroot(function(s) rate_for(s) - target, c(0, 40), tol = 1e-4)$root
}

study_tab <- wide %>%
  group_by(study) %>%
  summarise(n = n(), mg = mean(abs(zg) > 1.96), mh = mean(abs(zh) > 1.96),
            .groups = "drop") %>%
  mutate(s1 = vapply(mg, calib_s, numeric(1)),
         s2 = vapply(mh, calib_s, numeric(1)))

ss_draw <- function() {
  tot_n <- sum(study_tab$n); jj <- 0L; ss <- 0L
  for (i in seq_len(nrow(study_tab))) {
    nk <- study_tab$n[i]
    u  <- rnorm(nk)
    z1 <- study_tab$s1[i] * rnorm(nk) - u
    z2 <- study_tab$s2[i] * rnorm(nk) - u
    j  <- abs(z1) > 1.96 & abs(z2) > 1.96
    jj <- jj + sum(j)
    ss <- ss + sum(j & sign(z1) == sign(z2))
  }
  c(joint = jj / tot_n, same_side = ss / tot_n,
    direction = if (jj) ss / jj else NA_real_)
}

ss_dist <- t(replicate(R_SS, ss_draw()))

p_ss <- vapply(names(obs), function(k) {
  (1 + sum(ss_dist[, k] >= obs[[k]], na.rm = TRUE)) / (1 + R_SS)
}, numeric(1))

ss_ci <- apply(ss_dist, 2, quantile, c(.025, .975), na.rm = TRUE)

message("Study-specific null, finite-sample Monte Carlo (", R_SS, " archives):")
for (k in names(obs)) {
  message(sprintf("  %-9s observed %5.1f%%  null mean %5.1f%% [%4.1f, %4.1f]  P %s",
                  k, 100 * obs[[k]], 100 * mean(ss_dist[, k], na.rm = TRUE),
                  100 * ss_ci[1, k], 100 * ss_ci[2, k],
                  format.pval(p_ss[[k]], digits = 2)))
}

report_value("SI 1.5 (C4)", "P, same-side rate vs study-specific null", p_ss[["same_side"]], 6)
report_value("SI 1.5 (C4)", "P, direction share vs study-specific null", p_ss[["direction"]], 6)
report_value("SI 1.5 (C4)", "P, joint rate vs study-specific null", p_ss[["joint"]], 6)


# --- 3. Bootstrap exceedance: clustering respected (C4) ---------------------
# The Monte-Carlo tests above draw comparisons independently. As a companion
# that keeps the real dependence structure, resample STUDIES with replacement
# (as in 03_main_results.R) and ask how often the resampled rate falls at or
# below the corrected-null benchmark. This inverts the study-bootstrap
# confidence interval into a test.

set.seed(43)
B <- 4000
studies <- unique(wide$study)
idx_by_study <- split(seq_len(nrow(wide)), wide$study)

boot <- replicate(B, {
  s  <- sample(studies, replace = TRUE)
  ii <- unlist(idx_by_study[s], use.names = FALSE)
  j  <- joint[ii]; sm <- same[ii]
  c(joint = mean(j), same_side = mean(sm),
    direction = if (sum(j)) sum(sm) / sum(j) else NA_real_)
})

p_boot <- c(
  joint     = mean(boot["joint", ]     <= nulls[["joint"]]),
  same_side = mean(boot["same_side", ] <= nulls[["same"]]),
  direction = mean(boot["direction", ] <= nulls[["direction"]], na.rm = TRUE)
)

boot_ci <- apply(boot, 1, quantile, c(.025, .975), na.rm = TRUE)

message("Study-bootstrap exceedance (share of ", B, " draws at or below the null):")
for (k in names(p_boot)) {
  message(sprintf("  %-9s  P %s", k, format.pval(p_boot[[k]], digits = 2)))
}

report_value("Main text, para 9 (C4)", "bootstrap P, same-side rate", p_boot[["same_side"]], 6)
report_value("Main text, para 9 (C4)", "bootstrap P, direction share", p_boot[["direction"]], 6)

exceedance <- data.frame(
  metric        = names(obs),
  observed      = as.numeric(obs),
  null_global   = as.numeric(nulls[c("joint", "same", "direction")]),
  null_global_lo = null_ci[1, ],
  null_global_hi = null_ci[2, ],
  p_global      = as.numeric(p_mc),
  null_ss_lo    = ss_ci[1, ],
  null_ss_hi    = ss_ci[2, ],
  p_ss          = as.numeric(p_ss),
  p_bootstrap   = as.numeric(p_boot[names(obs)])
)
write.csv(exceedance, file.path(table_dir, "exceedance_tests.csv"), row.names = FALSE)
print(exceedance %>% mutate(across(where(is.numeric), ~ signif(., 3))))


# --- 4. Slopes in BOTH archives, same convention (C17) ----------------------
# Error-on-error slope: the other predictor's standardized error regressed on
# GPT-4's, after recalibration, with the confidence interval clustered on the
# study (primary archive) or megastudy family (secondary archive).

sb_primary <- slope_band(wide$zg, wide$zh, wide$study,
                         grid = seq(-12, 12, length.out = 241))

message(sprintf("Primary archive slope: %.2f (95%% CI %.2f-%.2f, clustered by study)",
                sb_primary$slope,
                sb_primary$slope - 1.96 * sb_primary$slope_se,
                sb_primary$slope + 1.96 * sb_primary$slope_se))

report_value("Main text, para 8 (C17)", "crowd-on-GPT-4 error slope", sb_primary$slope, 2)
report_value("Main text, para 8 (C17)", "slope CI lower (clustered)",
             sb_primary$slope - 1.96 * sb_primary$slope_se, 2)
report_value("Main text, para 8 (C17)", "slope CI upper (clustered)",
             sb_primary$slope + 1.96 * sb_primary$slope_se, 2)

a2_path <- file.path(derived_dir, "comparisons_archive2.csv")
sb_expert <- NULL
if (file.exists(a2_path)) {
  a2 <- read.csv(a2_path, check.names = FALSE)
  a2 <- a2[!is.na(a2$z_recal_gpt) & !is.na(a2$z_recal_expert), ]

  cl2 <- if ("family" %in% names(a2)) a2$family else seq_len(nrow(a2))
  sb_expert <- slope_band(a2$z_recal_gpt, a2$z_recal_expert, cl2)

  message(sprintf(
    "Secondary archive slope: %.2f (95%% CI %.2f-%.2f, clustered by megastudy)",
    sb_expert$slope,
    sb_expert$slope - 1.96 * sb_expert$slope_se,
    sb_expert$slope + 1.96 * sb_expert$slope_se))

  report_value("Main text, para 12 (C17)", "expert-on-GPT-4 error slope", sb_expert$slope, 2)
  report_value("Main text, para 12 (C17)", "expert slope CI lower (clustered)",
               sb_expert$slope - 1.96 * sb_expert$slope_se, 2)
  report_value("Main text, para 12 (C17)", "expert slope CI upper (clustered)",
               sb_expert$slope + 1.96 * sb_expert$slope_se, 2)
}


# --- 5. One reporting template for all three analyses (C16, C7, C8) ---------
# For every predictor paired with GPT-4, on the comparisons where both are
# available: n, both predictors' miss rates, share of GPT-4's misses the
# other predictor repeats in the same direction (with a Wilson interval), the
# pair's corrected-null benchmark for that share, and the binomial P against
# it. The 28% the draft quotes without a source (C8) is this benchmark for
# the human crowd; the 60-78% (C7) are these shares for the other models.

zwide <- read.csv(file.path(derived_dir, "comparisons_all_models.csv"),
                  check.names = FALSE) %>%
  select(study, outcome.name, reference_condition, condition.name, model, z_recal) %>%
  pivot_wider(names_from = model, values_from = z_recal)

gp <- zwide$`gpt-4`

others <- c("deepseek/deepseek-chat-v3-0324", "openai/gpt-oss-120b", "gpt-3.5-turbo",
            "google/gemma-3-27b-it", "babbage-002", "davinci-002", "human")

template <- map_dfr(others, function(m) {
  zo    <- zwide[[m]]
  ok    <- !is.na(zo) & !is.na(gp)
  gmiss <- ok & abs(gp) > 1.96
  k     <- sum((abs(zo) > 1.96 & sign(zo) == sign(gp))[gmiss])
  n     <- sum(gmiss)
  ci    <- wilson_ci(k, n)
  gmr   <- mean(abs(gp[ok]) > 1.96)
  ne    <- null_exact(gmr, mean(abs(zo[ok]) > 1.96))
  bench <- ne[["same"]] / gmr                       # same-direction share | GPT-4 miss
  data.frame(
    analysis        = "primary: vs GPT-4 misses",
    predictor       = m,
    n_overlap       = sum(ok),
    n_gpt4_misses   = n,
    repeats_same_dir = k / n,
    wilson_lo       = ci[["lo"]],
    wilson_hi       = ci[["hi"]],
    corrected_null  = bench,
    p_binomial      = stats::pbinom(k - 1, n, bench, lower.tail = FALSE)
  )
})

# The expert row of the same template, from the secondary archive.
if (file.exists(a2_path)) {
  gm2   <- abs(a2$z_recal_gpt)    > 1.96
  em2   <- abs(a2$z_recal_expert) > 1.96
  same_dir2 <- gm2 & em2 & (sign(a2$z_recal_gpt) == sign(a2$z_recal_expert))
  k2 <- sum(same_dir2); n2g <- sum(gm2)
  ci2 <- wilson_ci(k2, n2g)
  ne2 <- null_exact(mean(gm2), mean(em2))
  bench2 <- ne2[["same"]] / mean(gm2)
  template <- bind_rows(template, data.frame(
    analysis        = "secondary: vs GPT-4 misses",
    predictor       = "domain experts",
    n_overlap       = nrow(a2),
    n_gpt4_misses   = n2g,
    repeats_same_dir = k2 / n2g,
    wilson_lo       = ci2[["lo"]],
    wilson_hi       = ci2[["hi"]],
    corrected_null  = bench2,
    p_binomial      = stats::pbinom(k2 - 1, n2g, bench2, lower.tail = FALSE)
  ))

  # The nested percentages the draft reports as "82%, and 99% of those"
  # (C16), each with its own denominator and interval.
  both2 <- gm2 & em2
  ci_also  <- wilson_ci(sum(both2), sum(gm2))          # experts miss | GPT-4 miss
  ci_dir2  <- wilson_ci(sum(same_dir2), sum(both2))    # same side | both miss

  report_value("Main text, para 12 (C16)", "experts also miss | GPT-4 miss", sum(both2) / sum(gm2))
  report_value("Main text, para 12 (C16)", "  Wilson CI lower", ci_also[["lo"]])
  report_value("Main text, para 12 (C16)", "  Wilson CI upper", ci_also[["hi"]])
  report_value("Main text, para 12 (C16)", "same side | both miss", sum(same_dir2) / sum(both2))
  report_value("Main text, para 12 (C16)", "  Wilson CI lower ", ci_dir2[["lo"]])
  report_value("Main text, para 12 (C16)", "  Wilson CI upper ", ci_dir2[["hi"]])
  report_value("Main text, para 12 (C16)", "null: same-direction share | GPT-4 miss", bench2)
  report_value("Main text, para 12 (C16)", "P (binomial) vs that null",
               stats::pbinom(k2 - 1, n2g, bench2, lower.tail = FALSE), 6)
}

template <- template %>%
  mutate(across(c(repeats_same_dir, wilson_lo, wilson_hi,
                  corrected_null, p_binomial), ~ signif(., 3)))
write.csv(template, file.path(table_dir, "consistent_reporting.csv"), row.names = FALSE)
print(template)

for (i in seq_len(nrow(template))) {
  report_value("Fig. 2 / SI 2.5 (C16, C7, C8)",
               paste0(template$predictor[i], ": repeats same direction | GPT-4 miss"),
               template$repeats_same_dir[i], 2)
  report_value("Fig. 2 / SI 2.5 (C16, C7, C8)",
               paste0(template$predictor[i], ": corrected null"),
               template$corrected_null[i], 2)
}


# --- 6. Intervals for the concentration analysis (C28, C10) -----------------
# "Nearly three times as common among the larger effects" and the screening
# cell both compare two conditional proportions; give each side its interval
# and a two-proportion test so the sentence can carry its own statistics.

big  <- abs(wide$estimate) >= 0.05
mag  <- abs((wide$pg + wide$ph) / 2)
gap  <- abs(wide$pg - wide$ph)
cell <- mag >= quantile(mag, .8) & sign(wide$pg) == sign(wide$ph) & gap <= median(gap)

split_stats <- function(flag, label) {
  k1 <- sum(same[flag]);  n1 <- sum(flag)
  k0 <- sum(same[!flag]); n0 <- sum(!flag)
  ci1 <- wilson_ci(k1, n1); ci0 <- wilson_ci(k0, n0)
  pt  <- stats::prop.test(c(k1, k0), c(n1, n0))
  message(sprintf("%s: %.1f%% [%.1f, %.1f] vs %.1f%% [%.1f, %.1f], chi-sq P %s",
                  label, 100 * k1 / n1, 100 * ci1[["lo"]], 100 * ci1[["hi"]],
                  100 * k0 / n0, 100 * ci0[["lo"]], 100 * ci0[["hi"]],
                  format.pval(pt$p.value, digits = 2)))
  report_value(paste0("Main text, para 14 (", label, ")"), "rate inside", k1 / n1)
  report_value(paste0("Main text, para 14 (", label, ")"), "  CI lower", ci1[["lo"]])
  report_value(paste0("Main text, para 14 (", label, ")"), "  CI upper", ci1[["hi"]])
  report_value(paste0("Main text, para 14 (", label, ")"), "rate outside", k0 / n0)
  report_value(paste0("Main text, para 14 (", label, ")"), "  CI lower ", ci0[["lo"]])
  report_value(paste0("Main text, para 14 (", label, ")"), "  CI upper ", ci0[["hi"]])
  report_value(paste0("Main text, para 14 (", label, ")"), "two-proportion P", pt$p.value, 6)
}

split_stats(big,  "C28 larger vs smaller effects")
split_stats(cell, "C28 screening cell vs remainder")

# Corrected-null benchmark within each subset (SI 2.7): the exact null of
# 00_functions.R applied to the subset's own marginal miss rates.
subset_null <- function(flag) null_exact(mean(miss_g[flag]), mean(miss_h[flag]))[["same"]]
report_value("SI 2.7", "corrected null, larger effects",  subset_null(big))
report_value("SI 2.7", "corrected null, smaller effects", subset_null(!big))
report_value("SI 2.7", "corrected null, screening cell",  subset_null(cell))


# --- 7. Attribution bound (SI 2.8) -------------------------------------------
# How much of the same-direction overlap could reflect inaccurate original
# estimates rather than forecasting failure? If an estimate's t-interval
# contains the true effect, a forecast outside that interval is provably not
# the truth, so that miss cannot be attributed to the estimate alone. Under
# calibrated standard errors the truth lies outside the t-interval with
# probability 2 * (1 - pnorm(t)), which therefore caps the share of
# comparisons whose same-direction joint miss the estimate alone could
# explain. The floor is how much of the observed rate must involve forecast
# error: (observed - cap) / observed. Two sensitivity variants accompany it:
# a binomial 97.5% upper limit on the realized number of miscalibrated
# intervals (comparisons cluster within studies, so this is indicative, not
# exact), and the floor if intervals failed at ten times their nominal rate.
attribution_bound <- function(t) {
  ss  <- mean(abs(wide$zg) > t & abs(wide$zh) > t & (sign(wide$zg) == sign(wide$zh)))
  cap <- 2 * (1 - pnorm(t))
  ub  <- qbinom(0.975, nrow(wide), cap) / nrow(wide)
  lab <- sprintf("threshold %.2f", t)
  report_value("SI 2.8", paste0(lab, ": same-side joint-miss rate"), ss)
  report_value("SI 2.8", paste0(lab, ": calibration cap"), cap)
  report_value("SI 2.8", paste0(lab, ": observed / cap"), ss / cap, 1)
  report_value("SI 2.8", paste0(lab, ": floor share, expectation"), (ss - cap) / ss)
  report_value("SI 2.8", paste0(lab, ": floor share, binomial 97.5% cap"), (ss - ub) / ss)
  report_value("SI 2.8", paste0(lab, ": floor share, tenfold miscalibration"),
               max(0, (ss - 10 * cap) / ss))
  invisible(NULL)
}
attribution_bound(1.96)
attribution_bound(3.29)


# --- Save objects for 05_figures.R ------------------------------------------

saveRDS(list(
  obs         = obs,
  nulls       = nulls,
  null_ci     = null_ci,      # finite-sample interval of the global null
  ss_ci       = ss_ci,        # and of the study-specific null
  p_mc        = p_mc,
  p_ss        = p_ss,
  p_boot      = p_boot,
  boot_ci     = boot_ci,      # study-bootstrap CI for the observed rates
  sb_primary  = sb_primary,   # slope + band for Fig. 1a
  sb_expert   = sb_expert,
  template    = template      # per-predictor shares, Wilson CIs, nulls
), file.path(derived_dir, "inference_objects.rds"))

message("06_inference.R complete. Tables in output/tables/, objects in data/derived/.")
