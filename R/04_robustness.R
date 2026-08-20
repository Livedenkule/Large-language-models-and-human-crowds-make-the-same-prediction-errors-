# ---------------------------------------------------------------------------
# 04_robustness.R
#
# The checks reported in the Supplementary Information: the study-specific
# null, alternative recalibrations and thresholds, study weighting, the
# authors' own sampling scheme, and study-level heterogeneity in the coupling.
#
# INPUT   data/derived/comparisons_primary.csv
# OUTPUT  output/tables/supplementary_table_1.csv
#         output/tables/study_level_summary.csv
#
# Wilhelmsen, Esfandiari & Gollwitzer
# ---------------------------------------------------------------------------

library(dplyr)
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


# --- 1. A fairer benchmark: letting each study set its own difficulty -------
# The corrected null of 03_main_results.R gives each predictor one miss rate for
# the whole archive, in effect assuming every study is equally hard. Studies
# plainly differ: some test effects almost anyone would call correctly, others
# are genuinely obscure. If a study is simply hard, both predictors will miss it
# more often even with wholly unrelated errors, and that extra overlap would
# show up as "shared error" when it is really "hard study".
#
# The fix is to run the same null but calibrate the two spreads within each
# study to that study's own miss rates, then pool. Study difficulty is then
# built into the benchmark itself.
#
# Because the shared errors themselves contribute to the within-study miss
# rates, part of the quantity being tested is built into this benchmark; it is
# in that sense stricter than the global null.

set.seed(21)
POOL   <- 1e6
pool_b <- rnorm(POOL)
pool_u <- rnorm(POOL)

rate_for <- function(s) mean(abs(s * pool_b - pool_u) > 1.96)

# Calibrate the spread so the simulated miss rate matches a target. Two guards:
# an error-free predictor still "misses" 5% of the time -- that is the floor --
# and the spread is capped at 40 for numerical safety.
calib_s <- function(target) {
  if (target <= rate_for(0))  return(0)
  if (target >= rate_for(40)) return(40)
  uniroot(function(s) rate_for(s) - target, c(0, 40), tol = 1e-4)$root
}

REP <- 4000                       # simulated replicates per real comparison

ss_null <- map_dfr(unique(wide$study), function(k) {
  m  <- wide$study == k
  nk <- sum(m)
  s1 <- calib_s(mean(miss_g[m]))
  s2 <- calib_s(mean(miss_h[m]))
  uu <- rnorm(nk * REP); b1 <- rnorm(nk * REP); b2 <- rnorm(nk * REP)
  z1 <- s1 * b1 - uu
  z2 <- s2 * b2 - uu
  j  <- abs(z1) > 1.96 & abs(z2) > 1.96
  data.frame(n = length(z1), joint = sum(j),
             same = sum(j & sign(z1) == sign(z2)))
})

ss_rates <- c(joint     = sum(ss_null$joint) / sum(ss_null$n),
              same_side = sum(ss_null$same)  / sum(ss_null$n),
              direction = sum(ss_null$same)  / sum(ss_null$joint))

print(round(rbind(observed = c(mean(joint), mean(same), sum(same) / sum(joint)),
                  `study-specific null` = ss_rates), 3))

report_value("SI 1.5", "study-specific null: joint rate",     ss_rates[["joint"]])
report_value("SI 1.5", "study-specific null: same-side rate", ss_rates[["same_side"]])
report_value("SI 1.5", "study-specific null: direction share", ss_rates[["direction"]])

# A hard study makes both predictors miss the same items, but gives them no
# reason to miss on the same side. The joint-miss benchmark therefore rises
# while the same-side benchmark barely moves and the direction benchmark falls,
# so the directional result -- the one that decides whether averaging can help
# -- survives this objection.
#
# Caveat: in a few studies a predictor misses so rarely that no spread
# reproduces its rate, since an error-free predictor still lands outside the
# interval 5% of the time. The spread is set to zero there, which makes the two
# simulated predictors perfectly aligned and so inflates the null -- the
# conservative direction.
degenerate <- map_lgl(unique(wide$study), function(k) {
  m <- wide$study == k
  calib_s(mean(miss_g[m])) == 0 || calib_s(mean(miss_h[m])) == 0
})
message("Studies where a spread of zero was assigned: ", sum(degenerate),
        " (", sum(wide$study %in% unique(wide$study)[degenerate]), " comparisons, ",
        sum(joint[wide$study %in% unique(wide$study)[degenerate]]), " joint misses).")

report_value("SI 1.5", "studies assigned zero spread", sum(degenerate), 0)


# --- 2. Robustness to recalibration and threshold ---------------------------
# Three alternatives, same conclusion: no recalibration at all; a more flexible
# monotone (isotonic) recalibration fitted with the same leave-one-study-out
# scheme; and different cut-offs for what counts as a miss.

zg_iso <- (iso_loo(wide$pg, wide$estimate, wide$study) - wide$estimate) / wide$std.error
zh_iso <- (iso_loo(wide$ph, wide$estimate, wide$study) - wide$estimate) / wide$std.error

supp_table_1 <- bind_rows(
  `No recalibration, 1.96`       = summarise_shared(wide$zg_raw, wide$zh_raw),
  `Linear recalibration, 1.96`   = summarise_shared(wide$zg, wide$zh),
  `Isotonic recalibration, 1.96` = summarise_shared(zg_iso, zh_iso),
  `Linear, threshold 1.64`       = summarise_shared(wide$zg, wide$zh, 1.64),
  `Linear, threshold 2.58`       = summarise_shared(wide$zg, wide$zh, 2.58),
  `Linear, threshold 3.29`       = summarise_shared(wide$zg, wide$zh, 3.29),
  .id = "analysis"
) %>%
  mutate(across(where(is.numeric), ~ round(., 3)))

print(as.data.frame(supp_table_1))
write.csv(supp_table_1, file.path(table_dir, "supplementary_table_1.csv"), row.names = FALSE)

report_value("Main text, para 7", "GPT-4 miss rate (unrecalibrated)",
             supp_table_1$miss_g[supp_table_1$analysis == "No recalibration, 1.96"])
report_value("Main text, para 7", "crowd miss rate (unrecalibrated)",
             supp_table_1$miss_h[supp_table_1$analysis == "No recalibration, 1.96"])

# Setting standardized errors aside entirely, the two predictors' errors in the
# original percentage-point metric correlate directly.
er_g <- wide$zg * wide$std.error
er_h <- wide$zh * wide$std.error
report_value("SI 2.1", "error correlation, recalibrated (pp metric)", cor(er_g, er_h), 2)
report_value("SI 2.1", "error correlation, unrecalibrated (pp metric)",
             cor(wide$pg - wide$estimate, wide$ph - wide$estimate), 2)


# --- 3. Study weighting -----------------------------------------------------
# Studies contribute very unequal numbers of comparisons, so the headline rates
# are dominated by the condition-rich studies. Giving each study equal weight
# raises both rates, which means the reported figures are the conservative
# choice rather than the flattering one.

by_study <- wide %>%
  group_by(study) %>%
  summarise(n     = n(),
            joint = mean(abs(zg) > 1.96 & abs(zh) > 1.96),
            same  = mean(abs(zg) > 1.96 & abs(zh) > 1.96 & sign(zg) == sign(zh)),
            .groups = "drop")

message("Comparisons per study range from ", min(by_study$n), " to ", max(by_study$n), ".")
message("Study-weighted joint rate ", round(100 * mean(by_study$joint), 1),
        "%, same-side ", round(100 * mean(by_study$same), 1), "%.")

report_value("SI 2.2", "study-weighted joint rate", mean(by_study$joint))
report_value("SI 2.2", "study-weighted same-side rate", mean(by_study$same))
report_value("SI 2.2", "smallest study, comparisons", min(by_study$n), 0)
report_value("SI 2.2", "largest study, comparisons", max(by_study$n), 0)

write.csv(by_study, file.path(table_dir, "study_level_summary.csv"), row.names = FALSE)


# --- 4. The authors' own sampling scheme ------------------------------------
# This reanalysis uses all 1,678 available comparisons. The authors instead draw
# one outcome and one reference condition per study, compare every other
# condition against that reference, and repeat the draw 32 times. Reproducing
# their scheme answers the reasonable question of whether the result depends on
# using the full set.
#
# A draw yields about 235 unordered pairs, which is the authors' reported 469
# effects counted in both orientations.

set.seed(22)
idx_by_study <- split(seq_len(nrow(wide)), wide$study)

one_draw <- function() {
  keep <- unlist(lapply(idx_by_study, function(ii) {
    ocs <- unique(wide$outcome.name[ii])
    oc  <- ocs[sample.int(length(ocs), 1)]           # one outcome per study
    jj  <- ii[wide$outcome.name[ii] == oc]
    cds <- unique(c(wide$reference_condition[jj], wide$condition.name[jj]))
    ref <- cds[sample.int(length(cds), 1)]           # one reference condition
    jj[wide$reference_condition[jj] == ref | wide$condition.name[jj] == ref]
  }), use.names = FALSE)

  g <- wide$zg[keep]; h <- wide$zh[keep]
  j <- abs(g) > 1.96 & abs(h) > 1.96
  c(n = length(keep), joint = mean(j), same = mean(j & sign(g) == sign(h)))
}

draws32  <- t(replicate(32,  one_draw()))
draws500 <- t(replicate(500, one_draw()))

sampling <- rbind(
  `all 1,678 comparisons`   = c(NA, mean(joint), mean(same)),
  `their scheme, 32 draws`  = c(median(draws32[,  "n"]),
                                median(draws32[,  "joint"]), median(draws32[,  "same"])),
  `their scheme, 500 draws` = c(median(draws500[, "n"]),
                                median(draws500[, "joint"]), median(draws500[, "same"]))
)
colnames(sampling) <- c("effects_per_draw", "joint", "same_side")
print(round(sampling, 3))

report_value("SI 2.2", "authors' scheme, 32 draws: joint",  median(draws32[, "joint"]))
report_value("SI 2.2", "authors' scheme, 32 draws: same-side", median(draws32[, "same"]))
report_value("SI 2.2", "authors' scheme, 500 draws: joint", median(draws500[, "joint"]))
report_value("SI 2.2", "authors' scheme, 500 draws: same-side", median(draws500[, "same"]))
report_value("SI 2.2", "comparisons retained per draw (median)", median(draws500[, "n"]), 0)

# The spread across draws measures sensitivity to the sampling step. It is not
# a confidence interval; the intervals in the manuscript come from the study
# bootstrap in 03_main_results.R.
report_value("SI 2.2", "s.d. across draws, joint", sd(draws500[, "joint"]))
report_value("SI 2.2", "s.d. across draws, same-side", sd(draws500[, "same"]))


# --- 5. How strongly are the two errors coupled, study by study? ------------

fit    <- lm(zh ~ zg, data = wide)
se_cl  <- cluster_se(fit, wide$study)
se_iid <- summary(fit)$coefficients[2, 2]
slope  <- coef(fit)[2]

message("Slope ", round(slope, 2),
        "; cluster-robust 95% CI ", round(slope - 1.96 * se_cl, 2),
        " to ", round(slope + 1.96 * se_cl, 2),
        "; i.i.d. CI ", round(slope - 1.96 * se_iid, 2),
        " to ", round(slope + 1.96 * se_iid, 2), ".")

report_value("SI 2.3", "pooled slope of crowd error on GPT-4 error", slope, 2)
report_value("SI 2.3", "slope CI lower (clustered by study)", slope - 1.96 * se_cl, 2)
report_value("SI 2.3", "slope CI upper (clustered by study)", slope + 1.96 * se_cl, 2)
report_value("SI 2.3", "slope CI lower (i.i.d.)", slope - 1.96 * se_iid, 2)
report_value("SI 2.3", "slope CI upper (i.i.d.)", slope + 1.96 * se_iid, 2)

# Share of the co-error's variance lying between studies (one-way ANOVA).
w   <- wide$zg * wide$zh
gm  <- mean(w)
ssb <- sum(tapply(w, wide$study, function(x) length(x) * (mean(x) - gm)^2))
report_value("SI 2.3", "percent of co-error variance between studies",
             100 * ssb / sum((w - gm)^2), 1)

slopes <- wide %>%
  group_by(study) %>%
  filter(n() >= 5, sd(zg) > 0) %>%
  summarise(slope = coef(lm(zh ~ zg))[2], .groups = "drop")

report_value("SI 2.3", "studies with at least five comparisons", nrow(slopes), 0)
report_value("SI 2.3", "median per-study slope", median(slopes$slope), 2)
report_value("SI 2.3", "per-study slope, 25th percentile", quantile(slopes$slope, .25), 2)
report_value("SI 2.3", "per-study slope, 75th percentile", quantile(slopes$slope, .75), 2)

message("04_robustness.R complete.")
