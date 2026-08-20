# ---------------------------------------------------------------------------
# 03_main_results.R
#
# The results reported in the main text: miss rates, joint and same-side
# misses, the corrected null, the permutation test, the study bootstrap, the
# performance of the authors' combined predictor, the cross-model comparison,
# and the secondary archive of megastudies.
#
# INPUT   data/derived/comparisons_primary.csv
#         data/derived/comparisons_all_models.csv
#         data/derived/comparisons_archive2.csv
# OUTPUT  output/tables/main_results.csv
#         output/tables/supplementary_table_2.csv
#
# Runs from the derived files alone: the original capsule is not needed.
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


# --- 1. Miss rates and shared misses ----------------------------------------
# A miss is a prediction outside the observed effect's own 95% confidence
# interval. A joint miss is a comparison both predictors miss. A same-side miss
# is a joint miss on which both err in the same direction -- both too high or
# both too low relative to what the experiment found.
#
# The 1,678 are all available comparisons, not a random or independent sample:
# pairs within a study share participants and often share a condition. That is
# why every interval below resamples studies rather than comparisons, and why
# 04_robustness.R also reports study-weighted rates and the authors' own
# sampling scheme.

miss_g <- abs(wide$zg) > 1.96
miss_h <- abs(wide$zh) > 1.96
joint  <- miss_g & miss_h
same   <- joint & (sign(wide$zg) == sign(wide$zh))

core <- data.frame(
  gpt4_miss_rate     = mean(miss_g),
  crowd_miss_rate    = mean(miss_h),
  joint_misses_n     = sum(joint),
  joint_rate         = mean(joint),
  same_side_n        = sum(same),
  same_side_rate     = mean(same),
  direction_given_joint = sum(same) / sum(joint)
)
print(round(core, 3))

report_value("Main text, para 7", "GPT-4 miss rate (recalibrated)", mean(miss_g))
report_value("Main text, para 7", "crowd miss rate (recalibrated)", mean(miss_h))
report_value("Main text, para 8", "joint misses (n)", sum(joint), 0)
report_value("Main text, para 8", "joint miss rate", mean(joint))
report_value("Main text, para 8", "same-side misses (n)", sum(same), 0)
report_value("Main text, para 8", "same-side rate", mean(same))
report_value("Main text, para 8", "direction share given joint miss", sum(same) / sum(joint))
report_value("SI 2.1", "Spearman correlation of standardized errors",
             cor(wide$zg, wide$zh, method = "spearman"), 2)


# --- 2. Are the errors dependent? A test with no assumptions ----------------
# Multiply the two errors comparison by comparison: w = zg * zh. Same-side
# errors give a positive product, opposite-side a negative one. If the two
# predictors were unrelated the products would balance and the mean would be
# zero. Shuffling the pairing -- matching GPT-4's error on one comparison with
# the crowd's on a random other -- destroys any true dependence while leaving
# each predictor's own error distribution untouched, and so supplies the
# "unrelated" benchmark without assuming anything.

set.seed(1)
w_obs  <- mean(wide$zg * wide$zh)
w_perm <- replicate(2000, mean(wide$zg * sample(wide$zh)))

message("Observed mean co-error ", round(w_obs, 2),
        " against a permutation mean of ", round(mean(w_perm), 3),
        " (s.d. ", round(sd(w_perm), 3), ").")

report_value("SI 2.4", "observed mean co-error", w_obs, 2)
report_value("SI 2.4", "permutation mean co-error", mean(w_perm))
report_value("SI 2.4", "permutation s.d.", sd(w_perm))

# Not one anomalous study: the study-level mean is positive in most studies.
# The Wilcoxon test treats the study, not the comparison, as the unit -- the
# conservative choice, since comparisons inside a study share participants.
per_study <- wide %>% group_by(study) %>% summarise(w = mean(zg * zh), .groups = "drop")
wilcox_p  <- wilcox.test(per_study$w)$p.value

message("Mean co-error positive in ", sum(per_study$w > 0), " of ",
        nrow(per_study), " studies (Wilcoxon P = ", signif(wilcox_p, 2), ").")

report_value("SI 2.4", "studies with positive mean co-error", sum(per_study$w > 0), 0)
report_value("SI 2.4", "studies total", nrow(per_study), 0)
report_value("SI 2.4", "Wilcoxon P", wilcox_p, 10)


# --- 3. The corrected null --------------------------------------------------
# Both predictors are graded against the same estimated effect, and that
# estimate carries sampling noise of its own. When an experiment's estimate
# happens to land high by chance, every prediction is "too low" relative to it,
# through no fault of the predictors. So even wholly unrelated predictions show
# some same-side clustering, and naive benchmarks -- 50% direction, or the
# product of the two miss rates -- would flatter the result.
#
# The corrected null keeps that shared term and makes the predictions otherwise
# independent, calibrating each one's spread to reproduce the observed marginal
# miss rate. Whatever the data show beyond it is genuine dependence between the
# predictions. See null_exact() in 00_functions.R.

nulls <- null_exact(mean(miss_g), mean(miss_h))

message("Corrected null: joint ", round(100 * nulls[["joint"]], 1),
        "%, same-side ", round(100 * nulls[["same"]], 1),
        "%, direction ", round(100 * nulls[["direction"]], 1), "%.")

report_value("Main text, para 9", "null joint rate", nulls[["joint"]])
report_value("Main text, para 9", "null same-side rate", nulls[["same"]])
report_value("Main text, para 9", "null direction share", nulls[["direction"]])

# The benchmark that ignores the shared estimate, shown as the third bar of
# Fig. 1b, is half the product of the marginal rates.
naive_same <- mean(miss_g) * mean(miss_h) / 2
report_value("Fig. 1b", "naive same-side null (ignores shared estimate)", naive_same)

# Does the null's normality matter? The real errors have heavier tails, so the
# two decisive benchmarks are recomputed with a Student t(3) and with the
# empirical error shape. Heavier tails lower the benchmarks slightly, so the
# normal values quoted in the manuscript are the conservative ones.
eg <- (wide$zg - mean(wide$zg)) / sd(wide$zg)
eh <- (wide$zh - mean(wide$zh)) / sd(wide$zh)

tail_variants <- rbind(
  normal       = 100 * nulls[c("same", "direction")],
  `Student t3` = null_exact_tail(mean(miss_g), mean(miss_h), ptail = t3_tail),
  empirical    = null_exact_tail(
    mean(miss_g), mean(miss_h),
    ptail = function(q) (empirical_tail(eg)(q) + empirical_tail(eh)(q)) / 2
  )
)
print(round(tail_variants, 1))

report_value("SI 1.4", "same-side null, t(3)", tail_variants["Student t3", "same"], 1)
report_value("SI 1.4", "same-side null, empirical", tail_variants["empirical", "same"], 1)
report_value("SI 1.4", "direction null, t(3)", tail_variants["Student t3", "direction"], 1)
report_value("SI 1.4", "direction null, empirical", tail_variants["empirical", "direction"], 1)


# --- 4. Uncertainty that respects the studies -------------------------------
# Comparisons inside a study are not independent: they share participants and
# conditions. Rather than pretend to 1,678 independent observations, the 70
# studies are resampled with replacement and each rate recomputed.

set.seed(3)
studies <- unique(wide$study)
boot <- replicate(4000, {
  s  <- sample(studies, replace = TRUE)
  d  <- map_dfr(s, ~ wide[wide$study == .x, ])
  jg <- abs(d$zg) > 1.96 & abs(d$zh) > 1.96
  sm <- jg & sign(d$zg) == sign(d$zh)
  c(joint = mean(jg), same = mean(sm), dir = sum(sm) / sum(jg))
})
boot_ci <- apply(boot, 1, quantile, c(.025, .975))
print(round(boot_ci, 3))

report_value("SI 1.6", "joint rate 95% CI lower", boot_ci[1, "joint"])
report_value("SI 1.6", "joint rate 95% CI upper", boot_ci[2, "joint"])
report_value("Main text, para 8", "same-side rate 95% CI lower", boot_ci[1, "same"])
report_value("Main text, para 8", "same-side rate 95% CI upper", boot_ci[2, "same"])
report_value("Main text, para 8", "direction share 95% CI lower", boot_ci[1, "dir"])
report_value("Main text, para 8", "direction share 95% CI upper", boot_ci[2, "dir"])


# --- 5. Does averaging actually fail? ---------------------------------------
# No assumptions needed: the authors propose averaging model and crowd, that
# average was built and recalibrated identically in 01_build_comparisons.R, so
# its behaviour on the comparisons in question can simply be read off.

combined_perf <- data.frame(
  same_side_combined_still_outside = mean(abs(wide$zc[same]) > 1.96),
  opposite_side_combined_outside   = mean(abs(wide$zc[joint & !same]) > 1.96)
)
print(round(combined_perf, 3))

report_value("Main text, para 8", "combined recovers | opposite-side joint miss",
             1 - combined_perf$opposite_side_combined_outside)
report_value("Main text, para 8", "combined still outside | same-side joint miss",
             combined_perf$same_side_combined_still_outside)


# --- 6. Where the shared misses fall ----------------------------------------
# A user screening interventions or flagging effects for replication acts where
# the predictions look large and agreed. Scattered at random, shared misses
# would be an accuracy tax; concentrated where a user would act, they are a
# selection hazard.

big <- abs(wide$estimate) >= 0.05          # the authors' 5-percentage-point threshold

report_value("SI 2.7", "comparisons with larger observed effects", sum(big), 0)
report_value("SI 2.7", "comparisons with smaller observed effects", sum(!big), 0)
report_value("Main text, para 14", "same-side rate | larger effects", mean(same[big]))
report_value("Main text, para 14", "same-side rate | smaller effects", mean(same[!big]))

# The screening configuration: top quintile of predicted magnitude, the two
# predictions agreeing in sign, and disagreeing by less than the median.
mag  <- abs((wide$pg + wide$ph) / 2)
gap  <- abs(wide$pg - wide$ph)
cell <- mag >= quantile(mag, .8) & sign(wide$pg) == sign(wide$ph) & gap <= median(gap)

message("Screening cell: ", sum(cell), " comparisons, ",
        round(100 * mean(same[cell]), 1), "% same-side joint misses, against ",
        round(100 * mean(same[!cell]), 1), "% of the remainder.")

report_value("Main text, para 14", "screening cell size", sum(cell), 0)
report_value("Main text, para 14", "same-side rate inside screening cell", mean(same[cell]))
report_value("SI 2.7", "same-side rate outside screening cell", mean(same[!cell]))

# Direction of the shared errors: both predictors erring towards zero is shared
# conservatism about large effects, which linear recalibration does not remove.
toward_zero <- mean(sign(wide$zg[same]) != sign(wide$estimate[same]))
message("Both err towards zero on ", round(100 * toward_zero, 1),
        "% of same-side joint misses.")

report_value("Main text, para 15", "share of same-side misses erring towards zero",
             toward_zero)

report_value("SI 2.3", "studies with a joint miss", n_distinct(wide$study[joint]), 0)
report_value("SI 2.3", "studies with a same-side miss", n_distinct(wide$study[same]), 0)


# --- 7. Can another model, or an expert, be the independent check? ----------

zwide <- read.csv(file.path(derived_dir, "comparisons_all_models.csv"),
                  check.names = FALSE) %>%
  select(study, outcome.name, reference_condition, condition.name, model, z_recal) %>%
  pivot_wider(names_from = model, values_from = z_recal)

others <- c("deepseek/deepseek-chat-v3-0324", "openai/gpt-oss-120b", "gpt-3.5-turbo",
            "google/gemma-3-27b-it", "babbage-002", "davinci-002", "human")
gp <- zwide$`gpt-4`

cross_model <- map_dfr(others, function(m) {
  zo  <- zwide[[m]]
  ok  <- !is.na(zo) & !is.na(gp)
  gmr <- mean(abs(gp[ok]) > 1.96)          # GPT-4's miss rate on this overlap
  omr <- mean(abs(zo[ok]) > 1.96)
  ne  <- null_exact(gmr, omr)              # exact corrected null for this pair
  gmiss <- ok & abs(gp) > 1.96
  data.frame(
    predictor      = m,
    miss_rate      = omr,
    also_same_side = mean((abs(zo) > 1.96 & sign(zo) == sign(gp))[gmiss]),
    corrected_null = ne[["same"]] / gmr,
    error_corr     = cor(gp[ok], zo[ok])
  )
}) %>%
  mutate(across(where(is.numeric), ~ round(., 2)))

print(cross_model)
write.csv(cross_model, file.path(table_dir, "supplementary_table_2.csv"), row.names = FALSE)

for (i in seq_len(nrow(cross_model))) {
  report_value("SI 2.5",
               paste0(cross_model$predictor[i], ": also misses same side"),
               cross_model$also_same_side[i], 2)
}


# --- 8. The secondary archive of megastudies --------------------------------

a2_path <- file.path(derived_dir, "comparisons_archive2.csv")
if (file.exists(a2_path)) {
  a2 <- read.csv(a2_path, check.names = FALSE)
  a2 <- a2[!is.na(a2$z_recal_gpt) & !is.na(a2$z_recal_expert), ]

  gm2   <- abs(a2$z_recal_gpt)    > 1.96
  em2   <- abs(a2$z_recal_expert) > 1.96
  both2 <- gm2 & em2
  same2 <- both2 & (sign(a2$z_recal_gpt) == sign(a2$z_recal_expert))
  n2    <- null_exact(mean(gm2), mean(em2))

  message("Secondary archive: GPT-4 misses ", round(100 * mean(gm2), 1),
          "%, experts ", round(100 * mean(em2), 1),
          "%; experts also miss ", round(100 * sum(both2) / sum(gm2), 1),
          "% of GPT-4's misses, ", round(100 * sum(same2) / sum(both2), 1),
          "% of those on the same side.")

  report_value("Main text, para 12", "expert miss rate", mean(em2))
  report_value("Main text, para 12", "GPT-4 miss rate (secondary archive)", mean(gm2))
  report_value("Main text, para 12", "experts also miss | GPT-4 misses", sum(both2) / sum(gm2))
  report_value("Main text, para 12", "same side | both miss", sum(same2) / sum(both2))
  report_value("Main text, para 12", "slope of expert error on model error",
               coef(lm(a2$z_recal_expert ~ a2$z_recal_gpt))[2], 2)
  report_value("SI 2.6", "null: experts also miss | GPT-4 misses", n2[["joint"]] / mean(gm2))
  report_value("SI 2.6", "null: same side | both miss", n2[["direction"]])
} else {
  warning("comparisons_archive2.csv not found; skipping the secondary archive. ",
          "Run 02_build_archive2.R first.")
}


# --- Save -------------------------------------------------------------------

write.csv(cbind(core, combined_perf), file.path(table_dir, "main_results.csv"),
          row.names = FALSE)

# Objects used by 04_robustness.R and 05_figures.R when run in one session.
saveRDS(list(wide = wide, miss_g = miss_g, miss_h = miss_h,
             joint = joint, same = same, nulls = nulls, boot_ci = boot_ci,
             cross_model = cross_model),
        file.path(derived_dir, "main_results_objects.rds"))

message("03_main_results.R complete.")
