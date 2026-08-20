# ---------------------------------------------------------------------------
# 01_build_comparisons.R
#
# Builds the comparison-level dataset for the primary archive from the authors'
# Code Ocean capsule (9843791): observed effects, model and crowd predictions,
# standardized errors, and leave-one-study-out recalibration.
#
# INPUT   data/raw/rct_responses.RDS
#         data/raw/llm_responses.RDS
#         data/raw/forecasting_responses.RDS
# OUTPUT  data/derived/comparisons_primary.csv
#
# This is the only script that touches the original data. Everything
# downstream runs from the derived file alone, so the analysis can be
# reproduced without access to the capsule.
#
# Wilhelmsen, Esfandiari & Gollwitzer
# ---------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(here)

source(here("R", "00_functions.R"))

raw_dir <- here("data", "raw")
out_dir <- here("data", "derived")


# --- 1. Load ----------------------------------------------------------------
# rct_responses  : individual answers of the ~120,000 experiment participants,
#                  the ground truth both predictors are trying to guess.
# llm_responses  : one row per model call (about 120 demographic-profile
#                  prompts per condition, each value already the mean of five
#                  samples), for each of seven models.
# forecasting_responses : one row per forecaster x condition, median 29
#                  forecasters per question.

rct <- readRDS(file.path(raw_dir, "rct_responses.RDS"))
llm <- readRDS(file.path(raw_dir, "llm_responses.RDS"))
fc  <- readRDS(file.path(raw_dir, "forecasting_responses.RDS"))

message("Loaded: ", sum(map_int(rct$data, nrow)), " participant responses; ",
        nrow(llm), " model calls; ", nrow(fc), " forecasts from ",
        n_distinct(fc$PROLIFIC_PID), " forecasters.")

report_value("SI 1.1", "participant responses (all)", sum(map_int(rct$data, nrow)), 0)
report_value("SI 1.1", "participant responses (non-missing)",
             sum(map_int(rct$data, ~ sum(!is.na(.x$y)))), 0)
report_value("SI 1.1", "model responses", nrow(llm), 0)
report_value("SI 1.1", "forecasts", nrow(fc), 0)
report_value("SI 1.1", "forecasters", n_distinct(fc$PROLIFIC_PID), 0)
report_value("SI 1.1", "studies in raw archive", n_distinct(rct$study), 0)
report_value("SI 1.1", "study-outcome cells in raw archive", nrow(rct), 0)


# --- 2. Observed effects ----------------------------------------------------
# Every outcome is rescaled to the unit interval using its own scale bounds, so
# that effects measured on different answer scales are comparable. Within each
# study-outcome cell, lm(y ~ condition) returns every condition-versus-reference
# difference together with its standard error in one step; the coefficient is
# literally the difference in means. This follows the authors' own pipeline, so
# the observed effects here are theirs.
#
# Condition names are recoded to plain ASCII before the regression and mapped
# back afterwards. Some names contain typographic quotes, and on a machine with
# a non-UTF-8 locale R silently rewrites such characters inside coefficient
# names, which would break the merge with the prediction files. The recode makes
# the result identical on any machine.

rct_effects <- rct %>%
  pmap_dfr(function(study, outcome.name, outcome.min, outcome.max, data, ...) {
    d <- data %>%
      ungroup() %>%
      mutate(y = (as.numeric(y) - outcome.min) / (outcome.max - outcome.min))

    conds <- unique(d$condition.name)
    code  <- setNames(paste0("k", seq_along(conds)), conds)
    d$cc  <- code[d$condition.name]

    map_dfr(conds, function(ref) {
      fit <- lm(y ~ cc, data = d %>% mutate(cc = forcats::fct_relevel(cc, code[[ref]])))
      broom::tidy(fit) %>%
        transmute(cc = str_match(term, "^cc(k\\d+)$")[, 2], estimate, std.error) %>%
        filter(!is.na(cc)) %>%
        mutate(condition.name      = names(code)[match(cc, code)],
               study               = study,
               outcome.name        = outcome.name,
               reference_condition = ref,
               .before = 1) %>%
        select(-cc)
    })
  })


# --- 3. Predictions ---------------------------------------------------------
# Model: the predicted mean of a condition is the weighted average of its
# prompt-level responses. The authors over-sampled some demographic groups to
# study subgroups; their weights undo that over-sampling so the ensemble
# represents the US population. Reversed answer scales are un-flipped first.
# Both steps follow their code.
#
# Crowd: each forecaster predicted a condition's average answer on a 0-100
# slider; the crowd prediction is the unweighted mean, divided by 100.
#
# A predicted effect is then the predicted mean of one condition minus that of
# the other, the same subtraction the experiment performs with real people.

llm_means <- llm %>%
  mutate(
    expectation = ifelse(scale_flip,
                         outcome_scale_max - (expectation - outcome_scale_min),
                         expectation),
    y = (expectation - outcome_scale_min) / (outcome_scale_max - outcome_scale_min)
  ) %>%
  filter(!is.na(y)) %>%
  group_by(model, study, outcome.name, condition.name) %>%
  summarise(pred_mean = weighted.mean(y, w = weight), .groups = "drop")

fc_means <- fc %>%
  mutate(y = value / 100) %>%
  semi_join(distinct(llm, study, outcome.name), by = c("study", "outcome.name")) %>%
  group_by(study, outcome.name, condition.name) %>%
  summarise(pred_mean = mean(y), .groups = "drop") %>%
  mutate(model = "human")

cond_means <- bind_rows(llm_means, fc_means)

pred_contrasts <- cond_means %>%
  group_by(model, study, outcome.name) %>%
  group_modify(function(d, k) {
    crossing(
      d %>% select(condition.name, pred_mean),
      d %>% select(reference_condition = condition.name, ref_mean = pred_mean)
    ) %>%
      filter(condition.name != reference_condition) %>%
      transmute(condition.name, reference_condition,
                prediction = pred_mean - ref_mean)
  }) %>%
  ungroup()


# --- 4. One row per comparison ---------------------------------------------
# Comparing A with B and B with A is the same comparison with the sign flipped,
# so each unordered pair is kept once.
#
# Ordering names with '<' depends on the machine's locale: in a plain C locale R
# returns NA for names containing typographic quotes, silently dropping those
# comparisons. A fixed ICU collator (English) is identical on every machine and
# matches the orientation of the deposited comparison-level dataset.
#
# The rates are unaffected by orientation in any case, because flipping a pair
# flips the sign of the observed effect and of both predictions together,
# leaving every |z| and every sign agreement unchanged.

en <- stringi::stri_opts_collator(locale = "en_US")

comparisons <- rct_effects %>%
  inner_join(pred_contrasts,
             by = c("study", "outcome.name", "condition.name", "reference_condition")) %>%
  filter(stringi::stri_cmp_lt(reference_condition, condition.name, opts_collator = en))

# The authors' proposed combination of model and crowd is their simple average;
# it is built the same way here so that it can be tested on the same footing.
combined <- comparisons %>%
  filter(model %in% c("gpt-4", "human")) %>%
  select(study, outcome.name, condition.name, reference_condition,
         estimate, std.error, model, prediction) %>%
  pivot_wider(names_from = model, values_from = prediction) %>%
  filter(!is.na(`gpt-4`), !is.na(human)) %>%
  mutate(model = "combined", prediction = (`gpt-4` + human) / 2) %>%
  select(-`gpt-4`, -human)

comparisons <- bind_rows(comparisons, combined)


# --- 5. Standardized errors and recalibration -------------------------------
# A prediction five points off means little in a small noisy study and a lot in
# a large precise one, so each error is expressed in units of the observed
# effect's own standard error:
#
#     z = (prediction - estimate) / s.e.(estimate)
#
# |z| > 1.96 places the prediction outside the experiment's own 95% confidence
# interval: the experiment itself is saying the prediction is wrong. That is a
# "miss".
#
# The authors report that predictions overshoot true effects by roughly a factor
# of two and that a linear rescaling fixes it. Applying their remedy before
# measuring anything forestalls the objection that the shared error is just that
# one shared, already-known bias. The line is fitted leaving out one study at a
# time, so no study helps correct itself and every predictor is treated alike.
#
# The shared error is not created by this step: it is already present in the
# raw predictions (see 04_robustness.R). Recalibration sharpens the picture; it
# does not manufacture it.

comparisons <- comparisons %>%
  mutate(z_raw = (prediction - estimate) / std.error) %>%
  group_by(model) %>%
  group_modify(function(d, k) {
    d$pred_recal <- linear_loo(d$prediction, d$estimate, d$study)
    d
  }) %>%
  ungroup() %>%
  mutate(z_recal = (pred_recal - estimate) / std.error)

# Pooled recalibration lines, quoted in SI 1.3. The leave-one-study-out fits
# vary only slightly around these: intercepts essentially zero, slopes about
# 0.5, i.e. both predictors overstate effects by roughly a factor of two.
pooled_lines <- comparisons %>%
  filter(model %in% c("gpt-4", "human")) %>%
  group_by(model) %>%
  group_modify(~ broom::tidy(lm(estimate ~ prediction, data = .x))) %>%
  ungroup() %>%
  transmute(model, term, estimate = round(estimate, 4))

print(pooled_lines)


# --- 6. Wide comparison-level file ------------------------------------------
# One row per comparison, with GPT-4's, the crowd's and the combined
# predictor's errors side by side. This is the analysed set: all comparisons
# for which both GPT-4 and the crowd made a prediction.

wide <- comparisons %>%
  filter(model %in% c("gpt-4", "human", "combined")) %>%
  select(study, outcome.name, reference_condition, condition.name,
         estimate, std.error, prediction, model, z_raw, z_recal) %>%
  pivot_wider(names_from = model, values_from = c(prediction, z_raw, z_recal)) %>%
  filter(!is.na(`z_recal_gpt-4`), !is.na(z_recal_human)) %>%
  rename(zg = `z_recal_gpt-4`,  zh = z_recal_human,  zc = z_recal_combined,
         pg = `prediction_gpt-4`, ph = prediction_human,
         zg_raw = `z_raw_gpt-4`, zh_raw = z_raw_human)

stopifnot(nrow(wide) == 1678L)

message("Built ", nrow(wide), " comparisons across ", n_distinct(wide$study),
        " studies and ", nrow(distinct(wide, study, outcome.name)),
        " study-outcome cells.")

report_value("SI 1.2", "analysed comparisons", nrow(wide), 0)
report_value("SI 1.2", "studies", n_distinct(wide$study), 0)
report_value("SI 1.2", "study-outcome cells", nrow(distinct(wide, study, outcome.name)), 0)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(wide, file.path(out_dir, "comparisons_primary.csv"), row.names = FALSE)

# The full long-format file retains the other six models, needed for the
# cross-model analysis in 03_main_results.R.
write.csv(
  comparisons %>%
    select(study, outcome.name, reference_condition, condition.name,
           estimate, std.error, model, prediction, z_raw, z_recal),
  file.path(out_dir, "comparisons_all_models.csv"), row.names = FALSE
)

message("Wrote comparisons_primary.csv and comparisons_all_models.csv to data/derived/")
