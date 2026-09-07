# ---------------------------------------------------------------------------
# 05_figures.R
#
# R diagnostic versions of Figures 1 and 2 (F1_same_effects.png,
# F2_shared_blind_spot.png), drawn from the derived comparison-level data.
# These are NOT the manuscript figures: fig1-3.png in output/figures/ are
# produced by scripts/make_figures.py (the same code is embedded,
# eval-guarded, in docs/analysis_report.Rmd).
#
# Revised for the supervisor comments of 2026-08-26:
#   C5   Panel labels and axis wording ("forecasting error", identical units
#        on both axes of Fig. 1a, clearer Panel c title).
#   C6   Error bars on every bar; fitted slope with a confidence band in
#        Fig. 1a (cluster-robust, clustered by study).
#   C19  Fig. 2 says explicitly which bars are lay forecasters and which are
#        domain experts, and the text now reports the crowd's 51%.
#   C20  "inherits" dropped; wording is now "errs in the same direction".
#   C21  "Also wrong, on the same side as GPT-4" replaced.
#   C24  The expert direction result gets its own Panel c instead of being a
#        text tack-on in Panel b.
#
# INPUT   data/derived/comparisons_primary.csv
#         data/derived/comparisons_all_models.csv
#         data/derived/comparisons_archive2.csv
#         data/derived/inference_objects.rds   (optional; recomputed if absent)
# OUTPUT  output/figures/F1_same_effects.png
#         output/figures/F2_shared_blind_spot.png
#
# Wilhelmsen, Esfandiari & Gollwitzer
# ---------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(patchwork)
library(here)

source(here("R", "00_functions.R"))

derived_dir <- here("data", "derived")
fig_dir     <- here("output", "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

wide <- read.csv(file.path(derived_dir, "comparisons_primary.csv"), check.names = FALSE)

miss_g <- abs(wide$zg) > 1.96
miss_h <- abs(wide$zh) > 1.96
joint  <- miss_g & miss_h
same   <- joint & (sign(wide$zg) == sign(wide$zh))

nulls <- null_exact(mean(miss_g), mean(miss_h))

# Uncertainty inputs (C4, C6). 06_inference.R computes these at full
# precision and saves them; if it has not been run, a lighter version is
# recomputed here so the figures always carry their error bars.
inf_path <- file.path(derived_dir, "inference_objects.rds")
if (file.exists(inf_path)) {
  inf <- readRDS(inf_path)
} else {
  message("inference_objects.rds not found; recomputing intervals (reduced draws). ",
          "Run 06_inference.R for the versions quoted in the text.")
  set.seed(41)
  sd_for <- function(target) uniroot(function(s)
    2 * (1 - pnorm(1.96 / sqrt(s^2 + 1))) - target, c(0.01, 40))$root
  sg <- sd_for(mean(miss_g)); sh <- sd_for(mean(miss_h))
  nd <- t(replicate(4000, {
    u <- rnorm(nrow(wide))
    z1 <- sg * rnorm(nrow(wide)) - u; z2 <- sh * rnorm(nrow(wide)) - u
    j <- abs(z1) > 1.96 & abs(z2) > 1.96
    c(joint = mean(j), same_side = mean(j & sign(z1) == sign(z2)),
      direction = if (sum(j)) sum(j & sign(z1) == sign(z2)) / sum(j) else NA_real_)
  }))
  set.seed(43)
  idx_by_study <- split(seq_len(nrow(wide)), wide$study)
  bt <- replicate(2000, {
    ii <- unlist(idx_by_study[sample(names(idx_by_study), replace = TRUE)],
                 use.names = FALSE)
    j <- joint[ii]; sm <- same[ii]
    c(joint = mean(j), same_side = mean(sm),
      direction = if (sum(j)) sum(sm) / sum(j) else NA_real_)
  })
  inf <- list(
    null_ci    = apply(nd, 2, quantile, c(.025, .975), na.rm = TRUE),
    boot_ci    = apply(bt, 1, quantile, c(.025, .975), na.rm = TRUE),
    sb_primary = slope_band(wide$zg, wide$zh, wide$study,
                            grid = seq(-12, 12, length.out = 241))
  )
}

# Simulated draws from the corrected null. The exact integrals above and this
# simulation agree to three decimals; the simulation is kept because Fig. 1c
# needs the direction share conditional on a joint miss, which is read off the
# same simulated pair. (Fig. 2a's per-predictor benchmarks use null_exact.)
set.seed(2)
N  <- 2e6
u  <- rnorm(N); bg <- rnorm(N); bh <- rnorm(N)
calib <- function(b, target) {
  uniroot(function(s) mean(abs(s * b - u) > 1.96) - target, c(0.05, 20))$root
}
zgs <- calib(bg, mean(miss_g)) * bg - u
zhs <- calib(bh, mean(miss_h)) * bh - u
joint_s <- abs(zgs) > 1.96 & abs(zhs) > 1.96


# --- Shared style -----------------------------------------------------------

INK <- "#111111"; GREY <- "#5b5a57"; GRID <- "#e3e2dc"
RED <- "#c0392b"; BLUE <- "#1f6eb2"; PALE <- "#b9c9d9"
MUT <- "#c9c8c2"; GREEN <- "#1b8f63"; SAND <- "#d9a441"; FAINT <- "#dddcd6"

thm <- function(grid_y = TRUE, grid_x = FALSE) {
  theme_minimal(base_size = 9) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = if (grid_y) element_line(colour = GRID, linewidth = .3) else element_blank(),
      panel.grid.major.x = if (grid_x) element_line(colour = GRID, linewidth = .3) else element_blank(),
      axis.text   = element_text(colour = GREY, size = 7.5),
      axis.title  = element_text(colour = INK, size = 8),
      plot.title  = element_text(colour = INK, size = 8.4),
      plot.tag    = element_text(face = "bold", size = 11),
      plot.margin = margin(4, 6, 4, 4)
    )
}

err_w <- .14   # error-bar cap width on categorical axes


# --- Figure 1 ---------------------------------------------------------------
# a  the error plane, with the fitted error-on-error slope and its band
# b  same-side rate against both nulls, each with its interval
# c  direction given a joint miss, with its interval

n_hi <- sum(wide$zg >  1.96 & wide$zh >  1.96)
n_lo <- sum(wide$zg < -1.96 & wide$zh < -1.96)
n_op <- sum(joint & !same)

band <- inf$sb_primary$band %>% filter(x >= -12, x <= 12) %>%
  mutate(fit = pmax(pmin(fit, 12), -12),
         lo  = pmax(pmin(lo, 12), -12),
         hi  = pmax(pmin(hi, 12), -12))
message(sprintf("Fig. 1a slope %.2f (95%%%% CI %.2f-%.2f); quote it in the caption/text.",
                inf$sb_primary$slope,
                inf$sb_primary$slope - 1.96 * inf$sb_primary$slope_se,
                inf$sb_primary$slope + 1.96 * inf$sb_primary$slope_se))

pa <- wide %>%
  mutate(x   = pmax(pmin(zg, 12), -12),
         y   = pmax(pmin(zh, 12), -12),
         cls = case_when(same ~ "same", joint ~ "opp", TRUE ~ "rest")) %>%
  arrange(match(cls, c("rest", "opp", "same"))) %>%
  ggplot(aes(x, y)) +
  annotate("rect", xmin = 1.96, xmax = Inf, ymin = 1.96, ymax = Inf, fill = RED, alpha = .06) +
  annotate("rect", xmin = -Inf, xmax = -1.96, ymin = -Inf, ymax = -1.96, fill = RED, alpha = .06) +
  geom_hline(yintercept = 0, colour = "#9a9992", linewidth = .3) +
  geom_vline(xintercept = 0, colour = "#9a9992", linewidth = .3) +
  geom_hline(yintercept = c(-1.96, 1.96), colour = GREY, linetype = "22", linewidth = .3) +
  geom_vline(xintercept = c(-1.96, 1.96), colour = GREY, linetype = "22", linewidth = .3) +
  geom_point(aes(colour = cls, size = cls), alpha = .7, stroke = 0) +
  geom_ribbon(data = band, aes(x = x, ymin = lo, ymax = hi),
              inherit.aes = FALSE, fill = INK, alpha = .10) +
  geom_line(data = band, aes(x = x, y = fit),
            inherit.aes = FALSE, colour = INK, linewidth = .45) +
  scale_colour_manual(values = c(same = RED, opp = "#4c7ba6", rest = MUT), guide = "none") +
  scale_size_manual(values = c(same = .9, opp = .9, rest = .55), guide = "none") +
  annotate("text", x = 12.6, y = 13.1, label = sprintf("both overshoot\n%d", n_hi),
           hjust = 1, vjust = 1, size = 2.6, colour = RED, lineheight = .95) +
  annotate("text", x = -12.6, y = -13.1, label = sprintf("both undershoot\n%d", n_lo),
           hjust = 0, vjust = 0, size = 2.6, colour = RED, lineheight = .95) +
  annotate("text", x = -12.6, y = 13.1, label = sprintf("opposite directions\n%d", n_op),
           hjust = 0, vjust = 1, size = 2.5, colour = "#4c7ba6", lineheight = .95) +
  scale_x_continuous(limits = c(-12.8, 12.8), breaks = seq(-10, 10, 5)) +
  scale_y_continuous(limits = c(-13.3, 13.3), breaks = seq(-10, 10, 5)) +
  labs(x = "GPT-4 forecasting error\n(standard errors of the observed effect)",
       y = "Human forecasting error\n(standard errors of the observed effect)",
       title = "The two predictors' errors line up, comparison by comparison") +
  thm(grid_y = FALSE)

# Panel b: observed bar carries the study-bootstrap interval; the corrected
# null carries its finite-sample Monte-Carlo interval; the naive null a
# binomial interval at the same archive size.
naive_same <- mean(miss_g) * mean(miss_h) / 2
naive_ci   <- qbinom(c(.025, .975), nrow(wide), naive_same) / nrow(wide)

b1 <- data.frame(
  x = factor(c("Observed", "Corrected null", "Null ignoring\nshared estimate"),
             levels = c("Observed", "Corrected null", "Null ignoring\nshared estimate")),
  v  = 100 * c(mean(same), nulls[["same"]], naive_same),
  lo = 100 * c(inf$boot_ci[1, "same_side"], inf$null_ci[1, "same_side"], naive_ci[1]),
  hi = 100 * c(inf$boot_ci[2, "same_side"], inf$null_ci[2, "same_side"], naive_ci[2]),
  f  = c(RED, MUT, FAINT)
)

pb <- ggplot(b1, aes(x, v)) +
  geom_col(fill = b1$f, width = .55) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = err_w,
                colour = INK, linewidth = .35) +
  geom_text(aes(y = hi, label = sprintf("%.1f%%", v)), vjust = -.55, size = 2.9,
            fontface = "bold", colour = INK) +
  scale_y_continuous(limits = c(0, max(b1$hi) + 3), breaks = c(0, 10, 20),
                     labels = c("0", "10%", "20%")) +
  labs(x = NULL, y = NULL,
       title = sprintf("Share of all %s comparisons\nthat both predictors miss the same way",
                       format(nrow(wide), big.mark = ","))) +
  thm() +
  theme(axis.text.x = element_text(colour = INK, size = 6.9, lineheight = .95))

dirv     <- 100 * sum(same) / sum(joint)
dir_null <- 100 * mean(sign(zgs[joint_s]) == sign(zhs[joint_s]))

c1 <- data.frame(
  x = factor(c("Same\ndirection", "Opposite\ndirections"),
             levels = c("Same\ndirection", "Opposite\ndirections")),
  v  = c(dirv, 100 - dirv),
  lo = c(100 * inf$boot_ci[1, "direction"], 100 - 100 * inf$boot_ci[2, "direction"]),
  hi = c(100 * inf$boot_ci[2, "direction"], 100 - 100 * inf$boot_ci[1, "direction"]),
  f  = c(RED, PALE)
)

pc <- ggplot(c1, aes(x, v)) +
  geom_col(fill = c1$f, width = .55) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = err_w,
                colour = INK, linewidth = .35) +
  geom_text(aes(y = hi, label = sprintf("%.0f%%", v)), vjust = -.55, size = 2.9,
            fontface = "bold", colour = INK) +
  geom_hline(yintercept = dir_null, colour = GREY, linetype = "22", linewidth = .4) +
  annotate("text", x = 2.55, y = dir_null + 6.5,
           label = sprintf("corrected null, %.0f%%", dir_null),
           hjust = 1, size = 2.4, colour = GREY) +
  geom_hline(yintercept = 50, colour = "#b9b8b2", linetype = "12", linewidth = .35) +
  annotate("text", x = 2.55, y = 41, label = "chance, 50%", hjust = 1,
           size = 2.4, colour = "#9a9992") +
  scale_y_continuous(limits = c(0, 112), breaks = c(0, 50, 100), labels = c("0", "50%", "100%")) +
  labs(x = NULL, y = NULL,
       title = sprintf("Of the %d comparisons both miss,\nthe share missed in the same direction",
                       sum(joint))) +
  thm() +
  theme(axis.text.x = element_text(colour = INK, size = 7.2, lineheight = .95))

fig1 <- pa + (pb / pc) + plot_layout(widths = c(1.52, 1)) +
  plot_annotation(tag_levels = "a") & theme(plot.tag.position = c(0, 1))

ggsave(file.path(fig_dir, "F1_same_effects.png"), fig1,
       width = 7.2, height = 4.3, dpi = 400, bg = "white")


# --- Figure 2 ---------------------------------------------------------------
# a  every other predictor scored on GPT-4's misses (lay crowd marked out)
# b  miss rates in the secondary archive: GPT-4, domain experts, and experts
#    conditional on a GPT-4 miss
# c  direction agreement among the conditions both miss (was a text tack-on)

zwide <- read.csv(file.path(derived_dir, "comparisons_all_models.csv"),
                  check.names = FALSE) %>%
  select(study, outcome.name, reference_condition, condition.name, model, z_recal) %>%
  pivot_wider(names_from = model, values_from = z_recal)

gp <- zwide$`gpt-4`

mods <- c(`deepseek/deepseek-chat-v3-0324` = "DeepSeek-V3",
          `openai/gpt-oss-120b`            = "GPT-OSS-120B",
          `gpt-3.5-turbo`                  = "GPT-3.5",
          `google/gemma-3-27b-it`          = "Gemma-3-27B",
          `babbage-002`                    = "babbage-002",
          `davinci-002`                    = "davinci-002",
          human                            = "Human crowd (lay)")

tab <- map_dfr(names(mods), function(mc) {
  zo  <- zwide[[mc]]
  ok  <- !is.na(zo) & !is.na(gp)
  gm  <- ok & abs(gp) > 1.96
  k   <- sum((abs(zo) > 1.96 & sign(zo) == sign(gp))[gm])
  ci  <- wilson_ci(k, sum(gm))
  gmr <- mean(abs(gp[ok]) > 1.96)                    # GPT-4's miss rate on this overlap
  omr <- mean(abs(zo[ok]) > 1.96)                    # this predictor's miss rate
  ne  <- null_exact(gmr, omr)                        # exact corrected null for this pair
  data.frame(name  = mods[[mc]],
             v     = 100 * k / sum(gm),
             lo    = 100 * ci[["lo"]],
             hi    = 100 * ci[["hi"]],
             r     = cor(gp[ok], zo[ok]),
             bench = 100 * ne[["same"]] / gmr)       # same-direction share | GPT-4 misses
}) %>%
  mutate(name = factor(name, levels = rev(name)),
         col  = ifelse(grepl("Human", name), GREEN, BLUE))

# Per-predictor corrected-null benchmarks differ slightly (the crowd's miss
# rate is not the models'), so the benchmark is drawn as a band spanning their
# range rather than a single line.
bench_lo <- min(tab$bench)
bench_hi <- max(tab$bench)

p2a <- ggplot(tab, aes(v, name)) +
  annotate("rect", xmin = bench_lo, xmax = bench_hi, ymin = -Inf, ymax = Inf,
           fill = GREY, alpha = .15) +
  geom_col(fill = tab$col, width = .62) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = err_w,
                 colour = INK, linewidth = .35) +
  geom_text(aes(x = hi, label = sprintf("%.0f%%", v)), hjust = -.22, size = 2.9,
            fontface = "bold", colour = INK) +
  geom_text(aes(x = 118, label = sprintf("r = %.2f", r)), hjust = 1, size = 2.5,
            fontface = "italic", colour = GREY) +
  annotate("text", x = (bench_lo + bench_hi) / 2, y = 8.9,
           label = sprintf("%.0f–%.0f%% expected under the\ncorrected null (range across predictors)",
                           bench_lo, bench_hi),
           hjust = .5, vjust = 1, size = 2.35, colour = GREY, lineheight = .95) +
  scale_y_discrete(expand = expansion(add = c(.55, 1.7))) +
  scale_x_continuous(limits = c(0, 120), breaks = c(0, 25, 50, 75),
                     labels = c("0", "25%", "50%", "75%")) +
  labs(x = "Share of GPT-4's misses also missed\nin the same direction",
       y = NULL,
       title = "Other predictors err in the same direction as GPT-4") +
  coord_cartesian(clip = "off") +
  thm(grid_y = FALSE, grid_x = TRUE) +
  theme(axis.text.y = element_text(colour = INK, size = 7.6))

a2_path <- file.path(derived_dir, "comparisons_archive2.csv")

if (file.exists(a2_path)) {
  a2 <- read.csv(a2_path, check.names = FALSE)
  a2 <- a2[!is.na(a2$z_recal_gpt) & !is.na(a2$z_recal_expert), ]

  gm2   <- abs(a2$z_recal_gpt)    > 1.96
  em2   <- abs(a2$z_recal_expert) > 1.96
  both2 <- gm2 & em2
  same2 <- both2 & (sign(a2$z_recal_gpt) == sign(a2$z_recal_expert))

  ci_g    <- wilson_ci(sum(gm2), length(gm2))
  ci_e    <- wilson_ci(sum(em2), length(em2))
  ci_also <- wilson_ci(sum(both2), sum(gm2))
  ci_dir  <- wilson_ci(sum(same2), sum(both2))
  n2      <- null_exact(mean(gm2), mean(em2))

  b2 <- data.frame(
    x = factor(c("GPT-4\nmisses", "Experts\nmiss", "Experts miss,\namong GPT-4's\nmisses"),
               levels = c("GPT-4\nmisses", "Experts\nmiss", "Experts miss,\namong GPT-4's\nmisses")),
    v  = 100 * c(mean(gm2), mean(em2), sum(both2) / sum(gm2)),
    lo = 100 * c(ci_g[["lo"]], ci_e[["lo"]], ci_also[["lo"]]),
    hi = 100 * c(ci_g[["hi"]], ci_e[["hi"]], ci_also[["hi"]]),
    f  = c(BLUE, SAND, RED)
  )

  p2b <- ggplot(b2, aes(x, v)) +
    geom_col(fill = b2$f, width = .55) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = err_w,
                  colour = INK, linewidth = .35) +
    geom_text(aes(y = hi, label = sprintf("%.0f%%", v)), vjust = -.55, size = 2.9,
              fontface = "bold", colour = INK) +
    scale_y_continuous(limits = c(0, 118), breaks = c(0, 50, 100),
                       labels = c("0", "50%", "100%")) +
    labs(x = NULL, y = NULL,
         title = sprintf("Domain experts, secondary archive\n(%d conditions)", nrow(a2))) +
    thm() +
    theme(axis.text.x = element_text(colour = INK, size = 6.8, lineheight = .95))

  # Panel c (C24): the direction result as its own panel, mirroring Fig. 1c
  # so the two figures report the same quantities the same way.
  dir2      <- 100 * sum(same2) / sum(both2)
  dir_null2 <- 100 * n2[["direction"]]

  c2 <- data.frame(
    x = factor(c("Same\ndirection", "Opposite\ndirections"),
               levels = c("Same\ndirection", "Opposite\ndirections")),
    v  = c(dir2, 100 - dir2),
    lo = c(100 * ci_dir[["lo"]], 100 - 100 * ci_dir[["hi"]]),
    hi = c(100 * ci_dir[["hi"]], 100 - 100 * ci_dir[["lo"]]),
    f  = c(RED, PALE)
  )

  p2c <- ggplot(c2, aes(x, v)) +
    geom_col(fill = c2$f, width = .55) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = err_w,
                  colour = INK, linewidth = .35) +
    geom_text(aes(y = hi, label = sprintf("%.0f%%", v)), vjust = -.55, size = 2.9,
              fontface = "bold", colour = INK) +
    geom_hline(yintercept = dir_null2, colour = GREY, linetype = "22", linewidth = .4) +
    annotate("text", x = 2.55, y = dir_null2 - 7,
             label = sprintf("corrected null, %.0f%%", dir_null2),
             hjust = 1, size = 2.4, colour = GREY) +
    geom_hline(yintercept = 50, colour = "#b9b8b2", linetype = "12", linewidth = .35) +
    annotate("text", x = 2.55, y = 43, label = "chance, 50%", hjust = 1,
             size = 2.4, colour = "#9a9992") +
    scale_y_continuous(limits = c(0, 118), breaks = c(0, 50, 100),
                       labels = c("0", "50%", "100%")) +
    labs(x = NULL, y = NULL,
         title = sprintf("Of the %d conditions both miss,\nshare in the same direction", sum(both2))) +
    thm() +
    theme(axis.text.x = element_text(colour = INK, size = 7.2, lineheight = .95))

  fig2 <- p2a + p2b + p2c + plot_layout(widths = c(1.5, 1, 1)) +
    plot_annotation(tag_levels = "a") & theme(plot.tag.position = c(0, 1))

  ggsave(file.path(fig_dir, "F2_shared_blind_spot.png"), fig2,
         width = 9.6, height = 2.85, dpi = 400, bg = "white")
} else {
  warning("comparisons_archive2.csv not found; Fig. 2 panels b-c skipped. ",
          "Run 02_build_archive2.R first.")
}

message("05_figures.R complete. Figures written to output/figures/")
