# ---------------------------------------------------------------------------
# 05_figures.R
#
# Figures 1 and 2 of the main text, drawn from the derived comparison-level
# data so that the numbers in the figures and the numbers in the text cannot
# drift apart.
#
# INPUT   data/derived/comparisons_primary.csv
#         data/derived/comparisons_all_models.csv
#         data/derived/comparisons_archive2.csv
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

# Simulated draws from the corrected null, used for the benchmark lines. The
# exact integrals above and this simulation agree to three decimals; the
# simulation is kept because Fig. 2a needs the conditional share, which is read
# off the same simulated pair.
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


# --- Figure 1 ---------------------------------------------------------------
# a  the error plane
# b  same-side rate against both nulls
# c  direction given a joint miss

n_hi <- sum(wide$zg >  1.96 & wide$zh >  1.96)
n_lo <- sum(wide$zg < -1.96 & wide$zh < -1.96)
n_op <- sum(joint & !same)

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
  labs(x = "GPT-4 error  (standard errors of the observed effect)",
       y = "Human crowd error  (standard errors)",
       title = "The two predictors' errors line up, comparison by comparison") +
  thm(grid_y = FALSE)

b1 <- data.frame(
  x = factor(c("Observed", "Corrected null", "Null ignoring\nshared estimate"),
             levels = c("Observed", "Corrected null", "Null ignoring\nshared estimate")),
  v = 100 * c(mean(same), nulls[["same"]], mean(miss_g) * mean(miss_h) / 2),
  f = c(RED, MUT, FAINT)
)

pb <- ggplot(b1, aes(x, v)) +
  geom_col(fill = b1$f, width = .55) +
  geom_text(aes(label = sprintf("%.1f%%", v)), vjust = -.55, size = 2.9,
            fontface = "bold", colour = INK) +
  annotate("curve", x = 1.2, xend = 1.8,
           y = 100 * mean(same) + 1.5, yend = 100 * nulls[["same"]] + 2.5,
           curvature = -.35, colour = GREY, linewidth = .3) +
  annotate("text", x = 1.5, y = 25.4,
           label = sprintf("%.1f×", mean(same) / nulls[["same"]]),
           size = 2.9, fontface = "bold", colour = GREY) +
  scale_y_continuous(limits = c(0, 27), breaks = c(0, 10, 20), labels = c("0", "10%", "20%")) +
  labs(x = NULL, y = NULL,
       title = "Share of all 1,678 comparisons\nthat both predictors miss the same way") +
  thm() +
  theme(axis.text.x = element_text(colour = INK, size = 6.9, lineheight = .95))

dirv     <- 100 * sum(same) / sum(joint)
dir_null <- 100 * mean(sign(zgs[joint_s]) == sign(zhs[joint_s]))

c1 <- data.frame(
  x = factor(c("Same\ndirection", "Opposite\ndirections"),
             levels = c("Same\ndirection", "Opposite\ndirections")),
  v = c(dirv, 100 - dirv),
  f = c(RED, PALE)
)

pc <- ggplot(c1, aes(x, v)) +
  geom_col(fill = c1$f, width = .55) +
  geom_text(aes(label = sprintf("%.0f%%", v)), vjust = -.55, size = 2.9,
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
       title = "Of the 339 comparisons both miss,\nhow the two misses line up") +
  thm() +
  theme(axis.text.x = element_text(colour = INK, size = 7.2, lineheight = .95))

fig1 <- pa + (pb / pc) + plot_layout(widths = c(1.52, 1)) +
  plot_annotation(tag_levels = "a") & theme(plot.tag.position = c(0, 1))

ggsave(file.path(fig_dir, "F1_same_effects.png"), fig1,
       width = 7.2, height = 4.15, dpi = 400, bg = "white")


# --- Figure 2 ---------------------------------------------------------------
# a  every other predictor against GPT-4's misses
# b  domain experts in the secondary archive

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
          human                            = "Human crowd")

tab <- map_dfr(names(mods), function(mc) {
  zo <- zwide[[mc]]
  ok <- !is.na(zo) & !is.na(gp)
  gm <- ok & abs(gp) > 1.96
  data.frame(name = mods[[mc]],
             v = 100 * mean((abs(zo) > 1.96 & sign(zo) == sign(gp))[gm]),
             r = cor(gp[ok], zo[ok]))
}) %>%
  mutate(name = factor(name, levels = rev(name)),
         col  = ifelse(name == "Human crowd", GREEN, BLUE))

bench <- 100 * mean((abs(zhs) > 1.96 & sign(zhs) == sign(zgs))[abs(zgs) > 1.96])

p2a <- ggplot(tab, aes(v, name)) +
  geom_col(fill = tab$col, width = .62) +
  geom_text(aes(label = sprintf("%.0f%%", v)), hjust = -.25, size = 2.9,
            fontface = "bold", colour = INK) +
  geom_text(aes(x = 118, label = sprintf("r = %.2f", r)), hjust = 1, size = 2.5,
            fontface = "italic", colour = GREY) +
  geom_vline(xintercept = bench, colour = GREY, linetype = "22", linewidth = .4) +
  annotate("text", x = bench + 2, y = 8.9,
           label = sprintf("%.0f%% expected under the corrected null", bench),
           hjust = 0, vjust = 1, size = 2.35, colour = GREY) +
  scale_y_discrete(expand = expansion(add = c(.55, 1.7))) +
  scale_x_continuous(limits = c(0, 120), breaks = c(0, 25, 50, 75),
                     labels = c("0", "25%", "50%", "75%")) +
  labs(x = "Also wrong, on the same side as GPT-4", y = NULL,
       title = "Every other predictor inherits GPT-4’s misses") +
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

  b2 <- data.frame(
    x = factor(c("GPT-4\nmisses", "Experts\nmiss", "Experts also\nmiss, given\nGPT-4 misses"),
               levels = c("GPT-4\nmisses", "Experts\nmiss", "Experts also\nmiss, given\nGPT-4 misses")),
    v = 100 * c(mean(gm2), mean(em2), sum(both2) / sum(gm2)),
    f = c(BLUE, SAND, RED)
  )

  p2b <- ggplot(b2, aes(x, v)) +
    geom_col(fill = b2$f, width = .55) +
    geom_text(aes(label = sprintf("%.0f%%", v)), vjust = -.55, size = 2.9,
              fontface = "bold", colour = INK) +
    annotate("text", x = 2.42, y = 112,
             label = sprintf("%.0f%% of those on\nthe same side",
                             100 * sum(same2) / sum(both2)),
             hjust = 1, size = 2.4, colour = RED, lineheight = .95) +
    scale_y_continuous(limits = c(0, 128), breaks = c(0, 50, 100),
                       labels = c("0", "50%", "100%")) +
    labs(x = NULL, y = NULL, title = "Domain experts are not an\nindependent check") +
    thm() +
    theme(axis.text.x = element_text(colour = INK, size = 6.8, lineheight = .95))

  fig2 <- p2a + p2b + plot_layout(widths = c(1.35, 1)) +
    plot_annotation(tag_levels = "a") & theme(plot.tag.position = c(0, 1))

  ggsave(file.path(fig_dir, "F2_shared_blind_spot.png"), fig2,
         width = 7.2, height = 2.7, dpi = 400, bg = "white")
} else {
  warning("comparisons_archive2.csv not found; Fig. 2 panel b skipped. ",
          "Run 02_build_archive2.R first.")
}

message("05_figures.R complete. Figures written to output/figures/")
