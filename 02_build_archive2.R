# ---------------------------------------------------------------------------
# 02_build_archive2.R
#
# Builds the condition-level dataset for the secondary archive of megastudies,
# in which the original teams collected forecasts from domain experts.
#
# INPUT   data/raw/megastudies.RDS
# OUTPUT  data/derived/comparisons_archive2.csv
#
# ---------------------------------------------------------------------------
# ACTION REQUIRED BEFORE DEPOSIT
#
# The derivation of comparisons_archive2.csv lives in a separate script that
# was not part of shared_error_analysis_5.Rmd (the Rmd reads the finished CSV).
# Paste that script's body into the section marked below, so that the deposited
# repository rebuilds the file rather than shipping it as an unexplained input.
#
# Per Supplementary Results 2.6, that script must:
#
#   1. Read the authors' megastudies.RDS.
#   2. Apply their scale conversions to a common proportion-of-scale metric.
#   3. Apply their per-family scale ratios.
#   4. Apply their two dataset corrections.
#   5. Recalibrate within megastudy family, leaving out one condition at a time,
#      since outcome scales differ across families.
#   6. Emit one row per condition with, at minimum:
#         family              megastudy family identifier
#         condition           condition identifier
#         estimate, std.error published effect and its standard error
#         z_recal_gpt         GPT-4's recalibrated standardized error
#         z_recal_expert      the experts' recalibrated standardized error
#
# The checks at the foot of this script verify those expectations against the
# figures reported in the manuscript, so a correct paste will pass them and an
# incorrect one will fail loudly.
# ---------------------------------------------------------------------------

library(dplyr)
library(here)

source(here("R", "00_functions.R"))

raw_dir <- here("data", "raw")
out_dir <- here("data", "derived")

archive2_path <- file.path(out_dir, "comparisons_archive2.csv")


# --- Derivation -------------------------------------------------------------
# >>> PASTE THE ARCHIVE-2 BUILD SCRIPT HERE <<<
#
# It should end by writing archive2_path.

if (!file.exists(archive2_path)) {
  stop(
    "comparisons_archive2.csv not found and no derivation supplied.\n",
    "Paste the archive-2 build script into 02_build_archive2.R before deposit ",
    "(see the header of this file for what it must produce)."
  )
}


# --- Checks -----------------------------------------------------------------
# These reproduce the figures quoted in Supplementary Results 2.6 and so verify
# that the derivation above produced the intended file.

a2 <- read.csv(archive2_path, check.names = FALSE)

required <- c("family", "estimate", "z_recal_gpt", "z_recal_expert")
missing_cols <- setdiff(required, names(a2))
if (length(missing_cols)) {
  stop("comparisons_archive2.csv is missing required columns: ",
       paste(missing_cols, collapse = ", "))
}

message("Secondary archive: ", nrow(a2), " effects across ",
        n_distinct(a2$family), " megastudies.")

report_value("SI 2.6", "effects in secondary archive", nrow(a2), 0)
report_value("SI 2.6", "megastudies", n_distinct(a2$family), 0)

with_expert <- a2[!is.na(a2$z_recal_gpt) & !is.na(a2$z_recal_expert), ]
message("Conditions with expert forecasts: ", nrow(with_expert))
report_value("SI 2.6", "conditions with expert forecasts", nrow(with_expert), 0)

if (nrow(with_expert) != 376L) {
  warning("Expected 376 conditions with expert forecasts, found ", nrow(with_expert),
          ". Check the scale conversions and dataset corrections.")
}
