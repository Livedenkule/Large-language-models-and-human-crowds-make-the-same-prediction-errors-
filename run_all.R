# ---------------------------------------------------------------------------
# run_all.R
#
# Reproduces every quantity, table and figure in
#
#   Wilhelmsen, L. L., Esfandiari, R. & Gollwitzer, A.
#   Large language models and human forecasters make the same prediction errors.
#   Matters Arising on Ashokkumar, Hewitt, Ghezae & Willer, Nature (2026).
#
# Usage
#   From the project root:      Rscript run_all.R
#   Or in RStudio:              open shared-error-reanalysis.Rproj, then
#                               source("run_all.R")
#
#   Skip the two build scripts and run from the deposited derived data:
#                               Rscript run_all.R --from-derived
#
# Runtime is roughly 10-20 minutes end to end, dominated by the study-specific
# null and the bootstrap in 04_robustness.R.
# ---------------------------------------------------------------------------

args         <- commandArgs(trailingOnly = TRUE)
from_derived <- "--from-derived" %in% args

if (!requireNamespace("here", quietly = TRUE)) {
  stop("Package 'here' is required. Install it with install.packages('here').")
}
library(here)

# --- Dependency check -------------------------------------------------------

required <- c("dplyr", "tidyr", "purrr", "stringr", "stringi", "forcats",
              "broom", "ggplot2", "patchwork", "here")
absent   <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]

if (length(absent)) {
  stop("Missing packages: ", paste(absent, collapse = ", "), "\n",
       "Install them with:\n  install.packages(c(",
       paste0('"', absent, '"', collapse = ", "), "))\n",
       "Or restore the recorded versions with:  renv::restore()")
}


# --- Pipeline ---------------------------------------------------------------

build_scripts <- c("01_build_comparisons.R", "02_build_archive2.R")
analysis_scripts <- c("03_main_results.R", "04_robustness.R", "05_figures.R")

scripts <- if (from_derived) analysis_scripts else c(build_scripts, analysis_scripts)

if (from_derived) {
  message("Running from the deposited derived data; skipping the build scripts.\n")
} else if (!file.exists(here("data", "raw", "rct_responses.RDS"))) {
  stop("data/raw/ is empty. Either place the Code Ocean capsule files there ",
       "(see data/raw/README.md) or run:\n  Rscript run_all.R --from-derived")
}

source(here("R", "00_functions.R"))

for (s in scripts) {
  message("\n", strrep("=", 70), "\n  ", s, "\n", strrep("=", 70))
  t0 <- Sys.time()
  source(here("R", s), echo = FALSE)
  message("  done in ", round(difftime(Sys.time(), t0, units = "secs")), "s")
}


# --- Provenance -------------------------------------------------------------

dir.create(here("output", "tables"), showWarnings = FALSE, recursive = TRUE)
write_reported_values(here("output", "tables", "reported_values.csv"))

writeLines(capture.output(sessionInfo()), here("output", "sessionInfo.txt"))

message("\n", strrep("=", 70),
        "\n  Complete. Tables in output/tables/, figures in output/figures/.",
        "\n  Every quantity quoted in the manuscript is listed in",
        "\n  output/tables/reported_values.csv.\n", strrep("=", 70))
