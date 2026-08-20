# ---------------------------------------------------------------------------
# tests/test_functions.R
#
# Checks on the helper functions in R/00_functions.R. Requires no data: run it
# to confirm the analytical machinery behaves before committing to a full run.
#
#   Rscript tests/test_functions.R
#
# The first block is the important one. The corrected null is the least
# transparent step in the analysis, so it is verified against every benchmark
# quoted in the Supplementary Information, and against the theoretical identity
# E[z_g z_h] = Var(u) = 1 that holds exactly under that null.
# ---------------------------------------------------------------------------

if (requireNamespace("here", quietly = TRUE)) {
  source(here::here("R", "00_functions.R"))
} else if (file.exists("R/00_functions.R")) {
  source("R/00_functions.R")                 # run from the project root
} else {
  source("../R/00_functions.R")              # run from tests/
}

failures <- 0L

check <- function(label, ok) {
  failures <<- failures + as.integer(!isTRUE(ok))
  cat(sprintf("%-62s %s\n", label, if (isTRUE(ok)) "pass" else "FAIL"))
}

near <- function(a, b, tol = 0.002) all(abs(a - b) < tol)


# --- 1. The corrected null reproduces the published benchmarks --------------
# Rows of Supplementary Table 1: marginal miss rates, threshold, and the joint,
# same-side and direction nulls reported in parentheses there.

cat("\nCorrected null against Supplementary Table 1\n")

cases <- list(
  list(name = "No recalibration, 1.96",       g = 0.463, h = 0.498, t = 1.96,
       want = c(0.232, 0.141, 0.608)),
  list(name = "Linear recalibration, 1.96",   g = 0.357, h = 0.380, t = 1.96,
       want = c(0.141, 0.100, 0.712)),
  list(name = "Isotonic recalibration, 1.96", g = 0.381, h = 0.352, t = 1.96,
       want = c(0.140, 0.100, 0.714)),
  list(name = "Linear, threshold 1.64",       g = 0.413, h = 0.443, t = 1.64,
       want = c(0.189, 0.134, 0.710)),
  list(name = "Linear, threshold 2.58",       g = 0.246, h = 0.274, t = 2.58,
       want = c(0.071, 0.053, 0.739)),
  list(name = "Linear, threshold 3.29",       g = 0.169, h = 0.186, t = 3.29,
       want = c(0.034, 0.026, 0.755))
)

for (cs in cases) {
  got <- null_exact(cs$g, cs$h, t = cs$t)
  check(cs$name, near(unname(got), cs$want))
}


# --- 2. Heavier-tailed variants --------------------------------------------
# SI 1.4 reports a same-side benchmark of 9.6% and a direction benchmark of
# 67.6% under a Student t(3) error shape.

cat("\nTail variants\n")
t3 <- null_exact_tail(0.357, 0.380, ptail = t3_tail)
check("t(3) same-side benchmark = 9.6%",  near(t3[["same"]], 9.6, tol = 0.15))
check("t(3) direction benchmark = 67.6%", near(t3[["direction"]], 67.6, tol = 0.5))
check("heavier tails lower the same-side benchmark", t3[["same"]] < 10.0)


# --- 3. The theoretical identity -------------------------------------------
# Under the corrected null the shared term contributes covariance exactly 1.

cat("\nTheoretical identity\n")
set.seed(2)
N <- 5e5
u <- rnorm(N); bg <- rnorm(N); bh <- rnorm(N)
cal <- function(b, target) uniroot(function(s) mean(abs(s * b - u) > 1.96) - target,
                                   c(0.05, 20))$root
zgs <- cal(bg, 0.357) * bg - u
zhs <- cal(bh, 0.380) * bh - u
check("simulated E[zg*zh] = 1", near(mean(zgs * zhs), 1, tol = 0.02))

sim <- c(joint = mean(abs(zgs) > 1.96 & abs(zhs) > 1.96),
         same  = mean(abs(zgs) > 1.96 & abs(zhs) > 1.96 & sign(zgs) == sign(zhs)))
exact <- null_exact(0.357, 0.380)
check("simulation agrees with exact integration",
      near(sim[["joint"]], exact[["joint"]], tol = 0.003) &&
      near(sim[["same"]],  exact[["same"]],  tol = 0.003))


# --- 4. Recalibration helpers ----------------------------------------------

cat("\nRecalibration\n")
set.seed(99)
n <- 400
study <- rep(paste0("s", 1:20), each = 20)
pred  <- rnorm(n)
est   <- 0.5 * pred + rnorm(n, 0, 0.2)

# linear_loo must equal the explicit leave-one-study-out loop.
manual <- rep(NA_real_, n)
for (s in unique(study)) {
  fit <- lm(est ~ pred, data = data.frame(est = est, pred = pred)[study != s, ])
  manual[study == s] <- predict(fit, newdata = data.frame(pred = pred[study == s]))
}
check("linear_loo matches the explicit loop", near(manual, linear_loo(pred, est, study), 1e-10))

# No study may contribute to its own correction: perturbing one study's outcomes
# must leave that study's own recalibrated predictions unchanged.
est2 <- est
est2[study == "s1"] <- est2[study == "s1"] + 5
check("held-out study does not correct itself",
      near(linear_loo(pred, est, study)[study == "s1"],
           linear_loo(pred, est2, study)[study == "s1"], 1e-10))

iso <- iso_loo(pred, est, study)
check("iso_loo returns no missing values", !any(is.na(iso)))
check("iso_loo is monotone within each study",
      all(vapply(unique(study), function(s) {
        i <- study == s; o <- order(pred[i])
        all(diff(iso[i][o]) >= -1e-9)
      }, logical(1))))


# --- 5. Cluster-robust standard errors --------------------------------------

cat("\nClustering\n")
zg <- rnorm(n); zh <- 0.5 * zg + rnorm(n)
fit <- lm(zh ~ zg)
X <- model.matrix(fit); u2 <- residuals(fit); XtXi <- solve(crossprod(X))
meat <- Reduce(`+`, lapply(split(seq_len(nrow(X)), study), function(i) {
  v <- crossprod(X[i, , drop = FALSE], u2[i]); tcrossprod(v)
}))
manual_se <- sqrt(diag(XtXi %*% meat %*% XtXi))[2]
check("cluster_se matches the explicit sandwich",
      near(unname(manual_se), unname(cluster_se(fit, study)), 1e-10))


# --- 6. Summaries -----------------------------------------------------------

cat("\nSummaries\n")
s1 <- summarise_shared(zg, zh)
check("summarise_shared rates lie in [0, 1]",
      all(unlist(s1[c("miss_g", "miss_h", "joint", "same_side", "direction")]) >= 0) &&
      all(unlist(s1[c("miss_g", "miss_h", "joint", "same_side", "direction")]) <= 1))
check("same-side misses are a subset of joint misses", s1$same_side <= s1$joint)
check("raising the threshold lowers the miss rate",
      summarise_shared(zg, zh, 2.58)$joint <= s1$joint)


# --- Result -----------------------------------------------------------------

cat("\n", strrep("-", 70), "\n", sep = "")
if (failures == 0L) {
  cat("All checks passed.\n")
} else {
  cat(failures, "check(s) FAILED.\n")
  quit(status = 1L)
}
