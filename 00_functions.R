# ---------------------------------------------------------------------------
# 00_functions.R
#
# Helper functions used by more than one analysis script. Sourced by every
# other script in R/ via run_all.R; sourcing it directly has no side effects
# beyond defining functions.
#
# Wilhelmsen, Esfandiari & Gollwitzer
# Matters Arising on Ashokkumar, Hewitt, Ghezae & Willer (Nature, 2026)
# ---------------------------------------------------------------------------


# --- Corrected null ---------------------------------------------------------

#' Exact benchmarks for two predictors scored against the same estimate
#'
#' Both predictors are graded against the same estimated effect, whose sampling
#' noise enters every standardized error with the same sign. Writing
#' z_i = a_i - u, with a_i the predictor's own deviation from the true effect
#' and u = eps / s.e. the estimate's standardized noise (variance 1 by
#' construction), independent predictors still show correlated errors, joint
#' misses above the product of the marginal rates, and directional agreement
#' above 50%. This function returns those benchmarks exactly.
#'
#' Under this null the error pair is bivariate normal, so each benchmark is a
#' one-dimensional integral over the shared term u and needs no simulation.
#'
#' @param miss_g,miss_h Marginal miss rates of the two predictors.
#' @param t Threshold on |z| defining a miss (1.96 for a 95% interval).
#' @return Named vector: joint-miss rate, same-side rate, direction share.
null_exact <- function(miss_g, miss_h, t = 1.96) {
  sig <- function(p) t / stats::qnorm(1 - p / 2)   # total sd yielding this miss rate
  s2g <- max(sig(miss_g)^2 - 1, 0)
  s2h <- max(sig(miss_h)^2 - 1, 0)

  tail_hi <- function(s2, u) {
    if (s2 > 0) stats::pnorm((t + u) / sqrt(s2), lower.tail = FALSE) else as.numeric(-u > t)
  }
  tail_lo <- function(s2, u) {
    if (s2 > 0) stats::pnorm((-t + u) / sqrt(s2)) else as.numeric(-u < -t)
  }

  f <- function(u, kind) vapply(u, function(uu) {
    gh <- tail_hi(s2g, uu); gl <- tail_lo(s2g, uu)
    hh <- tail_hi(s2h, uu); hl <- tail_lo(s2h, uu)
    v <- if (kind == "joint") (gh + gl) * (hh + hl) else gh * hh + gl * hl
    v * stats::dnorm(uu)
  }, numeric(1))

  joint <- stats::integrate(f, -9, 9, kind = "joint",
                            rel.tol = 1e-9, subdivisions = 500L)$value
  same  <- stats::integrate(f, -9, 9, kind = "same",
                            rel.tol = 1e-9, subdivisions = 500L)$value
  c(joint = joint, same = same, direction = same / joint)
}


#' As null_exact(), but for a non-normal error shape
#'
#' The predictors' real errors have heavier tails than a normal. This variant
#' takes an arbitrary standardized error shape through its upper-tail function
#' and recalibrates each predictor's scale to reproduce its observed miss rate.
#' Heavier tails lower the benchmarks slightly, so the normal case reported in
#' the manuscript is the conservative one.
#'
#' @param miss_g,miss_h Marginal miss rates of the two predictors.
#' @param t Threshold on |z| defining a miss.
#' @param ptail Function giving P(b > q) for the standardized error shape.
#' @return Named vector: same-side rate and direction share, both as percentages.
null_exact_tail <- function(miss_g, miss_h, t = 1.96, ptail) {
  miss_of <- function(s) {
    stats::integrate(function(u) vapply(u, function(uu) {
      (ptail((t + uu) / s) + ptail((t - uu) / s)) * stats::dnorm(uu)
    }, numeric(1)), -9, 9, subdivisions = 500L, stop.on.error = FALSE)$value
  }
  cal <- function(p) stats::uniroot(function(s) miss_of(s) - p, c(0.05, 30))$root
  sg <- cal(miss_g); sh <- cal(miss_h)

  f <- function(u, kind) vapply(u, function(uu) {
    gh <- ptail((t + uu) / sg); gl <- ptail((t - uu) / sg)
    hh <- ptail((t + uu) / sh); hl <- ptail((t - uu) / sh)
    v <- if (kind == "joint") (gh + gl) * (hh + hl) else gh * hh + gl * hl
    v * stats::dnorm(uu)
  }, numeric(1))

  joint <- stats::integrate(f, -9, 9, kind = "joint",
                            subdivisions = 500L, stop.on.error = FALSE)$value
  same  <- stats::integrate(f, -9, 9, kind = "same",
                            subdivisions = 500L, stop.on.error = FALSE)$value
  c(same = 100 * same, direction = 100 * same / joint)
}


#' Student t(3) upper tail, scaled to unit variance
t3_tail <- function(q) stats::pt(q * sqrt(3), df = 3, lower.tail = FALSE)


#' Empirical upper-tail function for a vector of standardized errors
empirical_tail <- function(z) function(q) vapply(q, function(qq) mean(z > qq), numeric(1))


# --- Recalibration ----------------------------------------------------------

#' Leave-one-study-out isotonic recalibration
#'
#' The monotone counterpart to the linear recalibration used in the main
#' analysis, fitted with the same leave-one-study-out scheme so that no study
#' contributes to its own correction. Out-of-range predictions are held at the
#' training range (rule = 2).
#'
#' @param pred,est Predicted and observed effects.
#' @param study Study identifier, the unit left out.
#' @return Recalibrated predictions, in the order supplied.
iso_loo <- function(pred, est, study) {
  out <- rep(NA_real_, length(pred))
  for (s in unique(study)) {
    tr  <- study != s
    fit <- stats::isoreg(pred[tr], est[tr])
    f   <- stats::approxfun(sort(fit$x), fit$yf, method = "linear", rule = 2)
    out[!tr] <- f(pred[!tr])
  }
  out
}


#' Leave-one-study-out linear recalibration
#'
#' The authors' own remedy for the systematic overshooting they document,
#' applied here before any error is measured so that the shared error cannot be
#' attributed to that already-known bias. For each study the line is fitted on
#' the other studies only and applied to the held-out study.
#'
#' @param pred,est Predicted and observed effects.
#' @param study Study identifier, the unit left out.
#' @return Recalibrated predictions, in the order supplied.
linear_loo <- function(pred, est, study) {
  out <- rep(NA_real_, length(pred))
  d   <- data.frame(prediction = pred, estimate = est)
  for (s in unique(study)) {
    tr  <- study != s
    fit <- stats::lm(estimate ~ prediction, data = d[tr, , drop = FALSE])
    out[!tr] <- stats::predict(fit, newdata = d[!tr, , drop = FALSE])
  }
  out
}


# --- Summaries --------------------------------------------------------------

#' Joint-miss, same-side and direction rates with their exact nulls
#'
#' @param a,b Standardized errors of the two predictors.
#' @param thr Threshold on |z| defining a miss.
#' @return One-row data frame of observed rates and their corrected nulls.
summarise_shared <- function(a, b, thr = 1.96) {
  j <- abs(a) > thr & abs(b) > thr
  s <- j & sign(a) == sign(b)
  nulls <- null_exact(mean(abs(a) > thr), mean(abs(b) > thr), t = thr)
  data.frame(
    miss_g         = mean(abs(a) > thr),
    miss_h         = mean(abs(b) > thr),
    joint          = mean(j),
    null_joint     = nulls[["joint"]],
    same_side      = mean(s),
    null_same      = nulls[["same"]],
    direction      = sum(s) / sum(j),
    null_direction = nulls[["direction"]]
  )
}


#' Cluster-robust standard error for a slope, clustered by group
#'
#' Comparisons within a study share participants and often a reference
#' condition, so the conventional standard error understates uncertainty.
#'
#' @param fit A fitted lm object.
#' @param cluster Cluster identifier, one per row of the model frame.
#' @param term Index of the coefficient (2 = slope in a simple regression).
#' @return The cluster-robust standard error of that coefficient.
cluster_se <- function(fit, cluster, term = 2) {
  X    <- stats::model.matrix(fit)
  u    <- stats::residuals(fit)
  XtXi <- solve(crossprod(X))
  meat <- Reduce(`+`, lapply(split(seq_len(nrow(X)), cluster), function(i) {
    v <- crossprod(X[i, , drop = FALSE], u[i])
    tcrossprod(v)
  }))
  sqrt(diag(XtXi %*% meat %*% XtXi))[term]
}


# --- Reporting --------------------------------------------------------------

#' Accumulate a value quoted in the manuscript
#'
#' Every quantity that appears in the commentary or its Supplementary
#' Information is registered here, keyed to where it appears, and written to
#' output/tables/reported_values.csv at the end of the run. This lets a reader
#' check any number in the paper against the code that produced it without
#' executing the pipeline themselves.
#'
#' @param location Where the value appears, e.g. "Main text, para 5".
#' @param label What the value is.
#' @param value The value.
#' @param digits Rounding applied before writing out.
.reported <- new.env(parent = emptyenv())
.reported$rows <- list()

report_value <- function(location, label, value, digits = 3) {
  .reported$rows[[length(.reported$rows) + 1L]] <- data.frame(
    location = location,
    label    = label,
    value    = round(as.numeric(value), digits),
    stringsAsFactors = FALSE
  )
  invisible(value)
}

#' Write everything registered by report_value() to disk
write_reported_values <- function(path) {
  if (!length(.reported$rows)) {
    warning("No values registered; nothing written to ", path)
    return(invisible(NULL))
  }
  out <- do.call(rbind, .reported$rows)
  utils::write.csv(out, path, row.names = FALSE)
  message("Wrote ", nrow(out), " reported values to ", path)
  invisible(out)
}
