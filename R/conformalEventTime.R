#' Conformal lower prediction bound for the residual event time
#'
#' @description
#' For each test patient event-free at the landmark \eqn{s}, a lower bound
#' \eqn{\hat L} such that \eqn{P(T - s \ge \hat L \mid T > s) \ge 1 - \alpha}:
#' "with probability \eqn{1-\alpha} this patient stays event-free for at
#' least \eqn{\hat L} more years". The bound is a calibrated quantile of the
#' BJM predictive distribution of \eqn{T}, \eqn{\hat F(u) =}
#' \code{predictRisk(horizon = u)}, evaluated on a grid of horizons up to
#' \code{max_horizon} (\eqn{c_0}); \eqn{\hat L \le c_0}.
#'
#' The score of a calibration patient is \eqn{V = \hat F(T - s)} if
#' \eqn{T - s < c_0} and \eqn{V = 1} otherwise (a patient event-free beyond
#' \eqn{c_0} is covered by any bound \eqn{\le c_0}), and
#' \eqn{\hat L = \min(\hat F^{-1}(\hat\tau), c_0)}, with \eqn{\hat\tau}
#' the lower conformal \eqn{\alpha}-quantile of the scores, so that
#' \eqn{T - s < \hat L} exactly when \eqn{V < \hat\tau} (Candes, Lei and
#' Ren 2023 use the same construction with a quantile-regression model).
#' Censored calibration patients are handled in one of two ways:
#' \describe{
#'   \item{\code{"conservative"}}{Every calibration patient is used, with
#'   the observed time \eqn{X = \min(T, C)} in place of \eqn{T}. Since
#'   \eqn{X \le T} and \eqn{\hat F} is non-decreasing, the scores are too small, so the bound is too low:
#'   coverage \eqn{\ge 1-\alpha} holds in finite samples, for any censoring
#'   mechanism (even informative) and any model. Also valid for
#'   interval-censored data with the last event-free visit as \eqn{X}. The
#'   price is a looser bound when many patients are censored before
#'   \eqn{c_0}.}
#'   \item{\code{"ipcw"} (default)}{Only patients whose \eqn{\min(T - s, c_0)}
#'   is observed are used, weighted by \eqn{1/\hat G(\min(T - s, c_0))}, the
#'   inverse Kaplan-Meier probability of being uncensored by then (weighted
#'   conformal, Tibshirani et al. 2019). Coverage is then approximate: it
#'   holds as the calibration set grows if censoring is independent of the
#'   event time and history (the same assumption as the Kaplan-Meier
#'   observed risks in \code{BJM::calibrationPlot()}).}
#' }
#'
#' @inheritParams conformalLongitudinal
#' @param max_horizon \eqn{c_0}: the bound is capped here. Choose it where
#'   follow-up is still adequate.
#' @param n_grid Number of horizons in \eqn{(0, c_0]} at which
#'   \code{predictRisk()} is evaluated.
#' @param censoring \code{"ipcw"} or \code{"conservative"}.
#' @param bandcount1,bandcount2 Passed to \code{predictRisk()}.
#' @return An object of class \code{"conformalEventTime"}: \code{bounds}
#'   (one row per test patient: \code{id}, \code{lower} and
#'   \code{native_lower} on the residual scale, \code{time}, \code{status},
#'   \code{covered}: \code{TRUE}/\code{FALSE}, or \code{NA} when censored
#'   before the bound), \code{tau}, \code{coverage_ipcw} (IPCW estimate of
#'   the test coverage when test outcomes are present) and the settings.
#' @export
conformalEventTime <- function(data_calib, data_test, long_fit_all,
                               survival_fit_all, prediction_time, max_horizon,
                               time_variable, survival_variable_all,
                               survival_trans_function, alpha = 0.1,
                               censoring = c("ipcw", "conservative"),
                               n_grid = 20, bandcount1 = 20, bandcount2 = 40) {
  censoring <- match.arg(censoring)
  cols <- surv_columns(survival_fit_all)
  s <- prediction_time
  c0 <- max_horizon
  grid <- seq(0, c0, length.out = n_grid + 1)
  id_variable <- all.vars(nlme::getGroupsFormula(long_fit_all$lfit[[1]]))[1]
  n_out <- length(long_fit_all$lfit)

  cdf_matrix <- function(data) {
    dl <- as_data_list(data, n_out)
    hist <- landmark_data(dl, s, time_variable, id_variable, cols$time)
    risk <- lapply(grid[-1], function(u) suppressWarnings(BJM::predictRisk(
      hist, long_fit_all, survival_fit_all, prediction_time = s, horizon = u,
      time_variable = time_variable, survival_variable_all = survival_variable_all,
      survival_trans_function = survival_trans_function,
      bandcount1 = bandcount1, bandcount2 = bandcount2))$risk_prob_1)
    Fm <- cbind(0, do.call(cbind, risk))
    # enforce monotonicity in u (grid error can make it wiggle slightly)
    Fm <- t(apply(Fm, 1, cummax))
    rownames(Fm) <- names(risk[[1]])
    out <- subject_outcome(dl[[1]], id_variable, cols)
    out <- out[match(rownames(Fm), out$id), ]
    list(F = Fm, time = out$time - s, status = out$status)
  }

  F_at <- function(Frow, u) stats::approx(grid, Frow, xout = pmin(u, c0), rule = 2)$y
  F_inv <- function(Frow, p) {
    if (!is.finite(p) || p <= 0) return(0)
    if (Frow[length(Frow)] < p) return(c0)
    i <- which(Frow >= p)[1]
    if (i == 1) return(0)
    grid[i - 1] + (p - Frow[i - 1]) / (Frow[i] - Frow[i - 1]) * (grid[i] - grid[i - 1])
  }

  cal <- cdf_matrix(data_calib)
  x <- cal$time; d <- cal$status
  if (censoring == "conservative") {
    V <- vapply(seq_along(x), function(i)
      if (x[i] < c0) F_at(cal$F[i, ], x[i]) else 1, numeric(1))
    w <- NULL
  } else {
    known <- d == 1 | x >= c0
    G <- censoring_survival(x, d)
    tt <- pmin(x, c0)
    V <- vapply(which(known), function(i)
      if (x[i] < c0) F_at(cal$F[i, ], x[i]) else 1, numeric(1))
    w <- 1 / G(tt[known])
  }
  tau <- conformal_lower_quantile(V, alpha, w)

  tst <- cdf_matrix(data_test)
  lower <- vapply(seq_len(nrow(tst$F)), function(i) F_inv(tst$F[i, ], tau), numeric(1))
  native <- vapply(seq_len(nrow(tst$F)), function(i) F_inv(tst$F[i, ], alpha), numeric(1))
  bounds <- data.frame(id = rownames(tst$F), lower = lower, native_lower = native,
                       time = tst$time, status = tst$status, stringsAsFactors = FALSE)
  cov <- coverage_lower_bound(bounds$lower, tst$time, tst$status)
  bounds$covered <- cov$covered
  nat <- coverage_lower_bound(bounds$native_lower, tst$time, tst$status)

  structure(list(bounds = bounds, tau = tau, scores = V, weights = w,
                 coverage_ipcw = cov$ipcw, native_coverage_ipcw = nat$ipcw,
                 alpha = alpha, censoring = censoring, prediction_time = s,
                 max_horizon = c0, n_calib = length(V)),
            class = "conformalEventTime")
}

# Whether T - s >= L, from right-censored residual times: known when the
# event is seen before L (FALSE) or follow-up reaches L (TRUE). IPCW
# estimate of the coverage uses P(known) = G(min(T, L)).
coverage_lower_bound <- function(L, time, status) {
  if (all(is.na(time))) return(list(covered = rep(NA, length(L)), ipcw = NA_real_))
  covered <- ifelse(time >= L, TRUE, ifelse(status == 1, FALSE, NA))
  G <- censoring_survival(time, status)
  tt <- pmin(time, L)
  w <- ifelse(is.na(covered), 0, 1 / G(tt))
  list(covered = covered,
       ipcw = sum(w * (covered %in% TRUE)) / sum(w))
}

#' @export
print.conformalEventTime <- function(x, ...) {
  b <- x$bounds
  cat(sprintf("Conformal lower bound on residual event time at s = %g (cap %g, alpha = %g, %s)\n",
              x$prediction_time, x$max_horizon, x$alpha, x$censoring))
  cat(sprintf("  calibration scores: %d, tau = %.4g (native alpha = %g)\n",
              x$n_calib, x$tau, x$alpha))
  cat(sprintf("  test patients: %d, median bound %.3g years (BJM native %.3g)\n",
              nrow(b), stats::median(b$lower), stats::median(b$native_lower)))
  if (is.finite(x$coverage_ipcw)) {
    cat(sprintf("  IPCW-estimated test coverage: %.3f (native %.3f)\n",
                x$coverage_ipcw, x$native_coverage_ipcw))
  }
  invisible(x)
}
