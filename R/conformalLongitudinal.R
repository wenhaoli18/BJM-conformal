#' Split-conformal prediction intervals for a future biomarker value
#'
#' @description
#' Calibrates the predictive density returned by
#' \code{BJM::predictLongitudinal()} on held-out patients, so that the
#' interval for a new patient covers that patient's biomarker value at
#' \code{prediction_time + horizon} with probability at least
#' \code{1 - alpha}, whatever the joint model gets wrong.
#'
#' \strong{Calibration patients.} A patient in \code{data_calib} is used if
#' they are event-free at \code{prediction_time} (\eqn{s}), have at least one
#' measurement of the biomarker at or before \eqn{s}, and have a measurement
#' at a time \eqn{t^*} within \code{window} of \eqn{s + h} (the one closest
#' to \eqn{s + h}). Their score is computed from
#' \code{predictLongitudinal()} with horizon \eqn{t^* - s}, using only their
#' history up to \eqn{s}: so it uses the same conditioning as the target
#' (history up to \eqn{s}, alive at the time predicted). Each patient
#' contributes one score, so scores are exchangeable across patients.
#'
#' \strong{Scores.}
#' \describe{
#'   \item{\code{"hpd"} (default)}{The predictive mass of the values more
#'   likely than the observed one, \eqn{\hat P(\hat f(Y) > \hat f(y))}
#'   (HPD-split, Izbicki, Shimizu and Stern 2022). The interval is an HPD
#'   region of the BJM density with calibrated mass, so it follows the
#'   density's skewness.}
#'   \item{\code{"cdf"}}{\eqn{|\hat F(y) - 1/2|} (distributional conformal
#'   prediction, Chernozhukov, Wuthrich and Zhu 2021). The interval is
#'   equal-tailed: \eqn{[\hat F^{-1}(1/2 - q), \hat F^{-1}(1/2 + q)]}.}
#'   \item{\code{"abs"}}{\eqn{|y - \hat m| / \hat\sigma}, with \eqn{\hat m}
#'   the BJM mode and \eqn{\hat\sigma} the predictive standard deviation.}
#' }
#' If the model is right, both \code{"hpd"} and \code{"cdf"} give back the
#' model's own interval (about the \code{native_*} columns); if it is wrong,
#' the calibrated quantile \code{qhat} moves to restore coverage.
#'
#' \strong{What is guaranteed.} With exchangeable calibration and test
#' patients, \eqn{P(Y \in \hat C) \ge 1 - \alpha}, marginally over
#' patients, for the fixed landmark and horizon, among patients alive and
#' measured at the time predicted. It is not conditional on a patient's
#' history. Informative measurement (being measured near \eqn{s+h}
#' depending on how sick one is, beyond what the history explains) breaks
#' exchangeability.
#'
#' @param bio_i Index of the biomarker in \code{long_fit_all} (one only).
#' @param data_calib,data_test Calibration and test data: full long-format
#'   histories (a data frame, or a list with one per biomarker, as in
#'   \code{BJM}). Patients in \code{data_calib} must not have been used to fit
#'   the sub-models. Rows after \code{prediction_time} in \code{data_test}
#'   are used only to report coverage, never for prediction.
#' @param long_fit_all,survival_fit_all,time_variable,survival_variable_all,survival_trans_function
#'   As for \code{BJM::predictLongitudinal()}.
#' @param prediction_time Landmark \eqn{s}.
#' @param horizon Horizon \eqn{h}.
#' @param alpha Miscoverage level.
#' @param score \code{"hpd"}, \code{"cdf"} or \code{"abs"}; several may be
#'   given, in which case the BJM predictions are computed once and a named
#'   list of results, one per score, is returned.
#' @param window Half-width of the window around \eqn{s+h} in which a
#'   measurement counts as the outcome.
#' @param bandcount2,bandcount3 Passed to \code{predictLongitudinal()}; fixed
#'   by default, so each patient's call is fast and deterministic.
#' @return An object of class \code{"conformalLongitudinal"}: a list with
#'   \code{intervals} (one row per test patient: \code{id}, \code{horizon}
#'   used, \code{prediction} (BJM mode), \code{lower}, \code{upper},
#'   \code{native_lower}, \code{native_upper} (the model's own
#'   equal-tailed interval), \code{y_obs} and \code{covered},
#'   \code{native_covered} when the test data hold a measurement in the
#'   window, and \code{connected} for the \code{"hpd"} region), \code{qhat},
#'   the calibration \code{scores}, and the settings.
#' @export
conformalLongitudinal <- function(bio_i, data_calib, data_test, long_fit_all,
                                  survival_fit_all, prediction_time, horizon,
                                  time_variable, survival_variable_all,
                                  survival_trans_function, alpha = 0.1,
                                  score = c("hpd", "cdf", "abs"),
                                  window = 0.5, bandcount2 = 40,
                                  bandcount3 = 200) {
  score_all <- match.arg(score, several.ok = TRUE)
  if (missing(score)) score_all <- "hpd"
  stopifnot(length(bio_i) == 1, horizon > 0, window >= 0)
  n_out <- length(long_fit_all$lfit)
  lfit <- long_fit_all$lfit[[bio_i]]
  bio_name <- as.character(stats::formula(lfit)[[2]])
  id_variable <- all.vars(nlme::getGroupsFormula(lfit))[1]
  event_time <- surv_columns_any(survival_fit_all)
  s <- prediction_time

  predict_one <- function(hist_list, hz) {
    p <- suppressWarnings(BJM::predictLongitudinal(
      bio_i, hist_list, long_fit_all, survival_fit_all,
      prediction_time = s, horizon = hz, time_variable = time_variable,
      survival_variable_all = survival_variable_all,
      survival_trans_function = survival_trans_function,
      bandcount2 = bandcount2, bandcount3 = bandcount3))
    y <- p$Y_all
    f <- normalize_density(y, p$Y_density[, 1])
    Fy <- density_cdf(y, f)
    m <- sum(diff(y) * (utils::head(y * f, -1) + utils::tail(y * f, -1)) / 2)
    v <- sum(diff(y) * (utils::head((y - m)^2 * f, -1) + utils::tail((y - m)^2 * f, -1)) / 2)
    list(y = y, f = f, Fy = Fy, mode = unname(p$Y_predict[1]), sd = sqrt(v))
  }

  score_of <- function(pr, y_obs, sc) {
    switch(sc,
      hpd = hpd_mass(pr$y, pr$f, y_obs),
      cdf = abs(stats::approx(pr$y, pr$Fy, xout = y_obs, rule = 2)$y - 0.5),
      abs = abs(y_obs - pr$mode) / pr$sd)
  }

  interval_of <- function(pr, q, sc) {
    switch(sc,
      hpd = hpd_region(pr$y, pr$f, q),
      cdf = c(grid_quantile(pr$y, pr$Fy, c(0.5 - q, 0.5 + q)), connected = 1),
      abs = c(pr$mode - q * pr$sd, pr$mode + q * pr$sd, connected = 1))
  }

  # For each patient: history up to s (all biomarkers) and the outcome
  # measurement of bio_i closest to s + h within the window.
  prepare <- function(data) {
    dl <- as_data_list(data, n_out)
    db <- dl[[bio_i]]
    first <- db[!duplicated(db[[id_variable]]), , drop = FALSE]
    at_risk <- as.character(first[[id_variable]][is.na(first[[event_time]]) |
                                                  first[[event_time]] > s])
    hist <- landmark_data(dl, s, time_variable, id_variable, event_time)
    has_hist <- unique(as.character(hist[[bio_i]][[id_variable]][!is.na(hist[[bio_i]][[bio_name]])]))
    ids <- intersect(at_risk, has_hist)
    fut <- db[db[[time_variable]] > s & !is.na(db[[bio_name]]) &
                abs(db[[time_variable]] - (s + horizon)) <= window, , drop = FALSE]
    fut <- fut[order(abs(fut[[time_variable]] - (s + horizon))), , drop = FALSE]
    fut <- fut[!duplicated(fut[[id_variable]]), , drop = FALSE]
    target <- stats::setNames(fut[[bio_name]], as.character(fut[[id_variable]]))
    t_star <- stats::setNames(fut[[time_variable]], as.character(fut[[id_variable]]))
    list(hist = hist, ids = ids, target = target, t_star = t_star)
  }

  cal <- prepare(data_calib)
  cal_ids <- intersect(cal$ids, names(cal$target))
  if (length(cal_ids) < ceiling(1 / alpha) - 1) {
    warning(sprintf("Only %d calibration patients: the interval will be infinite at alpha = %g.",
                    length(cal_ids), alpha), call. = FALSE)
  }
  # BJM predictions are computed once and shared by every requested score
  cal_pred <- lapply(cal_ids, function(i)
    predict_one(subset_ids(cal$hist, i, id_variable), cal$t_star[[i]] - s))
  tst <- prepare(data_test)
  tst_pred <- lapply(tst$ids, function(i) {
    hz <- if (i %in% names(tst$target)) tst$t_star[[i]] - s else horizon
    pr <- predict_one(subset_ids(tst$hist, i, id_variable), hz)
    pr$horizon <- hz
    pr
  })

  one_score <- function(sc) {
    scores <- stats::setNames(vapply(seq_along(cal_ids), function(k)
      score_of(cal_pred[[k]], cal$target[[cal_ids[k]]], sc), numeric(1)), cal_ids)
    qhat <- conformalQuantile(scores, alpha)
    rows <- lapply(seq_along(tst$ids), function(k) {
      i <- tst$ids[k]; pr <- tst_pred[[k]]
      has_y <- i %in% names(tst$target)
      ci <- interval_of(pr, qhat, sc)
      nat <- grid_quantile(pr$y, pr$Fy, c(alpha / 2, 1 - alpha / 2))
      y_obs <- if (has_y) tst$target[[i]] else NA_real_
      covered <- if (!has_y) NA else if (sc == "hpd" && is.finite(qhat))
        score_of(pr, y_obs, sc) <= qhat else y_obs >= ci[[1]] && y_obs <= ci[[2]]
      data.frame(id = i, horizon = pr$horizon, prediction = pr$mode,
                 lower = ci[[1]], upper = ci[[2]],
                 native_lower = nat[1], native_upper = nat[2],
                 y_obs = y_obs, covered = covered,
                 native_covered = if (has_y) y_obs >= nat[1] && y_obs <= nat[2] else NA,
                 connected = as.logical(ci[[3]]), stringsAsFactors = FALSE)
    })
    intervals <- do.call(rbind, rows)
    rownames(intervals) <- NULL
    structure(list(intervals = intervals, qhat = qhat, scores = scores,
                   n_calib = length(scores), alpha = alpha, score = sc,
                   biomarker = bio_name, prediction_time = s, horizon = horizon,
                   window = window),
              class = "conformalLongitudinal")
  }
  if (length(score_all) == 1) return(one_score(score_all))
  stats::setNames(lapply(score_all, one_score), score_all)
}

# Event time column, also for interval-censored Surv(L, R, type = ...) fits.
surv_columns_any <- function(survival_fit_all) {
  all.vars(survival_fit_all$form_marginal_surv[[2]])[1]
}

#' @export
print.conformalLongitudinal <- function(x, ...) {
  iv <- x$intervals
  cat(sprintf("Split-conformal interval for %s at s = %g + %g (score '%s', alpha = %g)\n",
              x$biomarker, x$prediction_time, x$horizon, x$score, x$alpha))
  cat(sprintf("  calibration patients: %d, qhat = %.4g\n", x$n_calib, x$qhat))
  cat(sprintf("  test patients: %d, median width %.3g (BJM native %.3g)\n",
              nrow(iv), stats::median(iv$upper - iv$lower),
              stats::median(iv$native_upper - iv$native_lower)))
  if (any(!is.na(iv$covered))) {
    cat(sprintf("  empirical coverage on %d measured test patients: %.3f (native %.3f)\n",
                sum(!is.na(iv$covered)), mean(iv$covered, na.rm = TRUE),
                mean(iv$native_covered, na.rm = TRUE)))
  }
  invisible(x)
}
