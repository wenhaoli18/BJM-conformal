#' Venn-Abers interval-valued risks for BJM dynamic predictions
#'
#' @description
#' Recalibrates \code{BJM::predictRisk()} on held-out patients with
#' inductive Venn-Abers prediction (Vovk and Petej 2014). For each test
#' patient with BJM risk \eqn{p}, an isotonic regression of the observed
#' window status \eqn{Z = 1\{s < T \le s + h\}} on the calibration risks is
#' fitted twice, once with the test patient added as \eqn{(p, 0)} and once as
#' \eqn{(p, 1)}; its fitted values at \eqn{p} are \eqn{p_0 \le p_1}.
#'
#' \strong{What is guaranteed.} No method can give a distribution-free
#' interval that covers a patient's true risk (it is never observed). What
#' Venn-Abers does guarantee, without a censoring and for exchangeable
#' patients, is calibration: one of \eqn{p_0, p_1} is perfectly calibrated,
#' \eqn{E[Z \mid p_k] = p_k}. The width \eqn{p_1 - p_0} is the uncertainty
#' in the calibration itself, and shrinks like \eqn{1/n_{calib}}.
#' \code{p_merged} \eqn{= p_1 / (1 - p_0 + p_1)} is a single recalibrated risk.
#' Ranking (AUC) is unchanged within ties, as the map is monotone.
#'
#' \strong{Censoring.} \eqn{Z} is known for patients with the event in the
#' window or follow-up beyond \eqn{s + h}. With \code{"ipcw"} (default) these
#' are weighted by \eqn{1/\hat G(\min(T - s, h))} (Kaplan-Meier of the
#' censoring), so the guarantee becomes approximate (independent
#' censoring); \code{"complete"} drops censored patients unweighted and is
#' biased towards higher risk when censoring is heavy.
#'
#' @inheritParams conformalLongitudinal
#' @param censoring \code{"ipcw"} or \code{"complete"}.
#' @param bandcount1,bandcount2 Passed to \code{predictRisk()}.
#' @return An object of class \code{"vennAbersRisk"}: \code{risks} (one row
#'   per test patient: \code{id}, \code{p_bjm}, \code{p0}, \code{p1},
#'   \code{p_merged}, and when test outcomes are present \code{time},
#'   \code{status}, \code{z}) and the calibration data used.
#' @export
vennAbersRisk <- function(data_calib, data_test, long_fit_all, survival_fit_all,
                          prediction_time, horizon, time_variable,
                          survival_variable_all, survival_trans_function,
                          censoring = c("ipcw", "complete"),
                          bandcount1 = 20, bandcount2 = 40) {
  censoring <- match.arg(censoring)
  cols <- surv_columns(survival_fit_all)
  s <- prediction_time
  id_variable <- all.vars(nlme::getGroupsFormula(long_fit_all$lfit[[1]]))[1]
  n_out <- length(long_fit_all$lfit)

  risk_and_outcome <- function(data) {
    dl <- as_data_list(data, n_out)
    hist <- landmark_data(dl, s, time_variable, id_variable, cols$time)
    p <- suppressWarnings(BJM::predictRisk(
      hist, long_fit_all, survival_fit_all, prediction_time = s, horizon = horizon,
      time_variable = time_variable, survival_variable_all = survival_variable_all,
      survival_trans_function = survival_trans_function,
      bandcount1 = bandcount1, bandcount2 = bandcount2))$risk_prob_1
    out <- subject_outcome(dl[[1]], id_variable, cols)
    out <- out[match(names(p), out$id), ]
    x <- out$time - s
    z <- ifelse(out$status == 1 & x <= horizon, 1, ifelse(x > horizon, 0, NA))
    data.frame(id = names(p), p_bjm = unname(p), time = x, status = out$status,
               z = z, stringsAsFactors = FALSE)
  }

  cal <- risk_and_outcome(data_calib)
  known <- !is.na(cal$z)
  if (censoring == "ipcw") {
    G <- censoring_survival(cal$time, cal$status)
    w <- 1 / G(pmin(cal$time[known], horizon))
  } else {
    w <- rep(1, sum(known))
  }
  xc <- cal$p_bjm[known]; zc <- cal$z[known]; wt <- mean(w)

  tst <- risk_and_outcome(data_test)
  va <- t(vapply(tst$p_bjm, function(p) {
    f0 <- weightedIsotonic(c(xc, p), c(zc, 0), c(w, wt))
    f1 <- weightedIsotonic(c(xc, p), c(zc, 1), c(w, wt))
    c(f0[length(f0)], f1[length(f1)])
  }, numeric(2)))
  tst$p0 <- va[, 1]; tst$p1 <- va[, 2]
  tst$p_merged <- tst$p1 / (1 - tst$p0 + tst$p1)
  risks <- tst[, c("id", "p_bjm", "p0", "p1", "p_merged", "time", "status", "z")]

  structure(list(risks = risks, calib = cal[known, ], weights = w,
                 censoring = censoring, prediction_time = s, horizon = horizon),
            class = "vennAbersRisk")
}

#' @export
print.vennAbersRisk <- function(x, ...) {
  r <- x$risks
  cat(sprintf("Venn-Abers risks for (%g, %g] (%s, %d calibration patients with known status)\n",
              x$prediction_time, x$prediction_time + x$horizon, x$censoring, nrow(x$calib)))
  cat(sprintf("  test patients: %d, median width p1 - p0 = %.3f\n",
              nrow(r), stats::median(r$p1 - r$p0)))
  cat(sprintf("  mean risk: BJM %.3f, merged %.3f\n", mean(r$p_bjm), mean(r$p_merged)))
  invisible(x)
}
