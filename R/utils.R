# Internal helpers shared by the conformal wrappers. Nothing here calls BJM
# internals: only the public predictLongitudinal()/predictRisk() results and
# the survivalSub() object's formula are used.

#' Conformal quantile of nonconformity scores
#'
#' The \eqn{\lceil (1-\alpha)(n+1) \rceil}-th smallest of \code{n} scores,
#' i.e. the \eqn{1-\alpha} quantile of the scores together with a point mass
#' at \code{+Inf} for the test point. With \code{weights}, the weighted
#' version of Tibshirani et al. (2019): the test point carries weight
#' \code{test_weight} (default: the mean calibration weight).
#'
#' @param scores Numeric vector of calibration scores (\code{NA} dropped).
#' @param alpha Miscoverage level.
#' @param weights Optional non-negative weights, same length as \code{scores}.
#' @param test_weight Weight given to the test point when \code{weights} is
#'   supplied.
#' @return A single number, possibly \code{Inf} when there are too few
#'   calibration points for level \code{alpha}.
conformalQuantile <- function(scores, alpha, weights = NULL, test_weight = NULL) {
  stopifnot(length(alpha) == 1, alpha > 0, alpha < 1)
  if (is.null(weights)) weights <- rep(1, length(scores))
  ok <- !is.na(scores) & !is.na(weights) & weights > 0
  scores <- scores[ok]; weights <- weights[ok]
  n <- length(scores)
  if (n == 0) return(Inf)
  if (is.null(test_weight)) test_weight <- mean(weights)
  o <- order(scores)
  scores <- scores[o]; weights <- weights[o]
  cum <- cumsum(weights) / (sum(weights) + test_weight)
  k <- which(cum >= 1 - alpha - 1e-12)
  if (length(k) == 0) Inf else scores[k[1]]
}

# Lower conformal quantile: the largest threshold tau with
# P(V_test < tau) <= alpha, from calibration scores V (test mass at -Inf).
conformal_lower_quantile <- function(scores, alpha, weights = NULL, test_weight = NULL) {
  -conformalQuantile(-scores, alpha, weights, test_weight)
}

#' Weighted isotonic (non-decreasing) regression by pool-adjacent-violators
#'
#' @param x Ordering variable.
#' @param y Responses.
#' @param w Positive weights.
#' @return Fitted values in the original order of \code{x}.
weightedIsotonic <- function(x, y, w = rep(1, length(y))) {
  o <- order(x, -y)  # ties: order does not matter for the pooled fit
  yy <- y[o]; ww <- w[o]
  n <- length(yy)
  val <- numeric(n); wt <- numeric(n); len <- integer(n); m <- 0L
  for (i in seq_len(n)) {
    m <- m + 1L; val[m] <- yy[i]; wt[m] <- ww[i]; len[m] <- 1L
    while (m > 1L && val[m - 1L] > val[m]) {
      tot <- wt[m - 1L] + wt[m]
      val[m - 1L] <- (wt[m - 1L] * val[m - 1L] + wt[m] * val[m]) / tot
      wt[m - 1L] <- tot; len[m - 1L] <- len[m - 1L] + len[m]
      m <- m - 1L
    }
  }
  fit <- rep(val[seq_len(m)], len[seq_len(m)])
  # pool exact ties in x so tied inputs get one fitted value
  xs <- x[o]
  fit <- stats::ave(fit * ww, xs, FUN = sum) / stats::ave(ww, xs, FUN = sum)
  out <- numeric(n); out[o] <- fit
  out
}

# Event time and status column names from a survivalSub() fit.
surv_columns <- function(survival_fit_all) {
  v <- all.vars(survival_fit_all$form_marginal_surv[[2]])
  if (length(v) != 2) {
    stop("Only right-censored Surv(time, status) survival sub-models are ",
         "supported for the event-time and risk wrappers.", call. = FALSE)
  }
  list(time = v[1], status = v[2])
}

as_data_list <- function(data, n_outcomes) {
  if (is.data.frame(data)) rep(list(data), n_outcomes) else data
}

# One row per subject: observed time and status (first row of each id).
subject_outcome <- function(data, id_variable, cols) {
  d <- data[!duplicated(data[[id_variable]]), , drop = FALSE]
  data.frame(id = as.character(d[[id_variable]]),
             time = d[[cols$time]],
             status = d[[cols$status]],
             stringsAsFactors = FALSE)
}

# Landmark subset: rows with time_variable <= s, for subjects with an
# observed (or unknown) event/censoring time > s.
landmark_data <- function(data_list, s, time_variable, id_variable, event_time) {
  lapply(data_list, function(d) {
    keep <- d[[time_variable]] <= s & (is.na(d[[event_time]]) | d[[event_time]] > s)
    d[keep, , drop = FALSE]
  })
}

subset_ids <- function(data_list, ids, id_variable) {
  lapply(data_list, function(d) d[as.character(d[[id_variable]]) %in% ids, , drop = FALSE])
}

# Kaplan-Meier of the censoring distribution on the residual time scale,
# returned as a function giving G(t-) = P(C >= t).
censoring_survival <- function(time, status) {
  km <- survival::survfit(survival::Surv(time, 1 - status) ~ 1)
  tt <- km$time; ss <- km$surv
  function(t) {
    idx <- findInterval(t, tt, left.open = TRUE)  # number of jump times < t
    out <- ifelse(idx == 0, 1, ss[pmax(idx, 1)])
    pmax(out, 1e-8)
  }
}

# ---- predictive-density helpers (continuous biomarker on a grid) ----------

normalize_density <- function(y, f) {
  f <- pmax(f, 0)
  area <- sum(diff(y) * (utils::head(f, -1) + utils::tail(f, -1)) / 2)
  if (!is.finite(area) || area <= 0) return(rep(NA_real_, length(f)))
  f / area
}

density_cdf <- function(y, f) {
  c(0, cumsum(diff(y) * (utils::head(f, -1) + utils::tail(f, -1)) / 2))
}

grid_quantile <- function(y, Fy, p) {
  vapply(p, function(pp) {
    if (pp <= 0) return(-Inf)
    if (pp >= 1) return(Inf)
    i <- which(Fy >= pp)[1]
    if (is.na(i)) return(max(y))
    if (i == 1) return(y[1])
    y[i - 1] + (pp - Fy[i - 1]) / (Fy[i] - Fy[i - 1]) * (y[i] - y[i - 1])
  }, numeric(1))
}

# HPD mass score: posterior mass of {f > f(y_obs)}.
hpd_mass <- function(y, f, y_obs) {
  f_obs <- if (y_obs < min(y) || y_obs > max(y)) 0 else
    stats::approx(y, f, xout = y_obs)$y
  g <- ifelse(f > f_obs, f, 0)
  sum(diff(y) * (utils::head(g, -1) + utils::tail(g, -1)) / 2)
}

# Smallest-density threshold c with mass{f >= c} >= q, and the region it
# defines (range, and whether it is a single interval on the grid).
hpd_region <- function(y, f, q) {
  if (!is.finite(q) || q >= 1) return(c(lower = -Inf, upper = Inf, connected = 1))
  lev <- sort(unique(f), decreasing = TRUE)
  mass <- vapply(lev, function(cc) {
    g <- ifelse(f >= cc, f, 0)
    sum(diff(y) * (utils::head(g, -1) + utils::tail(g, -1)) / 2)
  }, numeric(1))
  cc <- lev[which(mass >= q)[1]]
  if (is.na(cc)) cc <- min(lev)
  inside <- which(f >= cc)
  c(lower = y[min(inside)], upper = y[max(inside)],
    connected = as.numeric(all(diff(inside) == 1)))
}
