# Data generator for the coverage simulations.
#
# Data are generated FORWARD from a shared-random-effects joint model
# (biomarker trajectory drives the hazard), so the BJM backward model
# Y | T is never the true model. Two scenarios:
#
#   "mild"   : linear trajectory, Gaussian errors, hazard on the current value.
#   "severe" : quadratic trajectory, heavy-tailed (t_3) heteroscedastic errors
#              whose SD grows with the current value, hazard on the current
#              value and the slope, plus a covariate x3 that affects both the
#              biomarker and the hazard but is omitted from both sub-models.
#
# Columns mirror pbc3: id, year (visit time), years (observed event time),
# status, y (biomarker), x1, x2, Tyears1..3 (|years - k|, unused by the
# fitted formulas but required by BJM's survival_trans interface).

simulate_jm <- function(n, scenario = c("mild", "severe"), id_offset = 0,
                        cens = TRUE, dt = 0.02, t_max = 20) {
  scenario <- match.arg(scenario)
  sev <- scenario == "severe"
  x1 <- stats::rnorm(n); x2 <- stats::rbinom(n, 1, 0.5)
  x3 <- if (sev) stats::rbinom(n, 1, 0.3) else rep(0, n)
  b0 <- stats::rnorm(n, 0, 0.8); b1 <- stats::rnorm(n, 0, 0.3)

  mean_fun <- function(i, t) {
    1 + b0[i] + (0.3 + b1[i]) * t + 0.3 * x1[i] - 0.3 * x2[i] +
      if (sev) 0.08 * t^2 + 0.8 * x3[i] else 0
  }
  slope_fun <- function(i, t) (0.3 + b1[i]) + if (sev) 0.16 * t else 0
  err <- function(m) {
    if (!sev) return(stats::rnorm(length(m), 0, 0.5))
    0.35 * exp(0.25 * pmin(m, 4)) * stats::rt(length(m), df = 3) / sqrt(3)
  }

  tg <- seq(dt, t_max, by = dt)
  T_true <- vapply(seq_len(n), function(i) {
    m <- mean_fun(i, tg)
    lp <- (if (sev) -4.6 else -3.6) + 0.7 * m + 0.3 * x2[i] +
      if (sev) 1.2 * slope_fun(i, tg) + 0.8 * x3[i] else 0
    H <- cumsum(exp(lp) * dt)
    k <- which(H >= stats::rexp(1))[1]
    if (is.na(k)) t_max else tg[k]
  }, numeric(1))

  C <- if (cens) pmin(stats::runif(n, 3, 10), stats::rexp(n, 0.05)) else rep(Inf, n)
  obs <- pmin(T_true, C)
  status <- as.numeric(T_true <= C)

  rows <- lapply(seq_len(n), function(i) {
    v <- seq(0, t_max, by = 0.5)
    v <- v + c(0, stats::runif(length(v) - 1, -0.1, 0.1))
    v <- v[v < obs[i]]
    m <- mean_fun(i, v)
    data.frame(id = id_offset + i, year = v, y = m + err(m),
               years = obs[i], status = status[i], x1 = x1[i], x2 = x2[i])
  })
  d <- do.call(rbind, rows)
  d <- add_trans(d)
  attr(d, "truth") <- list(T = T_true, mean_fun = mean_fun, err = err,
                           id = id_offset + seq_len(n))
  d
}

add_trans <- function(d) {
  for (k in 1:3) d[[paste0("Tyears", k)]] <- abs(d$years - c(1, 3, 5)[k])
  d
}

# Test set: true (uncensored) event times, history up to s, and one extra
# row holding the biomarker drawn exactly at s + h for patients alive then.
make_test <- function(n, scenario, s, h, id_offset) {
  d <- simulate_jm(n, scenario, id_offset = id_offset, cens = FALSE)
  tr <- attr(d, "truth")
  d <- d[d$year <= s, ]
  alive <- which(tr$T > s + h)
  extra <- d[!duplicated(d$id) & d$id %in% tr$id[alive], ]
  i_ext <- match(extra$id, tr$id)
  m <- vapply(i_ext, function(i) tr$mean_fun(i, s + h), numeric(1))
  extra$year <- s + h
  extra$y <- m + tr$err(m)
  d <- rbind(d, extra)
  d[order(d$id, d$year), ]
}
