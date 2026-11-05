# Coverage simulation: BJM native intervals vs. conformal wrappers.
#
#   Rscript sim/run_sim.R <scenario> <n_reps> [n_cores]
#
# Each replicate: fit BJM on n_fit patients, calibrate on n_cal, evaluate on
# n_test patients whose outcomes are fully known (no censoring in the test
# set, biomarker drawn exactly at s + h). Saves results/sim_<scenario>.rds.

suppressMessages({library(BJM); library(BJMconformal); library(parallel)})
source("sim/simulate.R")

args <- commandArgs(trailingOnly = TRUE)
scenario <- if (length(args) >= 1) args[1] else "mild"
n_reps <- if (length(args) >= 2) as.integer(args[2]) else 4
n_cores <- if (length(args) >= 3) as.integer(args[3]) else 4

s <- 2; h <- 2; alpha <- 0.1
n_fit <- 400; n_cal <- 300; n_test <- 300
tr <- survivalTrans(c(1, 3, 5))

one_rep <- function(r) {
  set.seed(1000 + r)
  fit <- simulate_jm(n_fit, scenario, id_offset = 0)
  cal <- simulate_jm(n_cal, scenario, id_offset = 1e5)
  tst <- make_test(n_test, scenario, s, h, id_offset = 2e5)

  sf <- survivalSub(fit[!duplicated(fit$id), ], Surv(years, status) ~ x1 + x2, NULL)
  lf <- longitudinalSub(fit[fit$status == 1, ],
                        list(l1 = y ~ year + x1 + x2 + years + years * year),
                        list(l1 = ~ year | id))
  args_bjm <- list(long_fit_all = lf, survival_fit_all = sf, prediction_time = s,
                   time_variable = "year", survival_variable_all = tr$survival_variable_all,
                   survival_trans_function = tr$survival_trans_function)

  # --- biomarker intervals
  bl <- do.call(conformalLongitudinal, c(list(bio_i = 1, data_calib = cal, data_test = tst,
          horizon = h, alpha = alpha, score = c("hpd", "cdf", "abs"), window = 0.25), args_bjm))
  long <- do.call(rbind, lapply(names(bl), function(sc) {
    iv <- bl[[sc]]$intervals; iv <- iv[!is.na(iv$covered), ]
    data.frame(method = sc, coverage = mean(iv$covered),
               width = median(iv$upper - iv$lower), n_calib = bl[[sc]]$n_calib)
  }))
  iv <- bl$hpd$intervals; iv <- iv[!is.na(iv$covered), ]
  long <- rbind(data.frame(method = "native", coverage = mean(iv$native_covered),
                           width = median(iv$native_upper - iv$native_lower),
                           n_calib = NA), long)
  # coverage by tertile of predicted value (a check of conditional coverage)
  tert <- cut(iv$prediction, quantile(iv$prediction, 0:3 / 3), include.lowest = TRUE, labels = 1:3)
  long_cond <- rbind(
    data.frame(method = "native", tertile = 1:3, coverage = tapply(iv$native_covered, tert, mean)),
    data.frame(method = "hpd", tertile = 1:3, coverage = tapply(iv$covered, tert, mean)))

  # --- event-time lower bounds (test outcomes uncensored: exact coverage)
  et <- lapply(c("ipcw", "conservative"), function(cm)
    do.call(conformalEventTime, c(list(data_calib = cal, data_test = tst, max_horizon = 4,
            alpha = alpha, censoring = cm), args_bjm)))
  b <- et[[1]]$bounds
  event <- data.frame(
    method = c("native", "ipcw", "conservative"),
    coverage = c(mean(b$time >= b$native_lower), mean(b$time >= b$lower),
                 mean(et[[2]]$bounds$time >= et[[2]]$bounds$lower)),
    median_bound = c(median(b$native_lower), median(b$lower), median(et[[2]]$bounds$lower)))

  # --- Venn-Abers risks
  va <- do.call(vennAbersRisk, c(list(data_calib = cal, data_test = tst, horizon = h), args_bjm))
  rk <- va$risks
  ece <- function(p, z) {
    g <- cut(p, unique(quantile(p, 0:10 / 10)), include.lowest = TRUE)
    sum(tapply(p, g, length) * abs(tapply(p, g, mean) - tapply(z, g, mean))) / length(p)
  }
  risk <- data.frame(
    method = c("bjm", "venn_abers"),
    citl = c(mean(rk$p_bjm) - mean(rk$z), mean(rk$p_merged) - mean(rk$z)),
    ece = c(ece(rk$p_bjm, rk$z), ece(rk$p_merged, rk$z)),
    brier = c(mean((rk$p_bjm - rk$z)^2), mean((rk$p_merged - rk$z)^2)),
    # does [p0, p1] contain the observed event rate of patients with similar p?
    width = c(NA, median(rk$p1 - rk$p0)))

  list(rep = r, long = long, long_cond = long_cond, event = event, risk = risk)
}

t0 <- Sys.time()
res <- mclapply(seq_len(n_reps), function(r) tryCatch(one_rep(r), error = function(e) e),
                mc.cores = n_cores)
err <- vapply(res, inherits, logical(1), "error")
if (any(err)) message(sum(err), " replicate(s) failed: ", conditionMessage(res[[which(err)[1]]]))
res <- res[!err]
dir.create("results", showWarnings = FALSE)
saveRDS(list(scenario = scenario, s = s, h = h, alpha = alpha, results = res),
        sprintf("results/sim_%s.rds", scenario))
message(sprintf("%s: %d replicates in %.1f min", scenario, length(res),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
