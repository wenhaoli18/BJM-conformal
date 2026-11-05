# pbc3: repeated random splits (fit 50% / calibrate 25% / test 25%).
#
#   Rscript analysis/pbc3_splits.R <n_splits> [n_cores]
#
# Model as in the BJM vignette and BJM-attribution: composite endpoint,
# Surv(years, status3) ~ age + sex; four biomarkers, each
# ~ year + age + sex + years + years * year, random ~ year | id, fitted on
# patients with an observed event. Landmark s = 3, horizon h = 2.
# Test outcomes are censored, so coverage is estimated: biomarker coverage
# on test patients measured within 0.5 years of s + h, event-time coverage
# by IPCW. Saves results/pbc3_splits.rds.

suppressMessages({library(BJM); library(BJMconformal); library(parallel)})
args <- commandArgs(trailingOnly = TRUE)
n_splits <- if (length(args) >= 1) as.integer(args[1]) else 2
n_cores <- if (length(args) >= 2) as.integer(args[2]) else 1

data(pbc3)
s <- 3; h <- 2; alpha <- 0.1
tr <- survivalTrans(c(1, 3, 5, 7))
fixed <- list(l1 = serBilir ~ year + age + sex + years + years * year,
              l2 = albumin ~ year + age + sex + years + years * year,
              l3 = prothrombin ~ year + age + sex + years + years * year,
              l4 = platelets ~ year + age + sex + years + years * year)
random <- rep(list(~ year | id), 4); names(random) <- names(fixed)
ids <- unique(as.character(pbc3$id))

one_split <- function(k) {
  set.seed(500 + k)
  g <- sample(rep(c(1, 1, 2, 3), length.out = length(ids))); names(g) <- ids
  gi <- g[as.character(pbc3$id)]
  fit <- pbc3[gi == 1, ]; cal <- pbc3[gi == 2, ]; tst <- pbc3[gi == 3, ]
  sf <- survivalSub(fit[!duplicated(fit$id), ], Surv(years, status3) ~ age + sex, NULL)
  lf <- longitudinalSub(fit[fit$status3 == 1, ], fixed, random)
  a <- list(long_fit_all = lf, survival_fit_all = sf, prediction_time = s,
            time_variable = "year", survival_variable_all = tr$survival_variable_all,
            survival_trans_function = tr$survival_trans_function)

  long <- do.call(rbind, lapply(1:2, function(b) {
    r <- do.call(conformalLongitudinal, c(list(bio_i = b, data_calib = cal, data_test = tst,
            horizon = h, alpha = alpha, score = c("hpd", "cdf"), window = 0.5), a))
    iv <- r$hpd$intervals; iv <- iv[!is.na(iv$covered), ]
    rbind(data.frame(biomarker = r$hpd$biomarker, method = "native",
                     coverage = mean(iv$native_covered), n_test = nrow(iv),
                     width = median(iv$native_upper - iv$native_lower), n_calib = NA),
          do.call(rbind, lapply(names(r), function(sc) {
            iv <- r[[sc]]$intervals; iv <- iv[!is.na(iv$covered), ]
            data.frame(biomarker = r[[sc]]$biomarker, method = sc,
                       coverage = mean(iv$covered), n_test = nrow(iv),
                       width = median(iv$upper - iv$lower), n_calib = r[[sc]]$n_calib)
          })))
  }))

  et <- lapply(c("ipcw", "conservative"), function(cm)
    do.call(conformalEventTime, c(list(data_calib = cal, data_test = tst, max_horizon = 4,
            alpha = alpha, censoring = cm), a)))
  event <- data.frame(method = c("native", "ipcw", "conservative"),
                      coverage_ipcw = c(et[[1]]$native_coverage_ipcw, et[[1]]$coverage_ipcw,
                                        et[[2]]$coverage_ipcw),
                      median_bound = c(median(et[[1]]$bounds$native_lower),
                                       median(et[[1]]$bounds$lower),
                                       median(et[[2]]$bounds$lower)))

  va <- do.call(vennAbersRisk, c(list(data_calib = cal, data_test = tst, horizon = h), a))
  rk <- va$risks
  km <- survival::survfit(survival::Surv(time, status) ~ 1, data = rk)
  obs <- 1 - summary(km, times = h, extend = TRUE)$surv
  risk <- data.frame(method = c("bjm", "venn_abers"),
                     mean_pred = c(mean(rk$p_bjm), mean(rk$p_merged)),
                     observed_km = obs, width = c(NA, median(rk$p1 - rk$p0)))
  list(split = k, long = long, event = event, risk = risk)
}

t0 <- Sys.time()
res <- mclapply(seq_len(n_splits), function(k) tryCatch(one_split(k), error = function(e) e),
                mc.cores = n_cores)
err <- vapply(res, inherits, logical(1), "error")
if (any(err)) message(sum(err), " split(s) failed: ", conditionMessage(res[[which(err)[1]]]))
dir.create("results", showWarnings = FALSE)
saveRDS(list(s = s, h = h, alpha = alpha, results = res[!err]), "results/pbc3_splits.rds")
message(sprintf("pbc3: %d splits in %.1f min", sum(!err),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
