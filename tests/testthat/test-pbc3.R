test_that("the wrappers run end to end on pbc3", {
  skip_on_cran()
  suppressMessages(library(BJM))
  data(pbc3, package = "BJM")
  ids <- unique(as.character(pbc3$id))
  set.seed(3); g <- sample(rep(1:3, length.out = length(ids))); names(g) <- ids
  gi <- g[as.character(pbc3$id)]
  fit <- pbc3[gi == 1, ]; cal <- pbc3[gi == 2, ]; tst <- pbc3[gi == 3, ]
  sf <- survivalSub(fit[!duplicated(fit$id), ], Surv(years, status3) ~ age + sex, NULL)
  lf <- longitudinalSub(fit[fit$status3 == 1, ],
                        list(l1 = serBilir ~ year + age + sex + years),
                        list(l1 = ~ year | id))
  tr <- survivalTrans(c(1, 3, 5, 7))
  a <- list(long_fit_all = lf, survival_fit_all = sf, prediction_time = 3,
            time_variable = "year", survival_variable_all = tr$survival_variable_all,
            survival_trans_function = tr$survival_trans_function)
  r <- do.call(conformalLongitudinal, c(list(bio_i = 1, data_calib = cal, data_test = tst,
                horizon = 2, alpha = 0.2, score = c("hpd", "cdf")), a))
  expect_named(r, c("hpd", "cdf"))
  iv <- r$cdf$intervals
  expect_true(is.finite(r$cdf$qhat))
  expect_true(all(iv$lower <= iv$prediction + 1e-8 & iv$prediction <= iv$upper + 1e-8))
  e <- do.call(conformalEventTime, c(list(data_calib = cal, data_test = tst,
                max_horizon = 3, alpha = 0.2, n_grid = 6), a))
  expect_true(all(e$bounds$lower >= 0 & e$bounds$lower <= 3))
  v <- do.call(vennAbersRisk, c(list(data_calib = cal, data_test = tst, horizon = 2), a))
  expect_true(all(v$risks$p0 <= v$risks$p1 + 1e-12))
})
