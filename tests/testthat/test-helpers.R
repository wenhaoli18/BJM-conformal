test_that("conformalQuantile is the ceiling((n+1)(1-alpha))-th order statistic", {
  x <- c(5, 1, 4, 2, 3, 9, 8, 7, 6, 10)
  expect_equal(conformalQuantile(x, 0.1), 10)           # ceiling(11 * 0.9) = 10
  expect_equal(conformalQuantile(x[1:8], 0.1), Inf)     # needs n >= 9 at 0.1
  expect_equal(conformalQuantile(x, 0.2), 9)            # ceiling(11 * 0.8) = 9
  expect_equal(conformalQuantile(x, 0.2, weights = rep(3, 10)), 9)
  expect_equal(BJMconformal:::conformal_lower_quantile(x, 0.2), 2)
})

test_that("weightedIsotonic matches isoreg with unit weights and is monotone", {
  set.seed(1)
  x <- runif(50); y <- rbinom(50, 1, x)
  f <- weightedIsotonic(x, y)
  expect_equal(f[order(x)], stats::isoreg(x, y)$yf)
  fw <- weightedIsotonic(x, y, runif(50, 0.5, 2))
  expect_true(all(diff(fw[order(x)]) >= -1e-12))
})

test_that("HPD mass and region agree on a Gaussian density", {
  y <- seq(-6, 6, length.out = 2001)
  f <- BJMconformal:::normalize_density(y, dnorm(y))
  expect_equal(BJMconformal:::hpd_mass(y, f, qnorm(0.95)), 0.9, tolerance = 1e-3)
  r <- BJMconformal:::hpd_region(y, f, 0.9)
  expect_equal(unname(r[1:2]), c(-1, 1) * qnorm(0.95), tolerance = 1e-2)
  Fy <- BJMconformal:::density_cdf(y, f)
  expect_equal(BJMconformal:::grid_quantile(y, Fy, 0.975), qnorm(0.975), tolerance = 1e-3)
})

test_that("conformal calibration restores coverage of a wrong predictive density", {
  # predictive density N(0, 1), truth t_2: native 90% interval undercovers
  set.seed(2)
  y <- seq(-30, 30, length.out = 3001)
  f <- BJMconformal:::normalize_density(y, dnorm(y))
  cover <- replicate(300, {
    cal <- rt(99, 2)
    q <- conformalQuantile(vapply(cal, function(v) BJMconformal:::hpd_mass(y, f, v), 1), 0.1)
    BJMconformal:::hpd_mass(y, f, rt(1, 2)) <= q
  })
  native <- mean(abs(rt(1e4, 2)) <= qnorm(0.95))
  expect_lt(native, 0.86)
  expect_gt(mean(cover), 0.86)
})

test_that("coverage_lower_bound handles censoring", {
  cv <- BJMconformal:::coverage_lower_bound(L = c(1, 1, 1), time = c(2, 0.5, 0.5),
                                            status = c(0, 1, 0))
  expect_equal(cv$covered, c(TRUE, FALSE, NA))
})
