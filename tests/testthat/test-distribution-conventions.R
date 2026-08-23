# The simplex family must follow the base-R d/p/q/r contract. Each expectation
# below is paired with the base-R behaviour it mirrors, so a future change that
# drifts away from the convention fails here rather than surprising a user.

test_that("NA propagates as NA and NaN as NaN, like dnorm/dbeta", {
  expect_identical(dsimplex(c(0.3, NA, NaN, 0.5), mu = 0.4, phi = 1)[2:3],
                   c(NA_real_, NaN))
  expect_identical(dsimplex(0.3, mu = NA, phi = 1), NA_real_)
  expect_identical(dsimplex(0.3, mu = 0.4, phi = NA), NA_real_)
  expect_identical(psimplex(c(0.3, NA), mu = 0.4, phi = 1)[2], NA_real_)
  expect_identical(qsimplex(c(0.3, NA), mu = 0.4, phi = 1)[2], NA_real_)
  # Reference behaviour in base R.
  expect_identical(dnorm(NA), NA_real_)
  expect_identical(dbeta(NaN, 2, 2), NaN)
})

test_that("an out-of-domain parameter gives NaN with the canonical warning", {
  expect_warning(v <- dsimplex(0.5, mu = 2, phi = 1), "NaNs produced")
  expect_identical(v, NaN)
  expect_warning(v <- dsimplex(0.5, mu = 0.4, phi = -1), "NaNs produced")
  expect_identical(v, NaN)
  expect_warning(v <- psimplex(0.5, mu = 0, phi = 1), "NaNs produced")
  expect_identical(v, NaN)
  expect_warning(v <- qsimplex(0.5, mu = 0.4, phi = 0), "NaNs produced")
  expect_identical(v, NaN)
  # r*() warns with "NAs produced", exactly as rbeta() does.
  expect_warning(v <- rsimplex(2L, mu = -1, phi = 1), "NAs produced")
  expect_true(all(is.nan(v)))
})

test_that("a probability outside [0, 1] gives NaN, like qbeta", {
  expect_warning(expect_identical(qsimplex(-0.1, 0.4, 1), NaN), "NaNs produced")
  expect_warning(expect_identical(qsimplex(1.1, 0.4, 1), NaN), "NaNs produced")
  expect_warning(expect_identical(qsimplex(0.5, 0.4, 1, log.p = TRUE), NaN),
                 "NaNs produced")
  expect_identical(qsimplex(0, 0.4, 1), 0)
  expect_identical(qsimplex(1, 0.4, 1), 1)
})

test_that("x outside the support has density 0, and is not a warning", {
  expect_silent(v <- dsimplex(c(-1, 0, 1, 2), mu = 0.4, phi = 1))
  expect_identical(v, rep(0, 4L))
  expect_identical(dsimplex(c(0, 1), mu = 0.4, phi = 1, log = TRUE),
                   rep(-Inf, 2L))
  expect_identical(psimplex(c(-1, 0), mu = 0.4, phi = 1), c(0, 0))
  expect_identical(psimplex(c(1, 2), mu = 0.4, phi = 1), c(1, 1))
})

test_that("psimplex reproduces a direct numerical integration of the density", {
  for (par in list(c(0.35, 0.8), c(0.5, 2), c(0.15, 0.5), c(0.8, 3),
                   c(0.5, 0.01))) {
    mu <- par[1L]; phi <- par[2L]
    qs <- c(0.001, 0.05, 0.2, mu, 0.6, 0.9, 0.999)
    ref <- vapply(qs, function(u) {
      stats::integrate(function(t) dsimplex(t, mu, phi), 1e-12, u,
                       rel.tol = 1e-12, subdivisions = 2000L)$value
    }, numeric(1))
    expect_equal(psimplex(qs, mu, phi), ref, tolerance = 1e-9,
                 info = paste(mu, phi))
    expect_equal(psimplex(1, mu, phi), 1, tolerance = 1e-10)
  }
})

test_that("lower.tail and log.p behave as in base R", {
  x <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  expect_equal(psimplex(x, 0.4, 1) + psimplex(x, 0.4, 1, lower.tail = FALSE),
               rep(1, length(x)))
  expect_equal(psimplex(x, 0.4, 1, log.p = TRUE), log(psimplex(x, 0.4, 1)))
  p <- c(0.1, 0.3, 0.7)
  expect_equal(qsimplex(log(p), 0.4, 1, log.p = TRUE), qsimplex(p, 0.4, 1))
  expect_equal(qsimplex(p, 0.4, 1, lower.tail = FALSE), qsimplex(1 - p, 0.4, 1))
})

test_that("qsimplex inverts psimplex across the tails", {
  ps <- c(1e-6, 1e-3, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 0.999, 1 - 1e-6)
  for (par in list(c(0.35, 0.8), c(0.5, 2), c(0.15, 0.5), c(0.8, 3),
                   c(0.5, 0.01))) {
    ys <- qsimplex(ps, par[1L], par[2L])
    expect_true(all(ys > 0 & ys < 1))
    expect_equal(psimplex(ys, par[1L], par[2L]), ps, tolerance = 1e-8,
                 info = paste(par, collapse = " "))
  }
})

test_that("qsimplex agrees with the empirical quantiles of rsimplex", {
  set.seed(11L)
  s <- rsimplex(2e5L, 0.35, 0.8)
  ps <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  expect_equal(unname(stats::quantile(s, ps)), qsimplex(ps, 0.35, 0.8),
               tolerance = 5e-3)
})
