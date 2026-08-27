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


# The CDF is a closed form from 0.2.4 on, not adaptive quadrature. These tests
# pin it against an INDEPENDENT route and against the tail behaviour that the
# quadrature could not deliver.
test_that("the closed-form CDF matches seeded numerical integration", {
  ref <- function(q, mu, phi) {
    # Knots seeded around the mean: a naive integrate() misses the needle when
    # phi is small and returns 6.2e-163 for a probability of 1.
    sd <- sqrt(phi * (mu * (1 - mu))^3)
    kn <- sort(unique(pmax(1e-14, pmin(q, mu + c(-12, -6, -3, -1, 0, 1, 3, 6, 12) * sd))))
    kn <- c(1e-14, kn[kn > 1e-14 & kn < q], q)
    sum(vapply(seq_len(length(kn) - 1L), function(i)
      stats::integrate(function(t) dsimplex(t, mu, phi), kn[i], kn[i + 1L],
                       rel.tol = 1e-13, subdivisions = 2000L)$value, numeric(1)))
  }
  for (mu in c(0.05, 0.5, 0.95)) {
    for (phi in c(1e-3, 1, 100)) {
      for (q in c(0.15, 0.5, 0.85)) {
        r <- ref(q, mu, phi)
        if (r < 1e-12) next
        expect_equal(psimplex(q, mu, phi), r, tolerance = 1e-9)
      }
    }
  }
})

test_that("log.p keeps full precision in the far lower tail", {
  # psimplex(0.15, 0.5, 0.01, log.p = TRUE) used to return -Inf, because it
  # computed log() of a linear value that had already underflowed to zero.
  expect_equal(psimplex(0.15, 0.5, 0.01, log.p = TRUE), -773.2159, tolerance = 1e-4)
  expect_true(is.finite(psimplex(0.4, 0.5, 1e-4, log.p = TRUE)))
  expect_lt(psimplex(0.4, 0.5, 1e-4, log.p = TRUE), -3000)
  # Consistent with the linear scale wherever the linear scale still works.
  q <- c(0.05, 0.2, 0.5, 0.8, 0.95)
  expect_equal(psimplex(q, 0.4, 1, log.p = TRUE), log(psimplex(q, 0.4, 1)),
               tolerance = 1e-12)
  expect_equal(psimplex(q, 0.4, 1, lower.tail = FALSE, log.p = TRUE),
               log(psimplex(q, 0.4, 1, lower.tail = FALSE)), tolerance = 1e-12)
})

test_that("both tails keep full relative accuracy", {
  # P(Y > y | mu) is computed by the exact reflection F(1-y | 1-mu), not as
  # 1 - F(y), so the upper tail does not lose precision to cancellation.
  for (mu in c(0.1, 0.5, 0.9)) {
    for (phi in c(1e-3, 1, 100)) {
      q <- seq(0.05, 0.95, by = 0.15)
      expect_equal(psimplex(q, mu, phi, lower.tail = FALSE),
                   psimplex(1 - q, 1 - mu, phi), tolerance = 1e-12)
      expect_equal(psimplex(q, mu, phi) + psimplex(q, mu, phi, lower.tail = FALSE),
                   rep(1, length(q)), tolerance = 1e-12)
    }
  }
})
