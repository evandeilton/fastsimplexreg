# Input-validation tests: malformed inputs must raise informative errors,
# exercised entirely through the public R API.

test_that("fastsimplexreg rejects responses outside (0, 1)", {
  dat <- data.frame(y = c(0.2, 1.5, 0.4), x1 = rnorm(3))
  expect_error(
    fastsimplexreg(y ~ x1, data = dat, n_threads = 1L),
    "strictly inside"
  )
  dat0 <- data.frame(y = c(0.2, 0, 0.4), x1 = rnorm(3))
  expect_error(
    fastsimplexreg(y ~ x1, data = dat0, n_threads = 1L),
    "strictly inside"
  )
})

test_that("fastsimplexreg rejects a formula with two response parts", {
  dat <- data.frame(y1 = runif(5, 0.1, 0.9), y2 = runif(5, 0.1, 0.9),
                    x1 = rnorm(5))
  expect_error(
    fastsimplexreg(y1 | y2 ~ x1, data = dat, n_threads = 1L),
    "one response"
  )
})

test_that("fastsimplexreg rejects a formula with three RHS parts", {
  dat <- data.frame(y = runif(6, 0.1, 0.9), x1 = rnorm(6),
                    x2 = rnorm(6), x3 = rnorm(6))
  expect_error(
    fastsimplexreg(y ~ x1 | x2 | x3, data = dat, n_threads = 1L),
    "one or two RHS"
  )
})

test_that("dsimplex rejects mismatched mu/phi lengths", {
  expect_error(
    dsimplex(c(0.2, 0.5, 0.8), mu = c(0.4, 0.6), phi = 1, n_threads = 1L),
    "length 1 or length"
  )
  expect_error(
    dsimplex(c(0.2, 0.5, 0.8), mu = 0.5, phi = c(1, 2), n_threads = 1L),
    "length 1 or length"
  )
})

test_that("rsimplex rejects mismatched mu/phi lengths", {
  expect_error(
    rsimplex(3L, mu = c(0.4, 0.6), phi = 1),
    "length 1 or n"
  )
  expect_error(
    rsimplex(3L, mu = 0.5, phi = c(1, 2)),
    "length 1 or n"
  )
})

test_that("fastsimplexreg rejects a start vector of the wrong length", {
  dat <- data.frame(y = runif(20, 0.1, 0.9), x1 = rnorm(20), z1 = rnorm(20))
  # Model needs ncol(X) + ncol(Z) = 2 + 1 = 3 parameters; supply 2.
  expect_error(
    fastsimplexreg(y ~ x1 | z1, data = dat, start = c(0, 0), n_threads = 1L),
    "length ncol"
  )
})
