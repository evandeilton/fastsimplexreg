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

# NOTE: the length contract of the distribution functions changed in 0.2.3.
# They used to raise an error unless mu/phi had length 1 or length(x); they now
# recycle to the common maximum length, which is the base-R d/p/q/r convention
# (dnorm(1:3, mean = 1:2) recycles silently). These tests encode the new
# contract deliberately -- see NEWS.md.
test_that("mu and phi are recycled to the common length, base-R style", {
  d <- dsimplex(c(0.2, 0.5, 0.8), mu = c(0.4, 0.6), phi = 1, n_threads = 1L)
  expect_length(d, 3L)
  # Recycling is positional: mu = 0.4, 0.6, 0.4.
  expect_equal(d, dsimplex(c(0.2, 0.5, 0.8), mu = c(0.4, 0.6, 0.4), phi = 1))

  expect_length(dsimplex(c(0.2, 0.5, 0.8), mu = 0.5, phi = c(1, 2)), 3L)
  # The result takes the maximum length of the three arguments.
  expect_length(dsimplex(0.5, mu = c(0.2, 0.4, 0.6), phi = 1), 3L)
  expect_length(psimplex(0.5, mu = c(0.2, 0.4, 0.6), phi = 1), 3L)
  expect_length(qsimplex(0.5, mu = c(0.2, 0.4, 0.6), phi = 1), 3L)

  set.seed(1)
  expect_length(rsimplex(3L, mu = c(0.4, 0.6), phi = 1), 3L)
  expect_length(rsimplex(3L, mu = 0.5, phi = c(1, 2)), 3L)
})

test_that("a zero-length argument gives a zero-length result", {
  expect_length(dsimplex(numeric(0), mu = 0.5, phi = 1), 0L)
  expect_length(dsimplex(0.5, mu = numeric(0), phi = 1), 0L)
  expect_length(psimplex(numeric(0), mu = 0.5, phi = 1), 0L)
})

test_that("rsimplex follows the base-R 'length(n) > 1' convention", {
  set.seed(1)
  expect_length(rsimplex(c(10, 20, 30), mu = 0.5, phi = 1), 3L)
  expect_length(rsimplex(0L, mu = 0.5, phi = 1), 0L)
  expect_error(rsimplex(-1L, mu = 0.5, phi = 1), "non-negative")
})

test_that("fastsimplexreg rejects a start vector of the wrong length", {
  dat <- data.frame(y = runif(20, 0.1, 0.9), x1 = rnorm(20), z1 = rnorm(20))
  # Model needs ncol(X) + ncol(Z) = 2 + 1 = 3 parameters; supply 2.
  expect_error(
    fastsimplexreg(y ~ x1 | z1, data = dat, start = c(0, 0), n_threads = 1L),
    "length ncol"
  )
})


# `subset` is non-standard-evaluated, as in stats::lm/glm. Until 0.2.4 it was a
# plain index vector: `subset = x1 > 0` silently used an unrelated `x1` from the
# caller when one existed, fitting the wrong rows with no symptom at all.
test_that("subset is evaluated inside data, not in the caller", {
  set.seed(1)
  n <- 200L
  dat <- data.frame(x1 = stats::rnorm(n))
  dat$y <- rsimplex(n, simplex_linkinv(0.2 + 0.7 * dat$x1, "logit"), 1)
  # A DIFFERENT vector of the same name in the calling frame. The old code used
  # this one; the correct behaviour is to prefer the data column.
  x1 <- stats::rnorm(n)
  stopifnot(sum(dat$x1 > 0) != sum(x1 > 0))   # the two must actually differ

  fit <- fastsimplexreg(y ~ x1, data = dat, subset = x1 > 0, n_threads = 1L)
  expect_equal(nobs(fit), sum(dat$x1 > 0))
  expect_equal(nobs(fit), nobs(stats::glm(stats::qlogis(y) ~ x1, data = dat,
                                          subset = x1 > 0)))
})

test_that("every pre-0.2.4 subset form still works", {
  set.seed(1)
  n <- 200L
  dat <- data.frame(x1 = stats::rnorm(n))
  dat$y <- rsimplex(n, simplex_linkinv(0.2 + 0.7 * dat$x1, "logit"), 1)
  f <- function(s) nobs(fastsimplexreg(y ~ x1, data = dat, subset = s,
                                       n_threads = 1L))

  expect_equal(f(1:50), 50L)                       # numeric index
  expect_equal(f(-(1:20)), 180L)                   # negative index
  expect_equal(f(dat$x1 > 0), sum(dat$x1 > 0))     # logical
  expect_equal(nobs(fastsimplexreg(y ~ x1, data = dat, n_threads = 1L)), n)

  lg <- dat$x1 > 0
  lg[1:5] <- NA                                    # NA drops the row, as in glm
  expect_equal(f(lg), nobs(stats::glm(stats::qlogis(y) ~ x1, data = dat,
                                      subset = lg)))
})

test_that("a logical subset of the wrong length is refused", {
  set.seed(1)
  n <- 50L
  dat <- data.frame(x1 = stats::rnorm(n))
  dat$y <- rsimplex(n, simplex_linkinv(0.2 * dat$x1, "logit"), 1)
  expect_error(fastsimplexreg(y ~ x1, data = dat, subset = c(TRUE, FALSE),
                              n_threads = 1L),
               "length nrow\\(data\\)")
})
