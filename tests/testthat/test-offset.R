# Offsets in the mean and dispersion submodels.
#
# Until 0.2.4 an `offset()` term in the formula was recorded by the `terms`
# object and then silently dropped from the design: the model that was fitted
# was not the model the user wrote, and nothing warned. These tests pin the
# offset to the only definition that makes it meaningful -- a column whose
# coefficient is fixed at one.

sim_offset <- function(n = 400L, seed = 99L) {
  set.seed(seed)
  d <- data.frame(x1 = stats::rnorm(n), off = stats::rnorm(n, 0, 0.8))
  d$y <- rsimplex(n, simplex_linkinv(0.2 + 0.7 * d$x1 + d$off, "logit"), 1)
  d
}

test_that("a mean offset enters the linear predictor with coefficient one", {
  d <- sim_offset()
  fit <- fastsimplexreg(y ~ x1 + offset(off), data = d, n_threads = 1L)
  th <- coef(fit)

  X <- cbind(`(Intercept)` = 1, x1 = d$x1)
  Z <- matrix(1, nrow(X), 1L)
  # Same theta, two routes: as an offset, and as an extra design column whose
  # coefficient is pinned to 1. They must agree exactly, value and gradient.
  via_offset <- fastsimplexreg:::simplex_eval_cpp(th, d$y, X, Z, 1L, 1L, d$off, NULL)
  via_column <- fastsimplexreg:::simplex_eval_cpp(c(th[1:2], 1, th[3]), d$y, cbind(X, d$off), Z, 1L, 1L)

  expect_equal(via_offset$value, via_column$value, tolerance = 1e-12)
  expect_equal(via_offset$gradient[1:2], via_column$gradient[1:2], tolerance = 1e-12)
})

test_that("the offset changes the fit rather than being ignored", {
  d <- sim_offset()
  with_off <- fastsimplexreg(y ~ x1 + offset(off), data = d, n_threads = 1L)
  without <- fastsimplexreg(y ~ x1, data = d, n_threads = 1L)

  expect_false(isTRUE(all.equal(unname(coef(with_off)), unname(coef(without)))))
  # The data were generated WITH the offset, so honouring it must fit better.
  expect_gt(with_off$logLik, without$logLik)
  # ... and must recover the generating coefficients.
  expect_equal(unname(coef(with_off)[1:2]), c(0.2, 0.7), tolerance = 0.1)
})

test_that("predict() rebuilds the offset from newdata", {
  d <- sim_offset()
  fit <- fastsimplexreg(y ~ x1 + offset(off), data = d, n_threads = 1L)
  expect_equal(unname(predict(fit, newdata = d[1:10, ])),
               unname(predict(fit)[1:10]))
  # Dropping the offset column from newdata must error, not silently zero it.
  expect_error(predict(fit, newdata = d[1:10, "x1", drop = FALSE]),
               "missing from 'newdata'")
})

test_that("a dispersion offset is kept separate from the mean offset", {
  set.seed(7)
  n <- 400L
  d <- data.frame(x1 = stats::rnorm(n), o1 = stats::rnorm(n, 0, 0.5),
                  o2 = stats::rnorm(n, 0, 0.5))
  d$y <- rsimplex(n, simplex_linkinv(0.3 + 0.6 * d$x1 + d$o1, "logit"),
                  exp(-0.2 + d$o2))
  fit <- fastsimplexreg(y ~ x1 + offset(o1) | offset(o2), data = d,
                        n_threads = 1L)
  # stats::model.offset() would return o1 + o2 for BOTH parts; the per-part
  # extraction must keep them apart.
  expect_equal(fit$offset$mean, d$o1)
  expect_equal(fit$offset$dispersion, d$o2)
})

test_that("an NA in the offset drops that row, like any other model variable", {
  d <- sim_offset(n = 50L)
  bad <- d
  bad$off[c(3L, 11L)] <- NA
  fit <- fastsimplexreg(y ~ x1 + offset(off), data = bad, n_threads = 1L)
  expect_equal(nobs(fit), 48L)
  expect_length(fit$offset$mean, 48L)
  expect_true(all(is.finite(fit$offset$mean)))
})

test_that("a non-finite offset that survives the model frame is refused", {
  d <- sim_offset(n = 50L)
  d$off[5] <- Inf
  expect_error(fastsimplexreg(y ~ x1 + offset(off), data = d, n_threads = 1L),
               "non-finite")
})

test_that("the mixed model honours a mean offset", {
  set.seed(5)
  J <- 40L; nj <- 8L; n <- J * nj
  d <- data.frame(g = factor(rep(seq_len(J), each = nj)),
                  x1 = stats::rnorm(n), off = stats::rnorm(n, 0, 0.6))
  b <- stats::rnorm(J, 0, 0.6)[d$g]
  d$y <- rsimplex(n, simplex_linkinv(0.3 - 0.5 * d$x1 + d$off + b, "logit"), 1)

  with_off <- fastsimplexregmixed(y ~ x1 + offset(off), random = ~ 1 | g,
                                  data = d, nAGQ = 7L, n_threads = 1L)
  without <- fastsimplexregmixed(y ~ x1, random = ~ 1 | g, data = d,
                                 nAGQ = 7L, n_threads = 1L)
  expect_gt(with_off$logLik, without$logLik)
  expect_equal(with_off$offset$mean, d$off)
  expect_equal(unname(predict(with_off, newdata = d[1:5, ])),
               unname(predict(with_off)[1:5]))
})

test_that("a model without an offset is numerically unchanged", {
  d <- sim_offset()
  fit <- fastsimplexreg(y ~ x1, data = d, n_threads = 1L)
  expect_null(fit$offset$mean)
  expect_null(fit$offset$dispersion)
  # The backend must treat "no offset" as an exact no-op, not as a zero vector
  # that has been added.
  X <- cbind(1, d$x1); Z <- matrix(1, nrow(X), 1L)
  th <- coef(fit)
  expect_equal(fastsimplexreg:::simplex_eval_cpp(th, d$y, X, Z, 1L, 1L)$value,
               fastsimplexreg:::simplex_eval_cpp(th, d$y, X, Z, 1L, 1L, rep(0, nrow(X)), NULL)$value,
               tolerance = 1e-14)
})
