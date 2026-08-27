# Inference must never present a missing standard error as a confident one.
# These tests pin the fail-safe behaviour introduced in 0.2.3 (see R/inference.R).

test_that("an exactly collinear column is dropped and reported as NA", {
  set.seed(3L)
  n <- 300L
  dat <- data.frame(x1 = rnorm(n))
  dat$x2 <- dat$x1                       # exact duplicate
  dat$y <- rsimplex(n, simplex_linkinv(0.3 + 0.5 * dat$x1, "logit"), 1)

  # From 0.2.4 the rank deficiency is caught in the DESIGN, before fitting, by
  # the same pivoted QR lm() uses -- rather than being left to surface as a
  # singular Hessian afterwards. Fitting the full design returned a finite
  # estimate for both columns, an arbitrary split of the one identified effect.
  expect_warning(fit <- fastsimplexreg(y ~ x1 + x2, data = dat, n_threads = 1L),
                 "rank deficient")
  expect_identical(fit$convergence, 0L)
  expect_true(fit$aliased[["x2"]])
  expect_false(fit$aliased[["x1"]])

  # The aliased column is NA in every user-facing place, never 0 and never a
  # share of the identified effect.
  expect_true(is.na(coef(fit)[["x2"]]))
  expect_true(is.na(fit$standard_errors[["x2"]]))
  expect_true(all(is.na(confint(fit)["x2", ])))
  expect_true(all(is.na(summary(fit)$coefficients$mean["x2", ])))

  # The surviving column carries the WHOLE identified effect, and now has a
  # usable standard error because the design it was fitted on is full rank.
  reduced <- fastsimplexreg(y ~ x1, data = dat, n_threads = 1L)
  expect_equal(coef(fit)[["x1"]], coef(reduced)[["x1"]], tolerance = 1e-6)
  expect_equal(fit$logLik, reduced$logLik, tolerance = 1e-8)
  expect_true(is.finite(fit$standard_errors[["x1"]]))
  expect_false(any(fit$standard_errors %in% 0))

  # Downstream methods survive the NA.
  expect_length(fitted(fit), n)
  expect_equal(unname(predict(fit, newdata = dat[1:5, ])),
               unname(predict(fit)[1:5]))
})

test_that(".simplex_vcov withholds the standard error of a singular direction", {
  # The design-level guard above means a singular observed information no longer
  # arises from exact collinearity -- but the fail-safe inversion must still
  # never turn a missing standard error into a confident one, so it is exercised
  # directly here.
  H <- matrix(c(4, 2, 2, 1), 2, 2)        # exactly singular: rank 1 of 2
  inf <- suppressWarnings(
    fastsimplexreg:::.simplex_vcov(H, c("a", "b"), what = "test"))
  expect_lt(inf$rank, 2L)
  expect_true(inf$pseudo)
  expect_true(all(is.na(inf$se)))
  expect_false(any(inf$se %in% 0))
  expect_warning(fastsimplexreg:::.simplex_vcov(H, c("a", "b"), what = "test"),
                 "rank deficient")

  # A direction of NEGATIVE curvature is not a variance either.
  Hneg <- matrix(c(4, 0, 0, -1), 2, 2)
  inf2 <- suppressWarnings(
    fastsimplexreg:::.simplex_vcov(Hneg, c("a", "b"), what = "test"))
  expect_true(is.na(inf2$se[["b"]]))
  expect_warning(fastsimplexreg:::.simplex_vcov(Hneg, c("a", "b"), what = "test"),
                 "not positive definite")
})

test_that("weak identification gives large standard errors, not NA", {
  # The fail-safe must distinguish "the data cannot identify this at all" (NA)
  # from "the data identify it poorly" (a large, honest number). Only the first
  # justifies withholding the standard error.
  set.seed(3L)
  n <- 300L
  dat <- data.frame(x1 = rnorm(n))
  dat$x2 <- dat$x1 + rnorm(n, sd = 1e-3)
  dat$y <- rsimplex(n, simplex_linkinv(0.3 + 0.5 * dat$x1, "logit"), 1)

  fit <- suppressWarnings(fastsimplexreg(y ~ x1 + x2, data = dat, n_threads = 1L))
  expect_identical(fit$vcov_rank, length(fit$par))
  expect_false(fit$vcov_pseudo)
  expect_true(all(is.finite(fit$standard_errors)))
  # Large, but a real number.
  expect_gt(fit$standard_errors[["x1"]], 1)
})

test_that("a well-conditioned fit is untouched by the fail-safe path", {
  set.seed(9L)
  n <- 800L
  dat <- data.frame(x1 = rnorm(n), z1 = rnorm(n))
  dat$y <- rsimplex(n, simplex_linkinv(0.2 + 0.6 * dat$x1, "logit"),
                    exp(-0.5 + 0.3 * dat$z1))
  fit <- fastsimplexreg(y ~ x1 | z1, data = dat, n_threads = 1L)

  expect_false(fit$vcov_pseudo)
  expect_identical(fit$vcov_rank, length(fit$par))
  expect_true(all(is.finite(fit$standard_errors) & fit$standard_errors > 0))
  # Every mean link stays comfortably clear of the rank cut-off, including
  # neglog, whose raw Hessian spans 13 orders of magnitude by scale alone.
  for (lk in c("logit", "probit", "cloglog", "neglog")) {
    f <- fastsimplexreg(y ~ x1 | z1, data = dat, link = lk, n_threads = 1L)
    expect_identical(f$vcov_rank, length(f$par), info = lk)
    expect_false(f$vcov_pseudo, info = lk)
  }
  # Identical to the plain inverse whenever the plain inverse is legitimate.
  expect_equal(unname(vcov(fit)), unname(solve(fit$hessian)), tolerance = 1e-9)
  expect_identical(fit$n_saturated, 0L)
})

test_that("saturation of the mean link is reported rather than applied silently", {
  set.seed(2L)
  n <- 200L
  dat <- data.frame(x1 = rnorm(n))
  dat$y <- rsimplex(n, simplex_linkinv(0.3 + 0.5 * dat$x1, "logit"), 1)
  # A start far into the flat region of the link saturates every observation:
  # their contribution to the mean score is exactly zero by construction, which
  # used to happen silently.
  msgs <- character(0)
  fit <- withCallingHandlers(
    fastsimplexreg(y ~ x1, data = dat, start = c(-80, 0, 0), maxit = 1L,
                   n_threads = 1L),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_identical(fit$n_saturated, 200L)
  expect_true(any(grepl("saturated at the numerical boundary", msgs)))
})

test_that("the unclamped reporting path keeps its dynamic range", {
  # The likelihood path floors mu at 1e-12; the reporting path must not.
  expect_equal(simplex_linkinv(-40, "logit"), stats::plogis(-40),
               tolerance = 1e-12)
  expect_lt(simplex_linkinv(-40, "logit"), 1e-15)
  # ... while still staying strictly inside the open support (0, 1), which the
  # simplex density and every residual formula require.
  eta <- c(-800, -100, -40, 0, 40, 100, 800)
  for (lk in c("logit", "probit", "cloglog", "neglog")) {
    mu <- simplex_linkinv(eta, lk)
    expect_true(all(mu > 0 & mu < 1), info = lk)
  }
})
