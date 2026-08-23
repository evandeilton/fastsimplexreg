# Inference must never present a missing standard error as a confident one.
# These tests pin the fail-safe behaviour introduced in 0.2.3 (see R/inference.R).

test_that("a rank-deficient Hessian yields NA standard errors, not zeros", {
  set.seed(3L)
  n <- 300L
  dat <- data.frame(x1 = rnorm(n))
  dat$x2 <- dat$x1 + rnorm(n, sd = 1e-7)          # numerically collinear
  dat$y <- rsimplex(n, simplex_linkinv(0.3 + 0.5 * dat$x1, "logit"), 1)

  # Whether the degenerate direction lands just below zero ("not positive
  # definite") or just above it ("rank deficient") is a rounding accident; both
  # must warn and both must withhold the standard error.
  expect_warning(fit <- fastsimplexreg(y ~ x1 + x2, data = dat, n_threads = 1L),
                 "not positive definite|rank deficient")
  expect_identical(fit$convergence, 0L)
  expect_true(fit$vcov_pseudo)
  expect_lt(fit$vcov_rank, length(fit$par))

  se <- fit$standard_errors
  # The collinear pair is unidentified: NA, never 0.
  expect_true(all(is.na(se[c("x1", "x2")])))
  expect_false(any(se %in% 0))
  # The identified parameters keep usable standard errors.
  expect_true(is.finite(se[["(Intercept)"]]))

  tab <- summary(fit)$coefficients$mean
  expect_true(all(is.na(tab[c("x1", "x2"), "Pr(>|z|)"])))
  expect_true(all(is.na(confint(fit)["x1", ])))
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
