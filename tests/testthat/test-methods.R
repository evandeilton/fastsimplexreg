# Tests for the standard S3 method surface of "simplex_fast" objects, beyond the
# estimator checks in test-simplexreg.R. Every fit uses n_threads = 1L.

make_fit <- function(n = 500L, store = TRUE) {
  set.seed(101L)
  dat <- data.frame(x1 = rnorm(n), x2 = rbinom(n, 1L, 0.4), z1 = rnorm(n))
  mu <- simplex_linkinv(-0.3 + 0.9 * dat$x1 - 0.5 * dat$x2, "logit")
  phi <- exp(-1.0 + 0.6 * dat$z1)
  dat$y <- rsimplex(n, mu, phi)
  fit <- fastsimplexreg(y ~ x1 + x2 | z1, data = dat, link = "logit",
                     n_threads = 1L, x = store, model = store)
  list(fit = fit, dat = dat)
}

test_that("structural accessors return the expected objects", {
  obj <- make_fit()
  fit <- obj$fit

  expect_s3_class(formula(fit), "formula")
  expect_identical(deparse(formula(fit)), "y ~ x1 + x2 | z1")

  expect_identical(attr(terms(fit, "mean"), "term.labels"), c("x1", "x2"))
  expect_identical(attr(terms(fit, "dispersion"), "term.labels"), "z1")

  mm_mean <- model.matrix(fit, "mean")
  mm_disp <- model.matrix(fit, "dispersion")
  expect_equal(dim(mm_mean), c(500L, 3L))
  expect_equal(dim(mm_disp), c(500L, 2L))
  expect_identical(colnames(mm_mean), c("(Intercept)", "x1", "x2"))

  expect_s3_class(model.frame(fit), "data.frame")
  expect_equal(nrow(model.frame(fit)), 500L)
})

test_that("model.matrix rebuilds from the model frame when x is not stored", {
  set.seed(7L)
  n <- 300L
  dat <- data.frame(x1 = rnorm(n))
  dat$y <- rsimplex(n, simplex_linkinv(0.3 + 0.5 * dat$x1, "logit"), 1)

  # No design stored, only the model frame.
  fit_mf <- fastsimplexreg(y ~ x1, data = dat, n_threads = 1L, x = FALSE, model = TRUE)
  expect_equal(dim(model.matrix(fit_mf, "mean")), c(300L, 2L))
  # Constant dispersion -> intercept-only dispersion design.
  expect_equal(dim(model.matrix(fit_mf, "dispersion")), c(300L, 1L))

  # Nothing stored -> informative error.
  fit_none <- fastsimplexreg(y ~ x1, data = dat, n_threads = 1L, x = FALSE, model = FALSE)
  expect_error(model.matrix(fit_none, "mean"), "Refit with")
})

test_that("deviance equals the sum of squared deviance residuals", {
  fit <- make_fit()$fit
  rdev <- residuals(fit, type = "deviance")
  expect_equal(deviance(fit), sum(rdev^2))
  expect_true(is.finite(deviance(fit)))
})

test_that("all residual types are finite and correctly signed", {
  fit <- make_fit()$fit
  y <- fitted(fit) + residuals(fit, "response")
  for (ty in c("response", "pearson", "deviance")) {
    r <- residuals(fit, type = ty)
    expect_length(r, nobs(fit))
    expect_true(all(is.finite(r)))
    # Every residual type shares the sign of (y - mu).
    expect_identical(sign(r), sign(y - fitted(fit)))
  }
})

test_that("simulate is reproducible, respects the support and records the seed", {
  fit <- make_fit()$fit
  s1 <- simulate(fit, nsim = 3L, seed = 42L)
  s2 <- simulate(fit, nsim = 3L, seed = 42L)
  expect_identical(s1, s2)
  expect_equal(dim(s1), c(nobs(fit), 3L))
  expect_identical(names(s1), c("sim_1", "sim_2", "sim_3"))
  vals <- unlist(s1, use.names = FALSE)
  expect_true(all(vals > 0 & vals < 1))
  expect_false(is.null(attr(s1, "seed")))
  expect_error(simulate(fit, nsim = 0L), "positive integer")
})

test_that("update refits with a modified multi-part formula", {
  # Build data in this scope so update() can resolve `data = dat` in its caller
  # frame, exactly as base-R update() requires.
  set.seed(101L)
  n <- 500L
  dat <- data.frame(x1 = rnorm(n), x2 = rbinom(n, 1L, 0.4), z1 = rnorm(n))
  mu <- simplex_linkinv(-0.3 + 0.9 * dat$x1 - 0.5 * dat$x2, "logit")
  dat$y <- rsimplex(n, mu, exp(-1.0 + 0.6 * dat$z1))
  fit <- fastsimplexreg(y ~ x1 + x2 | z1, data = dat, link = "logit", n_threads = 1L)

  fit2 <- update(fit, . ~ . - x2 | z1)
  expect_false("x2" %in% names(coef(fit2, "mean")))
  expect_identical(attr(terms(fit2, "mean"), "term.labels"), "x1")
  # Dispersion part preserved.
  expect_identical(attr(terms(fit2, "dispersion"), "term.labels"), "z1")
  # evaluate = FALSE returns an unevaluated call.
  expect_type(update(fit, . ~ . - x2, evaluate = FALSE), "language")
})

test_that("AIC, BIC and logLik are mutually consistent", {
  fit <- make_fit()$fit
  ll <- logLik(fit)
  k <- attr(ll, "df")
  expect_equal(attr(ll, "nobs"), nobs(fit))
  expect_equal(AIC(fit), -2 * as.numeric(ll) + 2 * k)
  expect_equal(BIC(fit), -2 * as.numeric(ll) + log(nobs(fit)) * k)
})

test_that("coefficient names are bare and summary tables match the standard", {
  fit <- make_fit()$fit

  # No "mean_"/"dispersion_" prefixes anywhere in the user-facing names.
  all_names <- names(coef(fit, "all"))
  expect_false(any(grepl("^mean_|^dispersion_", all_names)))
  expect_identical(all_names, c("(Intercept)", "x1", "x2", "(Intercept)", "z1"))
  expect_identical(rownames(vcov(fit)), all_names)
  expect_identical(rownames(confint(fit)), all_names)

  s <- summary(fit)
  # summary$coefficients is a list of mean/dispersion tables with bare rownames.
  expect_type(s$coefficients, "list")
  expect_named(s$coefficients, c("mean", "dispersion"))
  expect_identical(rownames(s$coefficients$mean), c("(Intercept)", "x1", "x2"))
  expect_identical(rownames(s$coefficients$dispersion), c("(Intercept)", "z1"))
  expect_identical(colnames(s$coefficients$mean),
                   c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))

  # confint selects the correct row by position even with duplicated names.
  ci_disp_int <- confint(fit, parm = 4L)
  expect_equal(unname(ci_disp_int[1, ]),
               unname(coef(fit, "all")[4] +
                        c(-1, 1) * qnorm(0.975) * fit$standard_errors[4]))
})

test_that("plot returns ggplot objects for single and multiple panels", {
  skip_if_not_installed("ggplot2")
  fit <- make_fit()$fit
  p_single <- plot(fit, which = 2L)
  expect_s3_class(p_single, "ggplot")

  p_multi <- plot(fit, which = 1:4)
  # A patchwork object when available, otherwise a list of ggplots.
  expect_true(inherits(p_multi, "patchwork") || is.list(p_multi))

  expect_error(plot(fit, which = 5L), "subset of 1:4")
})
