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

test_that("deviance is informative by default and scaled on request", {
  obj <- make_fit()
  fit <- obj$fit
  rdev <- residuals(fit, type = "deviance")

  # The SCALED deviance is the sum of squared deviance residuals -- and is
  # identically nobs whenever the dispersion submodel has an intercept, because
  # that is exactly what the score equation for gamma forces. It therefore
  # cannot distinguish two models, which is why it is no longer the default.
  expect_equal(deviance(fit, type = "scaled"), sum(rdev^2))
  expect_equal(deviance(fit, type = "scaled"), as.numeric(nobs(fit)),
               tolerance = 1e-6)

  # The default is the UNSCALED deviance, which does respond to the mean model.
  expect_equal(deviance(fit), sum(rdev^2 * fitted(fit, "dispersion")))
  expect_false(isTRUE(all.equal(deviance(fit), as.numeric(nobs(fit)),
                                tolerance = 1e-3)))

  dat <- obj$dat
  bigger <- fastsimplexreg(y ~ x1 + x2 + I(x1^2) | z1, data = dat,
                           n_threads = 1L)
  expect_false(isTRUE(all.equal(deviance(fit), deviance(bigger))))
  # ... where the scaled version would have been identical for both.
  expect_equal(deviance(fit, type = "scaled"),
               deviance(bigger, type = "scaled"), tolerance = 1e-5)
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

test_that("full-vector names are unambiguous and summary tables stay bare", {
  fit <- make_fit()$fit

  # The FULL parameter vector must be addressable by name. Before 0.2.4 the
  # dispersion block reused the mean's bare names, so with an intercept in both
  # submodels names(coef(fit)) repeated "(Intercept)": vcov()["(Intercept)",
  # "(Intercept)"] silently returned the MEAN intercept's variance whatever the
  # user meant, and confint(parm = "(Intercept)") returned two identical-looking
  # rows. betareg solves this the same way, with a "(phi)_" prefix.
  all_names <- names(coef(fit, "all"))
  expect_false(any(grepl("^mean_|^dispersion_", all_names)))
  expect_identical(all_names,
                   c("(Intercept)", "x1", "x2", "(phi)_(Intercept)", "(phi)_z1"))
  expect_false(anyDuplicated(all_names) > 0L)
  expect_identical(rownames(vcov(fit)), all_names)
  expect_identical(rownames(confint(fit)), all_names)

  # Name-based access now reaches the parameter it names.
  expect_equal(unname(vcov(fit)["(phi)_(Intercept)", "(phi)_(Intercept)"]),
               unname(fit$standard_errors[["(phi)_(Intercept)"]]^2))
  expect_identical(nrow(confint(fit, parm = "(phi)_(Intercept)")), 1L)

  s <- summary(fit)
  # The per-submodel tables keep the BARE names: within a table there is no
  # ambiguity, and the prefix would only add noise.
  expect_type(s$coefficients, "list")
  expect_named(s$coefficients, c("mean", "dispersion"))
  expect_identical(rownames(s$coefficients$mean), c("(Intercept)", "x1", "x2"))
  expect_identical(rownames(s$coefficients$dispersion), c("(Intercept)", "z1"))
  expect_identical(colnames(s$coefficients$mean),
                   c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
  # coef(model = ) likewise stays bare.
  expect_identical(names(coef(fit, "dispersion")), c("(Intercept)", "z1"))

  # Positional selection keeps working.
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


test_that("predict keeps the length of newdata when rows carry NA", {
  fit <- make_fit()$fit
  nd <- data.frame(x1 = c(0.1, NA, 0.3), x2 = c(0, 1, 1), z1 = c(0, 0.2, NA))
  p <- predict(fit, newdata = nd, type = "response")
  # A shorter, unaligned vector was the old behaviour; predict.lm keeps length.
  expect_length(p, 3L)
  expect_equal(unname(which(is.na(p))), c(2L, 3L))

  both <- predict(fit, newdata = nd, type = "both")
  expect_equal(nrow(both), 3L)

  # A variable the model needs but newdata lacks is an error, not a guess.
  expect_error(predict(fit, newdata = data.frame(x1 = 0.1)), "missing from 'newdata'")
})

test_that("simulate carries the row names of the model frame", {
  set.seed(5L)
  d <- data.frame(x1 = rnorm(50))
  rownames(d) <- paste0("obs", seq_len(50))
  d$y <- rsimplex(50, simplex_linkinv(0.3 + 0.5 * d$x1, "logit"), 1)
  fit <- fastsimplexreg(y ~ x1, data = d, n_threads = 1L, model = TRUE)
  s <- simulate(fit, nsim = 2L, seed = 1L)
  expect_identical(rownames(s), rownames(d))
})


# na.action = na.exclude was accepted and had no effect: fitted() and
# residuals() came back with the number of COMPLETE rows, so nothing could be
# aligned back to the source data without knowing which rows had been dropped.
test_that("na.exclude pads fitted, residuals and predict back to nrow(data)", {
  set.seed(11)
  n <- 300L
  dat <- data.frame(x1 = rnorm(n), z1 = rnorm(n))
  dat$y <- rsimplex(n, simplex_linkinv(0.3 + 0.8 * dat$x1, "logit"),
                    exp(-0.6 + 0.4 * dat$z1))
  dat$x1[c(3L, 7L, 10L)] <- NA

  omit <- fastsimplexreg(y ~ x1 | z1, data = dat, n_threads = 1L,
                         na.action = stats::na.omit)
  excl <- fastsimplexreg(y ~ x1 | z1, data = dat, n_threads = 1L,
                         na.action = stats::na.exclude)

  expect_equal(nobs(omit), nobs(excl))
  # na.omit is unchanged: complete rows only.
  expect_length(fitted(omit), 297L)
  expect_length(residuals(omit), 297L)
  # na.exclude pads, exactly as glm does.
  expect_length(fitted(excl), n)
  expect_length(residuals(excl), n)
  expect_length(predict(excl), n)
  expect_equal(unname(which(is.na(fitted(excl)))), c(3L, 7L, 10L))
  expect_equal(unname(which(is.na(residuals(excl, "deviance")))), c(3L, 7L, 10L))
  expect_equal(length(fitted(excl)),
               length(stats::fitted(stats::glm(stats::qlogis(y) ~ x1, data = dat,
                                               na.action = stats::na.exclude))))

  # The padded and unpadded values agree on the complete rows.
  expect_equal(unname(fitted(excl)[-c(3L, 7L, 10L)]), unname(fitted(omit)))

  # Internal consumers must stay on the UNPADDED vectors or their lengths would
  # no longer match fitted.values.
  expect_silent(invisible(summary(excl)))
  skip_if_not_installed("ggplot2")
  expect_s3_class(plot(excl, which = 1L), "ggplot")
})

test_that("fitted, residuals and predict carry the observation labels", {
  set.seed(4)
  n <- 40L
  dat <- data.frame(x1 = rnorm(n))
  dat$y <- rsimplex(n, simplex_linkinv(0.2 + 0.6 * dat$x1, "logit"), 1)
  rownames(dat) <- paste0("obs", seq_len(n))

  fit <- fastsimplexreg(y ~ x1, data = dat, n_threads = 1L)
  expect_identical(names(fitted(fit)), rownames(dat))
  expect_identical(names(residuals(fit)), rownames(dat))
  expect_identical(names(predict(fit, newdata = dat[5:8, ])),
                   rownames(dat)[5:8])
})


# Randomized quantile residuals (Dunn and Smyth, 1996). Under a correct model
# these are EXACTLY standard normal, which neither Pearson nor deviance
# residuals are once the fitted means move away from 1/2.
test_that("quantile residuals are standard normal under a correct model", {
  set.seed(505)
  rate <- function(b0, phi0, B = 60L, n = 300L) {
    rej <- vapply(seq_len(B), function(b) {
      d <- data.frame(x1 = rnorm(n))
      d$y <- rsimplex(n, simplex_linkinv(b0 + 0.8 * d$x1, "logit"), phi0)
      f <- suppressWarnings(fastsimplexreg(y ~ x1, data = d, n_threads = 1L,
                                           inference = FALSE))
      c(quantile = stats::shapiro.test(residuals(f, "quantile"))$p.value,
        pearson = stats::shapiro.test(residuals(f, "pearson"))$p.value) < 0.05
    }, logical(2))
    rowMeans(rej)
  }

  # Means near 1/2: both behave.
  centred <- rate(0.0, 1)
  expect_lt(centred[["quantile"]], 0.20)
  # Means away from 1/2: Pearson rejects almost every CORRECT model, quantile
  # residuals hold their nominal rate.
  skewed <- rate(-1.4, 1)
  expect_lt(skewed[["quantile"]], 0.20)
  expect_gt(skewed[["pearson"]], 0.70)
})

test_that("quantile residuals are finite, named and the plot default", {
  obj <- make_fit()
  fit <- obj$fit
  r <- residuals(fit)                       # default type
  expect_equal(r, residuals(fit, "quantile"))
  expect_true(all(is.finite(r)))
  expect_identical(names(r), names(fitted(fit)))
  # Standard normal: mean ~ 0, sd ~ 1 on a correct model.
  expect_lt(abs(mean(r)), 0.15)
  expect_equal(stats::sd(r), 1, tolerance = 0.1)

  expect_identical(
    eval(formals(fastsimplexreg:::residuals.simplex_fast)$type)[1L], "quantile")
  expect_identical(
    eval(formals(fastsimplexreg:::plot.simplex_fast)$type)[1L], "quantile")
})
