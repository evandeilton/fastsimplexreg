# Aliased columns in fastsimplexregmixed(): detection, NA reporting, reduced-model
# equivalence and the methods that consume an aliased fit.
# Fits are small (J = 24, nj = 6, nAGQ = 7, one thread) to keep the file fast.

sim_mixed_alias <- function(J = 24L, nj = 6L, seed = 20261016L) {
  set.seed(seed)
  n <- J * nj
  g <- factor(rep(seq_len(J), each = nj))
  x1 <- rnorm(n); z1 <- rnorm(n)
  b <- rnorm(J, 0, 0.6)[g]
  y <- rsimplex(n, simplex_linkinv(0.3 + 0.5 * x1 + b, "logit"),
                exp(-0.3 + 0.3 * z1))
  # x2 and z2 are exact multiples of x1 and z1, hence aliased when added.
  data.frame(g = g, x1 = x1, x2 = 2 * x1, z1 = z1, z2 = -3 * z1, y = y)
}

fit_mixed <- function(formula, data, ...) {
  fastsimplexregmixed(formula, random = ~ 1 | g, data = data, nAGQ = 7L,
                      n_threads = 1L, ...)
}

dat_a <- sim_mixed_alias()
fr <- fit_mixed(y ~ x1 | z1, dat_a)                    # reduced reference model

test_that("an aliased mean column is dropped and reported as NA", {
  expect_warning(fa <- fit_mixed(y ~ x1 + x2 | z1, dat_a, x = TRUE),
                 "rank deficient")
  expect_identical(fa$convergence, 0L)

  # Flags are named like `par` (fixed effects only, dispersion prefixed).
  expect_identical(names(fa$aliased), names(fa$par)[seq_along(fa$aliased)])
  expect_identical(names(fa$aliased),
                   c("(Intercept)", "x1", "x2", "(phi)_(Intercept)", "(phi)_z1"))
  expect_identical(unname(fa$aliased), c(FALSE, FALSE, TRUE, FALSE, FALSE))

  # NA, never 0, in every user-facing place; all else is finite.
  expect_identical(unname(which(is.na(fa$par))), 3L)
  expect_true(is.na(coef(fa, "mean")[["x2"]]))
  expect_true(is.na(fa$standard_errors[["x2"]]))
  est <- !is.na(fa$par)
  expect_true(all(is.finite(fa$standard_errors[est]) & fa$standard_errors[est] > 0))

  # The reduced design is what was fitted: same path as y ~ x1.
  expect_equal(unname(fa$par[est]), unname(fr$par), tolerance = 1e-6)
  expect_equal(fa$logLik, fr$logLik, tolerance = 1e-6)
  expect_equal(fa$AIC, fr$AIC, tolerance = 1e-6)
  expect_equal(fa$BIC, fr$BIC, tolerance = 1e-6)
  expect_equal(unname(fa$standard_errors[est]), unname(fr$standard_errors),
               tolerance = 1e-5)
  expect_equal(fitted(fa), fitted(fr), tolerance = 1e-6)
  expect_equal(residuals(fa, "response"), residuals(fr, "response"),
               tolerance = 1e-6)

  # df, vcov and hessian count the estimable parameters only.
  expect_identical(fa$df, fr$df)
  expect_identical(fa$df.residual, fr$df.residual)
  expect_identical(attr(logLik(fa), "df"), sum(est))
  expect_identical(dim(vcov(fa)), rep(sum(est), 2L))
  expect_identical(rownames(vcov(fa)), names(fr$par))
  expect_equal(vcov(fa), vcov(fr), tolerance = 1e-5)
  expect_identical(dimnames(fa$hessian), dimnames(vcov(fa)))

  # x = TRUE stores the FULL design, aliased column included.
  expect_identical(colnames(fa$x$mean), c("(Intercept)", "x1", "x2"))
  expect_identical(dim(fa$x$mean), c(nrow(dat_a), 3L))
})

test_that("predict, summary and confint survive an aliased fit", {
  fa <- suppressWarnings(fit_mixed(y ~ x1 + x2 | z1 + z2, dat_a))
  expect_identical(unname(fa$aliased),
                   c(FALSE, FALSE, TRUE, FALSE, FALSE, TRUE))
  expect_identical(names(fa$aliased)[fa$aliased], c("x2", "(phi)_z2"))
  expect_equal(fa$logLik, fr$logLik, tolerance = 1e-6)

  # In-sample prediction equals the fitted values; population level and
  # newdata agree with the reduced model, both drop the aliased columns.
  expect_equal(predict(fa), fitted(fa))
  expect_false(anyNA(predict(fa, re.form = NA)))
  expect_equal(predict(fa, re.form = NA), predict(fr, re.form = NA),
               tolerance = 1e-6)
  nd <- dat_a[1:7, ]
  expect_equal(predict(fa, newdata = nd), predict(fr, newdata = nd),
               tolerance = 1e-6)
  expect_equal(predict(fa, newdata = nd, re.form = NA),
               predict(fr, newdata = nd, re.form = NA), tolerance = 1e-6)
  both <- predict(fa, newdata = nd, type = "both")
  expect_false(anyNA(both))
  expect_equal(both$phi, predict(fr, newdata = nd, type = "dispersion"),
               tolerance = 1e-6)

  # summary(): NA rows for the aliased terms, a note in the printout.
  s <- summary(fa)
  expect_true(all(is.na(s$coefficients$mean["x2", ])))
  expect_true(all(is.na(s$coefficients$dispersion["z2", ])))
  expect_false(anyNA(s$coefficients$mean["x1", ]))
  expect_identical(s$npar, fa$df)
  out <- capture.output(print(s))
  expect_match(out, "Aliased (not estimable): x2, (phi)_z2", fixed = TRUE, all = FALSE)
  expect_false(any(grepl("Moore-Penrose", out)))
  expect_error(capture.output(print(fa)), NA)

  # confint(): positional over the full vector, NA rows for aliased terms.
  ci <- confint(fa)
  expect_identical(rownames(ci), names(fa$par))
  expect_true(all(is.na(ci["x2", ])))
  expect_true(all(is.na(ci["(phi)_z2", ])))
  expect_false(anyNA(ci[c("x1", "(phi)_z1", "logsd.(Intercept)"), ]))
  expect_true(all(is.na(confint(fa, parm = "x2"))))
})

test_that("a design with every column aliased is an error", {
  dat <- dat_a
  dat$zero <- 0
  expect_error(suppressWarnings(fit_mixed(y ~ 0 + zero, dat)),
               "Every column of the mean design is aliased")
  expect_error(suppressWarnings(fit_mixed(y ~ x1 | 0 + zero, dat)),
               "Every column of the dispersion design is aliased")
  expect_warning(try(fit_mixed(y ~ 0 + zero, dat), silent = TRUE),
                 "rank deficient")
})

test_that("start follows the reduced design", {
  # The length counts the non-aliased columns, as in fastsimplexreg().
  fa <- suppressWarnings(fit_mixed(y ~ x1 + x2 | z1, dat_a, start = fr$par,
                                   inference = FALSE))
  expect_equal(fa$logLik, fr$logLik, tolerance = 1e-6)
  expect_error(suppressWarnings(
    fit_mixed(y ~ x1 + x2 | z1, dat_a, start = c(fr$par, 0))), "'start'")
  # Without inference there are no standard errors, but the NA layout holds.
  expect_null(fa$vcov)
  expect_true(all(is.na(fa$standard_errors)))
  expect_identical(unname(which(is.na(fa$par))), 3L)
})

test_that("an unaliased mixed fit is unchanged and AIC/BIC match the fit", {
  expect_false(any(fr$aliased))
  expect_identical(names(fr$aliased), names(fr$par)[seq_along(fr$aliased)])
  expect_identical(fr$df, length(fr$par))
  expect_identical(attr(logLik(fr), "df"), length(fr$par))
  expect_equal(AIC(fr), fr$AIC)
  expect_equal(BIC(fr), fr$BIC)
  expect_identical(dim(vcov(fr)), rep(length(fr$par), 2L))
  expect_identical(summary(fr)$npar, length(fr$par))

  fa <- suppressWarnings(fit_mixed(y ~ x1 + x2 | z1, dat_a))
  expect_equal(AIC(fa), fa$AIC)
  expect_equal(BIC(fa), fa$BIC)
  expect_equal(AIC(fa), AIC(fr), tolerance = 1e-6)
})
