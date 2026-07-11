# Tests for fastsimplexregmixed() and the "simplex_fast_mixed" S3 methods.
# All fits use n_threads = 1L for determinism and modest sizes for speed.

sim_mixed <- function(J = 120L, nj = 9L, sigma = 0.6, seed = 20260711L, q = 1L) {
  set.seed(seed)
  n <- J * nj
  g <- factor(rep(seq_len(J), each = nj))
  x1 <- rnorm(n); z1 <- rnorm(n)
  if (q == 1L) {
    b <- rnorm(J, 0, sigma)[g]
    eta <- 0.4 - 0.7 * x1 + b
  } else {
    b0 <- rnorm(J, 0, sigma)[g]
    b1 <- rnorm(J, 0, 0.4)[g]
    eta <- 0.4 - 0.7 * x1 + b0 + b1 * x1
  }
  y <- rsimplex(n, simplex_linkinv(eta, "logit"), exp(-0.3 + 0.4 * z1))
  data.frame(g = g, x1 = x1, z1 = z1, y = y)
}

test_that("random-effects specification parsing rejects invalid inputs", {
  dat <- sim_mixed(J = 30L)
  expect_error(fastsimplexregmixed(y ~ x1, data = dat, n_threads = 1L), "'random' must be supplied")
  expect_error(fastsimplexregmixed(y ~ x1, random = y ~ 1 | g, data = dat, n_threads = 1L),
               "one-sided")
  expect_error(fastsimplexregmixed(y ~ x1, random = ~ 1, data = dat, n_threads = 1L),
               "grouping bar")
  expect_error(fastsimplexregmixed(y ~ x1, random = ~ 1 | g:x1, data = dat, n_threads = 1L),
               "single variable")
})

test_that("fastsimplexregmixed converges and recovers a random intercept", {
  dat <- sim_mixed(J = 150L, nj = 10L, sigma = 0.6)
  fit <- fastsimplexregmixed(y ~ x1 | z1, random = ~ 1 | g, data = dat,
                             nAGQ = 9L, n_threads = 1L)
  expect_s3_class(fit, "simplex_fast_mixed")
  expect_identical(fit$convergence, 0L)
  # Fixed effects near truth.
  expect_equal(unname(coef(fit, "mean")), c(0.4, -0.7), tolerance = 0.15)
  expect_equal(unname(coef(fit, "dispersion")), c(-0.3, 0.4), tolerance = 0.15)
  # Variance component near truth (sigma^2 = 0.36).
  expect_equal(as.numeric(VarCorr(fit)), 0.36, tolerance = 0.15)
  # True fixed effects inside Wald 99% intervals (compare by position).
  est <- fit$par[1:4]
  se <- fit$standard_errors[1:4]
  truth <- c(0.4, -0.7, -0.3, 0.4)
  inside <- truth >= est - 2.576 * se & truth <= est + 2.576 * se
  expect_true(all(inside))
})

test_that("the AGHQ marginal reduces to the fixed-effects fit as Sigma -> 0", {
  dat <- sim_mixed(J = 40L, nj = 8L)
  ns <- getNamespace("fastsimplexreg")
  # Build design matrices in group-contiguous order.
  X <- cbind(1, dat$x1); Z <- matrix(1, nrow(dat), 1); W <- matrix(1, nrow(dat), 1)
  gi <- as.integer(dat$g); ord <- order(gi)
  starts <- as.integer(c(0, cumsum(tabulate(gi, nlevels(dat$g)))))
  yo <- dat$y[ord]; Xo <- X[ord, , drop = FALSE]; Zo <- Z[ord, , drop = FALSE]; Wo <- W[ord, , drop = FALSE]
  beta <- c(0.3, -0.6); gamma <- -0.2
  fe_nll <- -sum(dsimplex(dat$y, simplex_linkinv(beta[1] + beta[2] * dat$x1, "logit"),
                          rep(exp(gamma), nrow(dat)), log = TRUE))
  th <- c(beta, gamma, -12)  # log-sd = -12 => sigma ~ 6e-6
  mix_nll <- ns$simplex_mixed_eval_cpp(th, yo, Xo, Zo, Wo, starts, 1L, 1L, 9L, 1L)$value
  expect_equal(mix_nll, fe_nll, tolerance = 1e-4)
})

test_that("analytic AGHQ gradient matches numerical differentiation", {
  skip_if_not_installed("numDeriv")
  dat <- sim_mixed(J = 30L, nj = 8L)
  ns <- getNamespace("fastsimplexreg")
  X <- cbind(1, dat$x1); Z <- matrix(1, nrow(dat), 1); W <- cbind(1, dat$z1)
  gi <- as.integer(dat$g); ord <- order(gi)
  starts <- as.integer(c(0, cumsum(tabulate(gi, nlevels(dat$g)))))
  yo <- dat$y[ord]; Xo <- X[ord, , drop = FALSE]; Zo <- Z[ord, , drop = FALSE]; Wo <- W[ord, , drop = FALSE]
  th <- c(0.3, -0.6, -0.2, 0.3, log(0.5))
  f <- function(p) ns$simplex_mixed_eval_cpp(p, yo, Xo, Zo, Wo, starts, 1L, 1L, 25L, 1L)$value
  ga <- as.numeric(ns$simplex_mixed_eval_cpp(th, yo, Xo, Zo, Wo, starts, 1L, 1L, 25L, 1L)$gradient)
  gn <- numDeriv::grad(f, th)
  expect_equal(ga, gn, tolerance = 1e-4)
})

test_that("a random intercept + slope model (q = 2) fits", {
  dat <- sim_mixed(J = 120L, nj = 12L, sigma = 0.6, q = 2L)
  fit <- fastsimplexregmixed(y ~ x1 | z1, random = ~ 1 + x1 | g, data = dat,
                             nAGQ = 5L, n_threads = 1L)
  expect_identical(fit$convergence, 0L)
  expect_equal(dim(fit$Sigma), c(2L, 2L))
  expect_equal(fit$q, 2L)
  expect_equal(dim(ranef(fit)), c(120L, 2L))
  # Sigma is a valid covariance (SPD).
  expect_true(all(eigen(fit$Sigma, only.values = TRUE)$values > 0))
})

test_that("S3 methods return the expected shapes", {
  dat <- sim_mixed(J = 80L, nj = 10L)
  fit <- fastsimplexregmixed(y ~ x1 | z1, random = ~ 1 | g, data = dat,
                             nAGQ = 7L, n_threads = 1L)

  expect_length(coef(fit, "all"), 4L)
  expect_length(coef(fit, "mean"), 2L)
  expect_equal(dim(vcov(fit)), c(5L, 5L))
  expect_equal(attr(logLik(fit), "df"), 5L)
  expect_identical(nobs(fit), 800L)
  expect_identical(ngrps(fit), 80L)

  re <- ranef(fit, postVar = TRUE)
  expect_equal(dim(re), c(80L, 1L))
  expect_equal(dim(attr(re, "postVar")), c(1L, 1L, 80L))

  vc <- VarCorr(fit)
  expect_s3_class(vc, "VarCorr.simplex_fast_mixed")
  expect_false(is.null(attr(vc, "stddev")))

  expect_length(fitted(fit), 800L)
  for (ty in c("response", "pearson", "deviance")) {
    expect_true(all(is.finite(residuals(fit, type = ty))))
  }
  # In-sample prediction equals the conditional fitted values.
  expect_equal(predict(fit, type = "response"), fitted(fit))
  # Population prediction (re.form = NA) differs.
  expect_false(isTRUE(all.equal(predict(fit, type = "response"),
                                predict(fit, type = "response", re.form = NA))))

  expect_output(print(fit), "simplex mixed model")
  sm <- summary(fit)
  expect_s3_class(sm, "summary.simplex_fast_mixed")
  expect_output(print(sm), "Random effects")
})

test_that("predict on new data handles known and unknown groups", {
  dat <- sim_mixed(J = 60L, nj = 8L)
  fit <- fastsimplexregmixed(y ~ x1 | z1, random = ~ 1 | g, data = dat,
                             nAGQ = 5L, n_threads = 1L)
  nd <- rbind(
    data.frame(g = factor("1", levels = levels(dat$g)), x1 = 0.5, z1 = 0.1),
    data.frame(g = factor(NA, levels = levels(dat$g)), x1 = 0.5, z1 = 0.1)
  )
  p <- predict(fit, newdata = nd, type = "response")
  expect_length(p, 2L)
  expect_true(all(p > 0 & p < 1))
})

test_that("input validation errors are raised", {
  dat <- sim_mixed(J = 30L)
  expect_error(fastsimplexregmixed(y ~ x1, random = ~ 1 | g, data = dat, nAGQ = 0L, n_threads = 1L),
               "positive integer")
  expect_error(fastsimplexregmixed(y ~ x1, random = ~ 1 | g, data = dat,
                                   start = c(0, 0), n_threads = 1L),
               "length")
  bad <- dat; bad$y[1] <- 1.5
  expect_error(fastsimplexregmixed(y ~ x1, random = ~ 1 | g, data = bad, n_threads = 1L),
               "strictly inside")
})

test_that("the AGHQ evaluation is invariant to the number of threads", {
  dat <- sim_mixed(J = 60L, nj = 8L)
  ns <- getNamespace("fastsimplexreg")
  X <- cbind(1, dat$x1); Z <- matrix(1, nrow(dat), 1); W <- cbind(1, dat$z1)
  gi <- as.integer(dat$g); ord <- order(gi)
  starts <- as.integer(c(0, cumsum(tabulate(gi, nlevels(dat$g)))))
  yo <- dat$y[ord]; Xo <- X[ord, , drop = FALSE]; Zo <- Z[ord, , drop = FALSE]; Wo <- W[ord, , drop = FALSE]
  th <- c(0.3, -0.6, -0.2, 0.3, log(0.5))
  r1 <- ns$simplex_mixed_eval_cpp(th, yo, Xo, Zo, Wo, starts, 1L, 1L, 11L, 1L)
  r2 <- ns$simplex_mixed_eval_cpp(th, yo, Xo, Zo, Wo, starts, 1L, 1L, 11L, 2L)
  expect_equal(r1$value, r2$value, tolerance = 1e-9)
  expect_equal(as.numeric(r1$gradient), as.numeric(r2$gradient), tolerance = 1e-9)
})
