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
  # An unseen (here: missing) group level falls back to a zero random effect,
  # and that substitution is now announced instead of being silent.
  expect_warning(p <- predict(fit, newdata = nd, type = "response"),
                 "not seen in the fit")
  expect_length(p, 2L)
  expect_true(all(p > 0 & p < 1))
})

test_that("conditional prediction refuses newdata without the grouping column", {
  dat <- sim_mixed(J = 40L, nj = 8L)
  fit <- fastsimplexregmixed(y ~ x1 | z1, random = ~ 1 | g, data = dat,
                             nAGQ = 5L, n_threads = 1L, inference = FALSE)
  nd_ok <- data.frame(g = factor("1", levels = levels(dat$g)), x1 = 0.5, z1 = 0.1)
  nd_no_g <- data.frame(x1 = 0.5, z1 = 0.1)

  # Silently returning a population-level value here would be a wrong answer
  # under the label of a conditional one.
  expect_error(predict(fit, newdata = nd_no_g), "missing from 'newdata'")
  # The population level is still reachable explicitly, and then g is not needed.
  expect_length(predict(fit, newdata = nd_no_g, re.form = NA), 1L)
  expect_false(isTRUE(all.equal(predict(fit, newdata = nd_ok),
                                predict(fit, newdata = nd_ok, re.form = NA))))
})

test_that("population prediction needs a stored design, and says so", {
  dat <- sim_mixed(J = 60L, nj = 8L)
  fit <- fastsimplexregmixed(y ~ x1 | z1, random = ~ 1 | g, data = dat,
                             nAGQ = 5L, n_threads = 1L, inference = FALSE,
                             model = FALSE, x = FALSE)
  # Without the model frame the design would otherwise be rebuilt from whatever
  # objects happen to be visible in the caller.
  expect_error(predict(fit, re.form = NA), "model = TRUE")
})

test_that("predict keeps the length of newdata when rows carry NA", {
  dat <- sim_mixed(J = 40L, nj = 8L)
  fit <- fastsimplexregmixed(y ~ x1 | z1, random = ~ 1 | g, data = dat,
                             nAGQ = 5L, n_threads = 1L, inference = FALSE)
  nd <- data.frame(g = factor(c("1", "2", "3"), levels = levels(dat$g)),
                   x1 = c(0.5, NA, 0.3), z1 = c(0.1, 0.2, NA))
  p <- predict(fit, newdata = nd, type = "response")
  expect_length(p, 3L)
  expect_equal(which(is.na(p)), c(2L, 3L))
})

test_that("VarCorr refuses a meaningless sigma instead of ignoring it", {
  dat <- sim_mixed(J = 60L, nj = 8L)
  # Convergence is not what this test is about; this small configuration is a
  # known hard case for the outer optimiser (unchanged from 0.2.2).
  fit <- suppressWarnings(
    fastsimplexregmixed(y ~ x1, random = ~ 1 | g, data = dat,
                        nAGQ = 5L, n_threads = 1L, inference = FALSE))
  expect_silent(v1 <- VarCorr(fit))
  expect_warning(v2 <- VarCorr(fit, sigma = 10), "no meaning")
  expect_equal(as.numeric(v1), as.numeric(v2))
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

test_that("malformed cluster offsets error cleanly instead of crashing", {
  ns <- getNamespace("fastsimplexreg")
  y <- runif(12, 0.1, 0.9); X <- matrix(1, 12, 1); Z <- matrix(1, 12, 1); W <- matrix(1, 12, 1)
  th <- rep(0, 3)
  expect_error(ns$simplex_mixed_eval_cpp(th, y, X, Z, W, as.integer(c(0, 5, 5, 12)), 1L),
               "strictly increasing")
  expect_error(ns$simplex_mixed_eval_cpp(th, y, X, Z, W, as.integer(c(0, 8, 4, 12)), 1L),
               "strictly increasing")
  expect_error(ns$simplex_mixed_eval_cpp(th, y, X, Z, W, integer(0), 1L),
               "0, ..., nrow", fixed = TRUE)
})

test_that("non-convergence is signalled and standard errors are withheld", {
  # This used to rely on cloglog happening to miss the gradient tolerance on
  # this data -- an accident that the 0.2.3 stopping rule removed. Drive the
  # failure deliberately instead: one iteration from a start far outside the
  # sensible region cannot take a single step, which is the case the stopping
  # rule must still refuse to call convergence.
  dat <- sim_mixed(J = 40L, nj = 6L, seed = 7L)
  # Such a start also saturates the mean link, which warns on its own; that is
  # correct but incidental here.
  expect_warning(
    fit <- withCallingHandlers(
      fastsimplexregmixed(y ~ x1 | z1, random = ~ 1 | g, data = dat,
                          nAGQ = 5L, n_threads = 1L, maxit = 1L,
                          start = c(-80, 0, 0, 0, 0)),
      warning = function(w) if (grepl("saturated", conditionMessage(w)))
        invokeRestart("muffleWarning")),
    "did not converge")
  expect_identical(fit$convergence, 2L)
  expect_true(all(is.na(fit$standard_errors)))
  expect_output(print(summary(fit)), "DID NOT CONVERGE")
})

test_that("the nAGQ^q node budget is enforced", {
  set.seed(1); n <- 200
  dat <- data.frame(g = factor(rep(1:20, each = 10)), x1 = rnorm(n), z1 = rnorm(n))
  dat$y <- rsimplex(n, simplex_linkinv(0.2 + 0.3 * dat$x1, "logit"), 1)
  expect_error(
    fastsimplexregmixed(y ~ x1, random = ~ 1 + x1 + z1 | g, data = dat, nAGQ = 60, n_threads = 1),
    "nAGQ")
})

test_that("inference = FALSE skips standard errors in the mixed model", {
  dat <- sim_mixed(J = 60L, nj = 8L)
  fit <- fastsimplexregmixed(y ~ x1 | z1, random = ~ 1 | g, data = dat,
                             nAGQ = 7L, n_threads = 1L, inference = FALSE)
  expect_true(all(is.na(fit$standard_errors)))
  expect_error(vcov(fit), "inference")
})

test_that("unbalanced clusters (including singletons) are handled", {
  set.seed(3)
  sizes <- sample(1:6, 80, replace = TRUE)          # includes singletons
  g <- factor(rep(seq_along(sizes), sizes)); n <- length(g)
  x1 <- rnorm(n)
  b <- rnorm(nlevels(g), 0, 0.5)[g]
  dat <- data.frame(g = g, x1 = x1,
                    y = rsimplex(n, simplex_linkinv(0.3 - 0.5 * x1 + b, "logit"), 1))
  # The fit must run and return well-formed output even if it does not fully
  # converge on this small, ragged data (it warns in that case).
  fit <- suppressWarnings(
    fastsimplexregmixed(y ~ x1, random = ~ 1 | g, data = dat, nAGQ = 7, n_threads = 1))
  expect_length(fitted(fit), n)
  expect_equal(ngrps(fit), nlevels(g))
})


# Convergence reporting at a modest nAGQ. Before 0.2.3 the optimiser terminated
# with code 2 ("line search failed") on ~26% of mixed fits over a grid of
# 240 (J, nj, nAGQ, seed) combinations -- and, because standard errors are only
# computed at a converged fit, those perfectly good fits silently lost their
# inference. Restarting from the reported stopping point gained ~1e-11 in
# log-likelihood: it was an optimum, not a failure.
#
# The cause is that the analytic score is the exact score of the TRUE marginal
# likelihood (Fisher's identity), not of its nAGQ-point quadrature
# approximation, so `grad_tol` is unreachable when nAGQ is small.
test_that("a modest nAGQ converges instead of reporting a spurious failure", {
  dat <- sim_mixed(J = 30L, nj = 4L, seed = 7L)
  fit <- fastsimplexregmixed(y ~ x1 | z1, random = ~ 1 | g, data = dat,
                             nAGQ = 5L, n_threads = 1L)
  expect_identical(fit$convergence, 0L)
  expect_match(fit$message, "floor|tolerance satisfied")
  # The point of converging is that inference is no longer withheld.
  expect_true(all(is.finite(fit$standard_errors)))
  expect_false(is.null(fit$vcov))

  # And it is a real optimum, not a stall: nudged away and refitted, the
  # optimiser comes back to the same log-likelihood.
  nudged <- suppressWarnings(
    fastsimplexregmixed(y ~ x1 | z1, random = ~ 1 | g, data = dat,
                        nAGQ = 5L, n_threads = 1L, inference = FALSE,
                        start = unname(fit$par) + 0.05))
  expect_equal(nudged$logLik, fit$logLik, tolerance = 1e-4)
})

test_that("a genuine optimiser failure is still reported as one", {
  # Softening the stopping rule must not turn every stall into a success: a run
  # that cannot take a single step keeps code 2.
  set.seed(2L)
  n <- 200L
  dat <- data.frame(x1 = rnorm(n))
  dat$y <- rsimplex(n, simplex_linkinv(0.3 + 0.5 * dat$x1, "logit"), 1)
  fit <- suppressWarnings(
    fastsimplexreg(y ~ x1, data = dat, start = c(-80, 0, 0), maxit = 1L,
                   n_threads = 1L))
  expect_identical(fit$convergence, 2L)
  expect_match(fit$message, "Line search failed")
  expect_true(all(is.na(fit$standard_errors)))
})


# Adaptive Gauss-Hermite node pruning. Until 0.2.4 build_tensor() dropped nodes
# whose product weight logW sat below a floor -- but the adaptive transformation
# undoes the e^{-t^2} factor, so a node's real multiplier is logW + t2. Pruning
# on logW alone discarded 48.7% of the effective quadrature mass at nAGQ = 21,
# and the AGHQ sequence stopped converging: raising nAGQ moved the marginal
# log-likelihood AWAY from its limit.
test_that("the AGHQ sequence converges as nAGQ grows (q = 2)", {
  set.seed(77)
  J <- 25L; nj <- 3L; n <- J * nj
  g <- rep(seq_len(J), each = nj)
  X <- cbind(1, stats::rnorm(n)); Z <- X; W <- matrix(1, n, 1L)
  B <- matrix(stats::rnorm(J * 2L), J, 2L) %*% diag(c(2, 1.5))
  mu <- simplex_linkinv(X %*% c(0.3, -0.5) + rowSums(Z * B[g, ]), "logit")
  y <- rsimplex(n, as.numeric(mu), rep(1, n))
  starts <- as.integer(c(0L, cumsum(rep(nj, J))))
  th <- c(0.3, -0.5, 0, log(2), 0, log(1.5))

  ll <- vapply(c(11L, 15L, 21L, 25L), function(m) {
    -fastsimplexreg:::simplex_mixed_eval_cpp(th, y, X, Z, W, starts, 2L, 1L, m,
                                             1L, 50L, 1e-8)$value
  }, numeric(1))

  # Successive increments must shrink towards zero. With the old pruning they
  # grew instead, reaching ~1e-3 in this regime.
  steps <- abs(diff(ll))
  expect_lt(steps[length(steps)], 1e-6)
  expect_true(all(diff(steps) < 0))
})

test_that("nAGQ below 5 warns that the standard errors are unreliable", {
  dat <- sim_mixed(J = 20L, nj = 5L, seed = 2L)
  expect_warning(
    fastsimplexregmixed(y ~ x1, random = ~ 1 | g, data = dat, nAGQ = 1L,
                        n_threads = 1L, inference = FALSE),
    "below the supported minimum")
  # nAGQ = 2.7 trips BOTH guards: it is truncated to 2, and 2 is below 5.
  w <- testthat::capture_warnings(
    fastsimplexregmixed(y ~ x1, random = ~ 1 | g, data = dat, nAGQ = 2.7,
                        n_threads = 1L, inference = FALSE))
  expect_match(w, "truncated", all = FALSE)
  expect_match(w, "below the supported minimum", all = FALSE)
})

test_that("inner_maxit below 10 is refused", {
  dat <- sim_mixed(J = 20L, nj = 5L, seed = 2L)
  # With inner_maxit in 1..3 the inner solver need not reach the posterior mode,
  # so the AGHQ expansion is taken around the wrong point and the marginal
  # likelihood is grossly wrong -- previously with no diagnostic at all.
  expect_error(
    fastsimplexregmixed(y ~ x1, random = ~ 1 | g, data = dat, nAGQ = 7L,
                        inner_maxit = 3L, n_threads = 1L),
    "at least 10")
})

test_that("degenerate cluster structures are refused or flagged", {
  dat <- sim_mixed(J = 20L, nj = 5L, seed = 2L)

  one_level <- dat
  one_level$g <- factor(rep("a", nrow(dat)))
  expect_error(
    fastsimplexregmixed(y ~ x1, random = ~ 1 | g, data = one_level, nAGQ = 7L,
                        n_threads = 1L),
    "at least 2 sampled levels")

  singletons <- dat
  singletons$g <- factor(seq_len(nrow(dat)))
  expect_warning(
    fastsimplexregmixed(y ~ x1, random = ~ 1 | g, data = singletons, nAGQ = 7L,
                        n_threads = 1L, inference = FALSE),
    "single observation")
})


# confint() on a mixed fit. Until 0.2.4 there was no method, so dispatch fell
# through to stats::confint.default, which indexes vcov() by NAME -- and the
# repeated "(Intercept)" made it report the MEAN intercept's interval for the
# DISPERSION intercept, an interval that need not contain its own estimate.
# The variance components were dropped entirely, coef() being shorter than par.
test_that("confint on a mixed fit reports each parameter's own interval", {
  dat <- sim_mixed(J = 40L, nj = 10L, seed = 21L)
  fit <- fastsimplexregmixed(y ~ x1, random = ~ 1 | g, data = dat, nAGQ = 7L,
                             n_threads = 1L)
  ci <- confint(fit)
  est <- fit$par
  se <- fit$standard_errors

  expect_identical(nrow(ci), length(est))
  expect_identical(rownames(ci), names(est))
  expect_true(all(est >= ci[, 1] & est <= ci[, 2]))
  expect_equal(unname(ci),
               unname(cbind(est - qnorm(0.975) * se, est + qnorm(0.975) * se)))

  # Selection by position and by the now-unambiguous names.
  expect_identical(nrow(confint(fit, parm = 2L)), 1L)
  expect_identical(nrow(confint(fit, parm = "(phi)_(Intercept)")), 1L)
  expect_error(confint(fit, level = 1.5), "strictly between 0 and 1")
})

test_that("the packed omega diagonal is labelled as a Cholesky factor when q >= 2", {
  skip_if_not_installed("MASS")
  set.seed(31)
  J <- 60L; nj <- 10L; n <- J * nj
  d <- data.frame(g = factor(rep(seq_len(J), each = nj)), x1 = stats::rnorm(n))
  B <- MASS::mvrnorm(J, c(0, 0), matrix(c(0.6, 0.45, 0.45, 0.5), 2, 2))
  d$y <- rsimplex(n, simplex_linkinv(0.3 - 0.6 * d$x1 + B[d$g, 1] +
                                       B[d$g, 2] * d$x1, "logit"), 1)
  fit <- fastsimplexregmixed(y ~ x1, random = ~ 1 + x1 | g, data = d,
                             nAGQ = 5L, n_threads = 1L, inference = FALSE)

  # The label must not promise a marginal standard deviation it is not: for
  # j >= 2, exp(omega_jj) is D[j, j], the CONDITIONAL sd. In a measured fit it
  # understated the marginal sd of the random slope by a factor of ~1.95.
  expect_false(any(grepl("^logsd\\.", names(fit$omega))))
  expect_match(names(fit$omega)[3], "^logchol\\.")
  expect_equal(unname(exp(fit$omega[3])), unname(fit$D[2, 2]))
  expect_gt(sqrt(diag(fit$Sigma))[2], exp(fit$omega[3]))

  # q = 1 keeps "logsd.", where it is exactly right.
  f1 <- fastsimplexregmixed(y ~ x1, random = ~ 1 | g, data = d, nAGQ = 5L,
                            n_threads = 1L, inference = FALSE)
  expect_match(names(f1$omega)[1], "^logsd\\.")
  expect_equal(unname(exp(f1$omega[1])), unname(sqrt(f1$Sigma[1, 1])))
})
