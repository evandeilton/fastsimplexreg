# Tests for the estimator fastsimplexreg() and its S3 methods.
# Data are simulated from known coefficients (logit mean submodel + log
# dispersion submodel with a covariate). All fits use n_threads = 1L.

# Shared simulated data set with known truth.
make_data <- function(n = 2000L, seed = 20260710L) {
  set.seed(seed)
  dat <- data.frame(
    x1 = rnorm(n),
    x2 = rbinom(n, 1L, 0.4),
    z1 = rnorm(n)
  )
  beta_true  <- c(`(Intercept)` = -0.3, x1 = 0.8, x2 = -0.5)
  gamma_true <- c(`(Intercept)` = -0.5, z1 = 0.6)
  eta_mu   <- beta_true[1] + beta_true[2] * dat$x1 + beta_true[3] * dat$x2
  mu_true  <- simplex_linkinv(eta_mu, "logit")
  phi_true <- exp(gamma_true[1] + gamma_true[2] * dat$z1)
  dat$y <- rsimplex(n, mu = mu_true, phi = phi_true)
  list(
    dat = dat,
    # True parameter vector in coefficient order: mean coefficients first, then
    # dispersion coefficients (matching coef()/confint() row order).
    truth = unname(c(beta_true, gamma_true))
  )
}

test_that("fastsimplexreg converges and returns a simplex_fast object", {
  sim <- make_data()
  fit <- fastsimplexreg(y ~ x1 + x2 | z1, data = sim$dat, link = "logit",
                     n_threads = 1L)
  expect_s3_class(fit, "simplex_fast")
  expect_identical(fit$convergence, 0L)
})

test_that("true coefficients lie within Wald 99% confidence intervals", {
  sim <- make_data()
  fit <- fastsimplexreg(y ~ x1 + x2 | z1, data = sim$dat, link = "logit",
                     n_threads = 1L)
  ci <- confint(fit, level = 0.99)
  # Compare by position: coefficient names are bare and the two intercepts share
  # the name "(Intercept)".
  truth <- sim$truth
  inside <- truth >= ci[, 1] & truth <= ci[, 2]
  expect_true(all(inside),
              info = paste("Parameters outside 99% CI:",
                           paste(rownames(ci)[!inside], collapse = ", ")))
})

test_that("coef() honours the model argument", {
  sim <- make_data()
  fit <- fastsimplexreg(y ~ x1 + x2 | z1, data = sim$dat, link = "logit",
                     n_threads = 1L)
  all_coef  <- coef(fit, model = "all")
  mean_coef <- coef(fit, model = "mean")
  disp_coef <- coef(fit, model = "dispersion")
  expect_length(all_coef, 5L)
  expect_length(mean_coef, 3L)
  expect_length(disp_coef, 2L)
  expect_named(mean_coef, c("(Intercept)", "x1", "x2"))
  expect_named(disp_coef, c("(Intercept)", "z1"))
  expect_equal(unname(c(mean_coef, disp_coef)), unname(all_coef))
})

test_that("vcov / logLik / nobs / AIC / BIC are coherent", {
  sim <- make_data()
  fit <- fastsimplexreg(y ~ x1 + x2 | z1, data = sim$dat, link = "logit",
                     n_threads = 1L)
  vc <- vcov(fit)
  expect_true(is.matrix(vc))
  expect_equal(dim(vc), c(5L, 5L))
  expect_equal(vc, t(vc), tolerance = 1e-10)          # symmetric
  expect_true(all(diag(vc) > 0))                       # positive variances

  ll <- logLik(fit)
  expect_s3_class(ll, "logLik")
  expect_equal(attr(ll, "df"), 5L)
  expect_equal(attr(ll, "nobs"), nobs(fit))
  expect_identical(nobs(fit), 2000L)

  k <- 5L; n <- 2000L
  expect_equal(AIC(fit), -2 * as.numeric(ll) + 2 * k, tolerance = 1e-8)
  expect_equal(BIC(fit), -2 * as.numeric(ll) + log(n) * k, tolerance = 1e-8)
})

test_that("fitted / residuals behave as documented", {
  sim <- make_data()
  fit <- fastsimplexreg(y ~ x1 + x2 | z1, data = sim$dat, link = "logit",
                     n_threads = 1L)
  mu  <- fitted(fit, model = "mean")
  phi <- fitted(fit, model = "dispersion")
  expect_length(mu, 2000L)
  expect_true(all(mu > 0 & mu < 1))
  expect_true(all(phi > 0))

  r_resp <- residuals(fit, type = "response")
  expect_equal(r_resp, sim$dat$y - mu, tolerance = 1e-10)

  r_pear <- residuals(fit, type = "pearson")
  expect_length(r_pear, 2000L)
  expect_true(all(is.finite(r_pear)))
})

test_that("predict() supports all types and matches in-sample fitted", {
  sim <- make_data()
  fit <- fastsimplexreg(y ~ x1 + x2 | z1, data = sim$dat, link = "logit",
                     n_threads = 1L)
  expect_equal(predict(fit, type = "response"), fitted(fit, "mean"))
  expect_equal(predict(fit, type = "mean"), fitted(fit, "mean"))
  expect_equal(predict(fit, type = "dispersion"), fitted(fit, "dispersion"))

  lk <- predict(fit, type = "link")
  expect_named(lk, c("mean", "dispersion"))
  expect_length(lk$mean, 2000L)

  both <- predict(fit, type = "both")
  expect_s3_class(both, "data.frame")
  expect_named(both, c("mu", "phi"))

  # Prediction on identical newdata reproduces the in-sample fit.
  p_new <- predict(fit, newdata = sim$dat, type = "response")
  expect_equal(p_new, fitted(fit, "mean"), tolerance = 1e-10)
})

test_that("print and summary produce output without error", {
  sim <- make_data()
  fit <- fastsimplexreg(y ~ x1 + x2 | z1, data = sim$dat, link = "logit",
                     n_threads = 1L)
  expect_output(print(fit), "simplex regression")
  sm <- summary(fit)
  expect_s3_class(sm, "summary.simplex_fast")
  expect_type(sm$coefficients, "list")
  expect_named(sm$coefficients, c("mean", "dispersion"))
  expect_equal(ncol(sm$coefficients$mean), 4L)        # Estimate, SE, z, p
  expect_equal(ncol(sm$coefficients$dispersion), 4L)
  expect_output(print(sm), "Coefficients \\(mean model")
})

test_that("inference = FALSE yields NA std errors and vcov() errors", {
  sim <- make_data()
  fit <- fastsimplexreg(y ~ x1 + x2 | z1, data = sim$dat, link = "logit",
                     n_threads = 1L, inference = FALSE)
  expect_identical(fit$convergence, 0L)
  expect_true(all(is.na(fit$standard_errors)))
  expect_error(vcov(fit), "inference")
})

test_that("all four mean links fit and converge", {
  sim <- make_data(n = 1000L, seed = 555L)
  for (lk in c("logit", "probit", "cloglog", "neglog")) {
    fit <- fastsimplexreg(y ~ x1 + x2 | z1, data = sim$dat, link = lk,
                       n_threads = 1L)
    expect_identical(fit$convergence, 0L, info = lk)
    expect_equal(fit$link$mean, lk, info = lk)
  }
})





# The suite otherwise checks INTERNAL coherence -- that deviance matches the
# residuals, that AIC matches logLik. Nothing anchored the estimates outside the
# package. This test maximises an independently written pure-R likelihood with
# stats::optim and requires the two to agree.
test_that("the MLE matches an independent pure-R likelihood, for every link", {
  linkinvs <- list(
    logit = stats::plogis,
    probit = stats::pnorm,
    cloglog = function(e) -expm1(-exp(e)),
    neglog = function(e) exp(-exp(-e))
  )
  nll_R <- function(th, y, X, Z, linkinv) {
    p <- ncol(X)
    mu <- linkinv(drop(X %*% th[seq_len(p)]))
    ph <- exp(drop(Z %*% th[-seq_len(p)]))
    d <- (y - mu)^2 / (y * (1 - y) * mu^2 * (1 - mu)^2)
    -sum(-0.5 * (log(2 * pi) + log(ph)) - 1.5 * (log(y) + log(1 - y)) - 0.5 * d / ph)
  }

  set.seed(42)
  n <- 600L
  dat <- data.frame(x1 = rnorm(n), x2 = rbinom(n, 1L, 0.5), z1 = rnorm(n))
  for (lk in names(linkinvs)) {
    mu <- linkinvs[[lk]](0.2 + 0.6 * dat$x1 - 0.3 * dat$x2)
    dat$y <- rsimplex(n, mu, exp(-0.4 + 0.5 * dat$z1))
    fit <- fastsimplexreg(y ~ x1 + x2 | z1, data = dat, link = lk, n_threads = 1L)

    X <- cbind(1, dat$x1, dat$x2)
    Z <- cbind(1, dat$z1)
    opt <- stats::optim(coef(fit), nll_R, y = dat$y, X = X, Z = Z,
                        linkinv = linkinvs[[lk]], method = "BFGS",
                        control = list(reltol = 1e-14, maxit = 2000L))
    expect_equal(fit$logLik, -opt$value, tolerance = 1e-8,
                 info = paste("link:", lk))
    expect_equal(unname(coef(fit)), unname(opt$par), tolerance = 1e-5,
                 info = paste("link:", lk))
  }
})

test_that("subset selects the rows the fit is actually computed on", {
  set.seed(8)
  n <- 300L
  dat <- data.frame(x1 = rnorm(n), grp = rep(c("a", "b"), length.out = n))
  dat$y <- rsimplex(n, simplex_linkinv(0.2 + 0.7 * dat$x1, "logit"), 1)

  sub <- fastsimplexreg(y ~ x1, data = dat, subset = grp == "a", n_threads = 1L)
  manual <- fastsimplexreg(y ~ x1, data = dat[dat$grp == "a", ], n_threads = 1L)
  expect_equal(coef(sub), coef(manual))
  expect_equal(nobs(sub), sum(dat$grp == "a"))
})


# Badly scaled designs. Until 0.2.4 the first trial step was 1.0 in the units of
# the RAW gradient, whose magnitude is set by the units of the design. With a
# covariate scaled by 1e5 the gradient norm reaches ~1e8, a unit step overshoots
# by eight orders of magnitude, and the 40 available halvings bottom out at
# 2^-40 without reaching a decrease -- so the fit stopped at iteration 1 with
# code 2 and a log-likelihood 155.5 units below the optimum that nlminb()
# reaches on the identical objective.
test_that("a badly scaled design converges to the same optimum as nlminb", {
  mk <- function(scale, seed = 1L) {
    set.seed(seed)
    n <- 400L
    x <- rnorm(n) * scale
    d <- data.frame(x1 = x)
    d$y <- rsimplex(n, simplex_linkinv(0.3 + (0.5 / scale) * x, "logit"), 1)
    d
  }
  for (scale in c(1e2, 1e4, 1e5, 1e6, 1e8)) {
    d <- mk(scale)
    fit <- fastsimplexreg(y ~ x1, data = d, n_threads = 1L, inference = FALSE)
    X <- cbind(1, d$x1)
    Z <- matrix(1, nrow(d), 1L)
    obj <- function(th) fastsimplexreg:::simplex_eval_cpp(th, d$y, X, Z, 1L, 1L)$value
    gr <- function(th) fastsimplexreg:::simplex_eval_cpp(th, d$y, X, Z, 1L, 1L)$gradient
    ref <- stats::nlminb(c(0, 0, 0), obj, gr,
                         control = list(rel.tol = 1e-14, iter.max = 1000L))

    expect_identical(fit$convergence, 0L, info = format(scale))
    expect_equal(fit$logLik, -ref$objective, tolerance = 1e-6,
                 info = format(scale))
  }
})

test_that("scaling the first step does not cost evaluations on ordinary fits", {
  # The scaling is restricted to the state where H is still the identity.
  # Applying it on every iteration would be a serious regression; restricted, it
  # is cheaper, because the very first step no longer wastes halvings.
  total <- 0L
  for (s in 1:10) {
    set.seed(s)
    n <- 500L
    d <- data.frame(x1 = rnorm(n), x2 = rbinom(n, 1L, 0.4), z1 = rnorm(n))
    d$y <- rsimplex(n, simplex_linkinv(0.2 + 0.7 * d$x1 - 0.4 * d$x2, "logit"),
                    exp(-0.5 + 0.4 * d$z1))
    f <- fastsimplexreg(y ~ x1 + x2 | z1, data = d, n_threads = 1L,
                        inference = FALSE)
    expect_identical(f$convergence, 0L)
    total <- total + f$function_evaluations
  }
  # Measured at 1379 evaluations before the change and 1073 after, over 20 fits;
  # this half of the grid must stay comfortably under the old budget.
  expect_lt(total, 900L)
})
