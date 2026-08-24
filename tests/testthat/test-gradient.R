# End-to-end validation of the model / analytic-score consistency.
#
# We reconstruct the negative log-likelihood (NLL) purely in R from the public
# dsimplex() density and simplex_linkinv(), then verify:
#   (a) the R NLL at the fitted optimum equals -logLik(fit) (the compiled
#       likelihood and the R reconstruction agree);
#   (b) a finite-difference gradient of the R NLL is essentially zero at the
#       fitted optimum (the analytic score used by the optimiser is consistent
#       with the likelihood surface);
#   (c) at a perturbed parameter, two independent finite-difference gradients
#       (different step sizes) agree, confirming the reconstruction is smooth
#       and well-behaved.
#
# No C++ symbol is referenced directly. n is kept small (300) for speed.

# Reconstruct the NLL(theta) = -sum log f(y; mu(theta), phi(theta)) in pure R.
make_nll <- function(y, X, Z, link) {
  p <- ncol(X)
  q <- ncol(Z)
  function(theta) {
    beta  <- theta[seq_len(p)]
    gamma <- theta[p + seq_len(q)]
    eta_mu  <- as.numeric(X %*% beta)
    mu      <- simplex_linkinv(eta_mu, link)
    eta_phi <- as.numeric(Z %*% gamma)
    phi     <- exp(eta_phi)
    -sum(dsimplex(y, mu = mu, phi = phi, log = TRUE, n_threads = 1L))
  }
}

# Self-coded central-difference gradient (no external dependency).
central_grad <- function(fn, theta, rel_step = 1e-5) {
  g <- numeric(length(theta))
  for (j in seq_along(theta)) {
    h <- rel_step * max(1, abs(theta[j]))
    tp <- theta; tp[j] <- tp[j] + h
    tm <- theta; tm[j] <- tm[j] - h
    g[j] <- (fn(tp) - fn(tm)) / (2 * h)
  }
  g
}

simulate_fit <- function(n = 300L, seed = 42L) {
  set.seed(seed)
  dat <- data.frame(
    x1 = rnorm(n),
    x2 = rbinom(n, 1L, 0.4),
    z1 = rnorm(n)
  )
  eta_mu   <- -0.3 + 0.8 * dat$x1 - 0.5 * dat$x2
  mu_true  <- simplex_linkinv(eta_mu, "logit")
  phi_true <- exp(-0.5 + 0.4 * dat$z1)
  dat$y <- rsimplex(n, mu = mu_true, phi = phi_true)
  fit <- fastsimplexreg(y ~ x1 + x2 | z1, data = dat, link = "logit",
                     n_threads = 1L)
  list(fit = fit, dat = dat)
}

test_that("R-reconstructed NLL matches -logLik at the optimum", {
  sim <- simulate_fit()
  fit <- sim$fit; dat <- sim$dat
  X <- model.matrix(~ x1 + x2, data = dat)
  Z <- model.matrix(~ z1, data = dat)
  nll <- make_nll(dat$y, X, Z, "logit")
  theta_hat <- coef(fit, model = "all")
  expect_equal(nll(theta_hat), -as.numeric(logLik(fit)), tolerance = 1e-6)
})

test_that("finite-difference gradient of the R NLL is ~0 at the optimum", {
  sim <- simulate_fit()
  fit <- sim$fit; dat <- sim$dat
  X <- model.matrix(~ x1 + x2, data = dat)
  Z <- model.matrix(~ z1, data = dat)
  nll <- make_nll(dat$y, X, Z, "logit")
  theta_hat <- coef(fit, model = "all")
  g <- central_grad(nll, theta_hat)
  # The optimiser converges on the analytic score; the FD gradient of the R
  # NLL therefore inherits a comparably small infinity norm. A generous
  # threshold guards against finite-difference noise on n = 300 data.
  expect_lt(max(abs(g)), 1e-2)
})

test_that("two independent FD gradients agree at a perturbed theta", {
  sim <- simulate_fit()
  fit <- sim$fit; dat <- sim$dat
  X <- model.matrix(~ x1 + x2, data = dat)
  Z <- model.matrix(~ z1, data = dat)
  nll <- make_nll(dat$y, X, Z, "logit")
  theta <- coef(fit, model = "all") + 0.1  # move away from the optimum
  g1 <- central_grad(nll, theta, rel_step = 1e-4)
  g2 <- central_grad(nll, theta, rel_step = 1e-6)
  # Both central-difference schemes approximate the same true gradient to
  # O(h^2); they should agree closely away from the optimum.
  expect_equal(g1, g2, tolerance = 1e-4)
  # And the gradient is genuinely non-zero away from the optimum.
  expect_gt(max(abs(g1)), 1e-3)
})


# Direct check of the ANALYTIC score against numerical differentiation, for all
# four mean links. The tests above validate the optimiser end to end but never
# compare the analytic gradient with a numerical one; a wrong analytic score
# typically still converges, just to the wrong point, so it deserves its own
# direct test.
test_that("the analytic score matches numerical differentiation for every link", {
  skip_if_not_installed("numDeriv")
  ns <- asNamespace("fastsimplexreg")
  set.seed(7L)
  n <- 400L
  X <- cbind(1, rnorm(n), rbinom(n, 1L, 0.4))
  Z <- cbind(1, rnorm(n))
  links <- c(logit = 1L, probit = 2L, cloglog = 3L, neglog = 4L)

  for (lk in names(links)) {
    set.seed(11L)
    mu <- simplex_linkinv(as.numeric(X %*% c(-0.3, 0.8, -0.5)), lk)
    phi <- exp(as.numeric(Z %*% c(-0.5, 0.4)))
    y <- rsimplex(n, mu, phi)
    theta <- c(-0.2, 0.6, -0.3, -0.4, 0.3)   # deliberately away from the optimum

    f <- function(t) ns$simplex_eval_cpp(t, y, X, Z, links[[lk]], 1L)$value
    analytic <- as.numeric(ns$simplex_eval_cpp(theta, y, X, Z, links[[lk]], 1L)$gradient)
    numeric_g <- numDeriv::grad(f, theta)
    # Compare on a relative scale: the neglog gradient is O(1e10) at this theta.
    expect_lt(max(abs(analytic - numeric_g)) / max(abs(numeric_g)), 1e-7)
  }
})
