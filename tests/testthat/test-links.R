# Tests for the mean-link inverse transformation simplex_linkinv().
# The four supported links and their integer codes are frozen by the contract:
#   logit = 1, probit = 2, cloglog = 3, neglog = 4.
# We check simplex_linkinv against the closed-form inverses and round-trip the
# analytic forward link g(mu) for interior mu.

eta <- seq(-4, 4, by = 0.25)

# Closed-form inverse links (mu = g^{-1}(eta)).
inv_logit   <- function(eta) 1 / (1 + exp(-eta))
inv_probit  <- function(eta) stats::pnorm(eta)
inv_cloglog <- function(eta) 1 - exp(-exp(eta))
inv_neglog  <- function(eta) exp(-exp(-eta))  # neglog inverse, per Zhang et al. (2016)

# Forward links (eta = g(mu)).
g_logit   <- function(mu) stats::qlogis(mu)
g_probit  <- function(mu) stats::qnorm(mu)
g_cloglog <- function(mu) log(-log1p(-mu))
g_neglog  <- function(mu) -log(-log(mu))

test_that("simplex_linkinv reproduces the logit inverse", {
  expect_equal(simplex_linkinv(eta, "logit"), inv_logit(eta), tolerance = 1e-9)
})

test_that("simplex_linkinv reproduces the probit inverse", {
  expect_equal(simplex_linkinv(eta, "probit"), inv_probit(eta), tolerance = 1e-9)
})

test_that("simplex_linkinv reproduces the cloglog inverse", {
  expect_equal(simplex_linkinv(eta, "cloglog"), inv_cloglog(eta), tolerance = 1e-9)
})

test_that("simplex_linkinv reproduces the neglog inverse", {
  expect_equal(simplex_linkinv(eta, "neglog"), inv_neglog(eta), tolerance = 1e-9)
})

test_that("simplex_linkinv output stays inside (0, 1)", {
  for (lk in c("logit", "probit", "cloglog", "neglog")) {
    mu <- simplex_linkinv(eta, lk)
    expect_true(all(mu > 0 & mu < 1), info = lk)
  }
})

test_that("simplex_linkinv round-trips against the forward link", {
  mu_interior <- seq(0.05, 0.95, by = 0.05)
  fwd <- list(
    logit   = g_logit,
    probit  = g_probit,
    cloglog = g_cloglog,
    neglog  = g_neglog
  )
  for (lk in names(fwd)) {
    eta_from_mu <- fwd[[lk]](mu_interior)
    mu_back <- simplex_linkinv(eta_from_mu, lk)
    expect_equal(mu_back, mu_interior, tolerance = 1e-7, info = lk)
  }
})

test_that("simplex_linkinv defaults to logit and matches match.arg", {
  expect_equal(simplex_linkinv(eta), simplex_linkinv(eta, "logit"))
})
