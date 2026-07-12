# Tests for the simplex density (dsimplex) and random generator (rsimplex).
# All tests use n_threads = 1L and modest sample sizes for speed and
# reproducibility. Correctness is verified purely through the documented public
# R API described in the consensus contract (Section 3).

# Closed-form reference for the simplex log-density, used to cross-check the
# compiled backend independently of any C++ symbol.
ref_log_dsimplex <- function(y, mu, phi) {
  dev <- (y - mu)^2 / (y * (1 - y) * (mu * (1 - mu))^2)
  -0.5 * (log(2 * pi) + log(phi)) - 1.5 * (log(y) + log(1 - y)) - 0.5 * dev / phi
}

test_that("dsimplex integrates to ~1 over (0, 1) for several (mu, phi)", {
  grid <- expand.grid(
    mu  = c(0.2, 0.5, 0.75),
    phi = c(0.5, 1, 2)
  )
  for (k in seq_len(nrow(grid))) {
    mu  <- grid$mu[k]
    phi <- grid$phi[k]
    integrand <- function(x) dsimplex(x, mu = mu, phi = phi, n_threads = 1L)
    val <- stats::integrate(
      integrand,
      lower = 0, upper = 1,
      rel.tol = 1e-8, subdivisions = 500L
    )$value
    expect_equal(
      val, 1,
      tolerance = 1e-4,
      info = sprintf("integral for mu = %g, phi = %g was %g", mu, phi, val)
    )
  }
})

test_that("log vs non-log densities are consistent", {
  x <- c(0.05, 0.2, 0.4, 0.6, 0.85, 0.95)
  d_lin <- dsimplex(x, mu = 0.45, phi = 1.3, log = FALSE, n_threads = 1L)
  d_log <- dsimplex(x, mu = 0.45, phi = 1.3, log = TRUE,  n_threads = 1L)
  expect_equal(exp(d_log), d_lin, tolerance = 1e-12)
  # Independent closed-form cross-check.
  expect_equal(d_log, ref_log_dsimplex(x, 0.45, 1.3), tolerance = 1e-10)
})

test_that("length-1 mu and phi are recycled against x", {
  x <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  d_scalar <- dsimplex(x, mu = 0.5, phi = 0.8, n_threads = 1L)
  d_recyc  <- dsimplex(x, mu = rep(0.5, length(x)), phi = rep(0.8, length(x)),
                       n_threads = 1L)
  expect_equal(d_scalar, d_recyc, tolerance = 1e-12)
  expect_length(d_scalar, length(x))
})

test_that("vectorised mu and phi are applied element-wise", {
  x   <- c(0.2, 0.4, 0.6)
  mu  <- c(0.3, 0.5, 0.7)
  phi <- c(0.5, 1.0, 2.0)
  got <- dsimplex(x, mu = mu, phi = phi, log = TRUE, n_threads = 1L)
  want <- ref_log_dsimplex(x, mu, phi)
  expect_equal(got, want, tolerance = 1e-10)
})

test_that("out-of-support x gives 0 (or -Inf for log)", {
  x <- c(-0.1, 0, 0.5, 1, 1.2)
  d_lin <- dsimplex(x, mu = 0.5, phi = 1, log = FALSE, n_threads = 1L)
  d_log <- dsimplex(x, mu = 0.5, phi = 1, log = TRUE,  n_threads = 1L)
  # Only the interior point 0.5 is on the support.
  expect_equal(d_lin[c(1, 2, 4, 5)], rep(0, 4))
  expect_true(all(is.infinite(d_log[c(1, 2, 4, 5)]) & d_log[c(1, 2, 4, 5)] < 0))
  expect_true(d_lin[3] > 0)
  expect_true(is.finite(d_log[3]))
})

test_that("rsimplex draws lie strictly inside (0, 1)", {
  set.seed(101)
  x <- rsimplex(5000L, mu = 0.4, phi = 0.7)
  expect_length(x, 5000L)
  expect_true(all(x > 0 & x < 1))
  expect_true(all(is.finite(x)))
})

test_that("rsimplex sample mean approximates mu", {
  # Large-ish but still fast; the simplex mean equals mu exactly, so with
  # n = 2e4 a generous Monte-Carlo tolerance is comfortably satisfied.
  set.seed(202)
  for (mu in c(0.3, 0.55, 0.8)) {
    x <- rsimplex(2e4L, mu = mu, phi = 1)
    expect_equal(
      mean(x), mu,
      tolerance = 0.02,
      info = sprintf("sample mean %g vs mu %g", mean(x), mu)
    )
  }
})

test_that("rsimplex is reproducible across two seeded calls", {
  set.seed(999)
  a <- rsimplex(2000L, mu = 0.5, phi = 1.2)
  set.seed(999)
  b <- rsimplex(2000L, mu = 0.5, phi = 1.2)
  expect_identical(a, b)
})

test_that("rsimplex recycles length-1 mu and phi", {
  set.seed(7)
  x <- rsimplex(1000L, mu = 0.5, phi = 1)
  expect_length(x, 1000L)
  expect_true(all(x > 0 & x < 1))
})

test_that("dsimplex matches a direct implementation of the density formula", {
  # Independent, dependency-free check of the C++ density against the closed-form
  # simplex log-density: an alternative code path evaluated in pure R.
  ld_ref <- function(y, mu, phi) {
    dev <- (y - mu)^2 / (y * (1 - y) * mu^2 * (1 - mu)^2)
    -0.5 * (log(2 * pi) + log(phi)) - 1.5 * (log(y) + log(1 - y)) - 0.5 * dev / phi
  }
  x_grid   <- seq(0.05, 0.95, by = 0.05)
  mu_grid  <- c(0.3, 0.5, 0.7)
  phi_grid <- c(0.5, 1, 2)
  for (mu in mu_grid) {
    for (phi in phi_grid) {
      ours <- dsimplex(x_grid, mu = mu, phi = phi, log = TRUE, n_threads = 1L)
      expect_equal(
        ours, ld_ref(x_grid, mu, phi),
        tolerance = 1e-10,
        info = sprintf("mu = %g, phi = %g", mu, phi)
      )
    }
  }
})
