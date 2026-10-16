# Tail accuracy of psimplex() and qsimplex().
#
# The upper tail used to lose digits to cancellation of two numbers of size
# 2 / (phi mu (1 - mu)) when mu is tiny (abs. error 1e-5 at the mode for
# mu = 1e-8), and qsimplex() formed 1 - p, so p = 1e-20 returned 1. Every
# reference below is independent of the compiled code: identities (F + S = 1),
# a Gauss-Legendre quadrature of the density coded here, numerical derivatives,
# and round trips.

sd_simplex <- function(mu, phi) sqrt(phi * (mu * (1 - mu))^3)

# 64-point Gauss-Legendre rule (Golub-Welsch); the integrand is analytic on each
# segment used below, so the quadrature error is ~1e-15.
gl_rule <- local({
  n <- 64L
  i <- seq_len(n - 1L)
  b <- i / sqrt(4 * i^2 - 1)
  J <- matrix(0, n, n)
  J[cbind(i, i + 1L)] <- b
  J[cbind(i + 1L, i)] <- b
  e <- eigen(J, symmetric = TRUE)
  list(x = e$values, w = 2 * e$vectors[1L, ]^2)
})
gl_int <- function(f, a, b) {
  h <- 0.5 * (b - a)
  sum(gl_rule$w * f(0.5 * (a + b) + h * gl_rule$x)) * h
}

# Lower and upper tail probabilities at z0 standard deviations from mu, by
# integrating the density in z = (y - mu) / sd. The offset y - mu = z * sd is
# formed exactly, so the reference does not suffer the rounding of y = mu + z sd
# that limits integrate(dsimplex) to ~1e-11 when mu = 1e-10.
tails_by_quadrature <- function(z0, mu, phi) {
  sd <- sd_simplex(mu, phi)
  u <- mu * (1 - mu)
  dens <- function(z) {
    y <- mu + z * sd
    yb <- 1 - y
    exp(-0.5 * log(2 * pi * phi) - 1.5 * log(y * yb) -
          0.5 * (z * sd)^2 / (y * yb * u^2) / phi) * sd
  }
  zl <- max(-12, -0.999999 * mu / sd)
  zr <- min(12, 0.999999 * (1 - mu) / sd)
  seg <- function(a, b) {
    k <- c(a, c(-12, -8, -5, -3, -1, 0, 1, 3, 5, 8, 12)[
      c(-12, -8, -5, -3, -1, 0, 1, 3, 5, 8, 12) > a &
        c(-12, -8, -5, -3, -1, 0, 1, 3, 5, 8, 12) < b], b)
    sum(vapply(seq_len(length(k) - 1L),
               function(j) gl_int(dens, k[j], k[j + 1L]), numeric(1)))
  }
  c(lower = seg(zl, z0), upper = seg(z0, zr))
}

test_that("F + S = 1 to rounding, including tiny mu and mu near 1", {
  mus <- c(1e-12, 1e-8, 1e-4, 0.3, 0.5, 1 - 1e-4, 1 - 1e-8)
  z <- seq(-6, 6, by = 0.5)
  for (mu in mus) {
    for (phi in c(1e-3, 1)) {
      q <- c(mu + z * sd_simplex(mu, phi), 10^seq(-12, -0.01, length.out = 25))
      q <- q[q > 0 & q < 1]
      Fq <- psimplex(q, mu, phi)
      Sq <- psimplex(q, mu, phi, lower.tail = FALSE)
      keep <- abs(Fq - 0.5) < 0.499
      expect_gt(sum(keep), 5L)
      expect_lt(max(abs(Fq[keep] + Sq[keep] - 1)), 1e-13)
    }
  }
})

test_that("both tails match an independent quadrature, down to mu = 1e-10", {
  cases <- expand.grid(mu = c(1e-10, 1e-8, 1e-6, 0.01), phi = c(0.1, 1),
                       z = c(-3, 0, 3))
  for (k in seq_len(nrow(cases))) {
    mu <- cases$mu[k]; phi <- cases$phi[k]
    q <- mu + cases$z[k] * sd_simplex(mu, phi)
    ref <- tails_by_quadrature((q - mu) / sd_simplex(mu, phi), mu, phi)
    info <- sprintf("mu = %g, phi = %g, z = %g", mu, phi, cases$z[k])
    # Absolute error: both tails are probabilities of O(1) at these points.
    expect_lt(abs(psimplex(q, mu, phi, lower.tail = FALSE) - ref[["upper"]]), 1e-13,
              label = paste("upper:", info))
    expect_lt(abs(psimplex(q, mu, phi) - ref[["lower"]]), 1e-13,
              label = paste("lower:", info))
  }
})

test_that("d/dq of both tails is -/+ the density (numDeriv), tiny mu included", {
  skip_if_not_installed("numDeriv")
  for (mu in c(1e-10, 1e-8, 1e-4, 0.05, 0.3, 0.8)) {
    for (phi in c(0.05, 1)) {
      for (z in c(-2, 0, 2)) {
        sd <- sd_simplex(mu, phi)
        q0 <- mu + z * sd
        # Step of 0.05 sd; zero.tol = 0 stops numDeriv adding its absolute
        # floor, which for |q0| < 2e-5 would dwarf the sd of a tiny-mu law.
        ma <- list(d = 0.05 * min(sd, q0) / q0, zero.tol = 0)
        f0 <- dsimplex(q0, mu, phi)
        gS <- numDeriv::grad(function(q) psimplex(q, mu, phi, lower.tail = FALSE),
                             q0, method.args = ma)
        gF <- numDeriv::grad(function(q) psimplex(q, mu, phi), q0, method.args = ma)
        info <- sprintf("mu = %g, phi = %g, z = %g", mu, phi, z)
        expect_lt(abs(gS + f0) / f0, 1e-6, label = paste("upper:", info))
        expect_lt(abs(gF - f0) / f0, 1e-6, label = paste("lower:", info))
      }
    }
  }
})

test_that("log tails stay finite, non-positive and ordered at extreme parameters", {
  q <- c(1e-300, 1e-12, 1e-6, 0.01, 0.3, 0.5, 0.7, 0.99, 1 - 1e-6, 1 - 1e-12)
  for (par in list(c(1e-12, 1e-6), c(1e-12, 1e6), c(1 - 1e-12, 1e-6),
                   c(1 - 1e-12, 1e6), c(0.5, 1e-6), c(0.5, 1e6))) {
    lF <- psimplex(q, par[1L], par[2L], log.p = TRUE)
    lS <- psimplex(q, par[1L], par[2L], lower.tail = FALSE, log.p = TRUE)
    info <- paste(par, collapse = ", ")
    expect_false(anyNA(c(lF, lS)), info = info)
    expect_true(all(lF <= 0 & lS <= 0), info = info)
    # F is non-decreasing and S non-increasing in q (to rounding of the log).
    expect_true(all(diff(lF) >= -1e-12), info = info)
    expect_true(all(diff(lS) <= 1e-12), info = info)
  }
})

test_that("qsimplex upper tail keeps tiny probabilities (no 1 - p)", {
  # Reference values: S(q) = p solved by the reflection q = 1 - q_lower(p; 1 - mu)
  # at mu = 1/2, and checked by a round trip below.
  q20 <- qsimplex(1e-20, 0.5, 1, lower.tail = FALSE)
  expect_equal(q20, 0.9590248069, tolerance = 1e-9)
  expect_equal(psimplex(q20, 0.5, 1, lower.tail = FALSE, log.p = TRUE), log(1e-20),
               tolerance = 1e-12)
  # p = 1e-16 used to give S(q) = 1.6e-16.
  q16 <- qsimplex(1e-16, 0.5, 1, lower.tail = FALSE)
  expect_equal(psimplex(q16, 0.5, 1, lower.tail = FALSE), 1e-16, tolerance = 1e-10)
  # log.p = TRUE: -100 used to return 1; -1000 underflows exp() on the lower tail.
  for (lp in c(-100, -1000)) {
    for (lower in c(TRUE, FALSE)) {
      q <- qsimplex(lp, 0.5, 1, lower.tail = lower, log.p = TRUE)
      expect_true(q > 0 && q < 1, info = paste(lp, lower))
      expect_equal(psimplex(q, 0.5, 1, lower.tail = lower, log.p = TRUE), lp,
                   tolerance = 1e-10, info = paste(lp, lower))
    }
  }
  # The tails mirror each other at mu = 1/2: q_upper(p) = 1 - q_lower(p).
  expect_equal(qsimplex(1e-8, 0.5, 2, lower.tail = FALSE),
               1 - qsimplex(1e-8, 0.5, 2), tolerance = 1e-12)
})

test_that("qsimplex boundary conventions hold for both tails and log.p", {
  expect_identical(qsimplex(c(0, 1), 0.4, 1, lower.tail = FALSE), c(1, 0))
  expect_identical(qsimplex(c(-Inf, 0), 0.4, 1, log.p = TRUE), c(0, 1))
  expect_identical(qsimplex(c(-Inf, 0), 0.4, 1, lower.tail = FALSE, log.p = TRUE),
                   c(1, 0))
  expect_warning(expect_identical(
    qsimplex(0.1, 0.4, 1, lower.tail = FALSE, log.p = TRUE), NaN), "NaNs produced")
})

test_that("qsimplex inverts psimplex for p from 1e-300 to 1, both tails", {
  # Moving q by one ulp moves log p by roughly (y f / T) * 1e-16, i.e. up to
  # ~1e-13 for these laws at p = 1e-300, so 1e-10 is the honest bound there and
  # 1e-12 holds away from the extreme tail.
  for (par in list(c(0.3, 1), c(0.5, 0.01), c(0.05, 2))) {
    mu <- par[1L]; phi <- par[2L]
    for (lower in c(TRUE, FALSE)) {
      lp_deep <- log(10^seq(-300, -8, length.out = 40))
      q <- qsimplex(exp(lp_deep), mu, phi, lower.tail = lower)
      expect_true(all(q > 0 & q < 1))
      expect_lt(max(abs(psimplex(q, mu, phi, lower.tail = lower, log.p = TRUE) - lp_deep)),
                1e-10)
      p_mid <- 10^seq(-8, -0.01, length.out = 40)
      q <- qsimplex(p_mid, mu, phi, lower.tail = lower)
      expect_lt(max(abs(psimplex(q, mu, phi, lower.tail = lower, log.p = TRUE) -
                          log(p_mid))), 1e-12)
      # log.p input over a range no linear p can represent.
      lp <- seq(-1000, -1e-3, length.out = 40)
      q <- qsimplex(lp, mu, phi, lower.tail = lower, log.p = TRUE)
      expect_lt(max(abs(psimplex(q, mu, phi, lower.tail = lower, log.p = TRUE) - lp)),
                1e-10)
    }
  }
})

test_that("qsimplex resolves a tiny-mu law to far below one standard deviation", {
  # At mu = 1e-10 the sd is 1e-15: an absolute step tolerance of 1e-15 (the old
  # rule) stops at about one sd, so the round trip was off by O(0.1) in p.
  mu <- 1e-10; phi <- 1
  p <- c(0.01, 0.1, 0.5, 0.9, 0.99)
  for (lower in c(TRUE, FALSE)) {
    q <- qsimplex(p, mu, phi, lower.tail = lower)
    expect_equal(psimplex(q, mu, phi, lower.tail = lower), p, tolerance = 1e-9)
  }
  # A probability close to 1 is solved on the opposite tail, with 1 - p exact.
  q <- qsimplex(1 - 1e-10, 0.3, 1)
  expect_equal(psimplex(q, 0.3, 1, lower.tail = FALSE), 1e-10, tolerance = 1e-8)
})
