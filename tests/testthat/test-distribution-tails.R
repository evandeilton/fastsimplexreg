# Tail accuracy of psimplex() and qsimplex().
#
# The upper tail used to lose digits to cancellation of two numbers of size
# 2 / (phi mu (1 - mu)) when mu is tiny (abs. error 1e-5 at the mode for
# mu = 1e-8), and qsimplex() formed 1 - p, so p = 1e-20 returned 1. Every
# reference below is independent of the compiled code: identities (F + S = 1),
# a Gauss-Legendre quadrature of the density coded here, numerical derivatives,
# and round trips.

sd_simplex <- function(mu, phi) sqrt(phi * (mu * (1 - mu))^3)

# Relative error. expect_equal() would NOT do here: like all.equal(), it switches
# to an absolute difference once |target| is below the tolerance, which is the
# regime (log T ~ -1e-20) these tests are about.
rel_err <- function(got, ref) abs(got - ref) / abs(ref)

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
  expect_lt(rel_err(psimplex(q16, 0.5, 1, lower.tail = FALSE), 1e-16), 1e-10)
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
  # A probability close to 1 is solved on the opposite tail, with 1 - p exact
  # (p itself is a double, so the target is 1 - p = 1.00000008e-10, not 1e-10).
  p1 <- 1 - 1e-10
  q <- qsimplex(p1, 0.3, 1)
  expect_lt(rel_err(psimplex(q, 0.3, 1, lower.tail = FALSE), 1 - p1), 1e-9)
})

# A tail above 1/2 has log T = -(the other tail) to first order, which a log of a
# sum cannot resolve once the other tail is below ~1e-16 (psimplex(0.9, 0.1, 10,
# log.p = TRUE) was -3.6e-21 for a true -7.97e-22, and the relative error of log F
# was 8.5e-5 at 1e-12). The complement of the small tail is used instead.
test_that("log F and log S near 0 are exact: 140-digit values", {
  # Columns: q, mu, phi, lower (1/0), log T. The reference is the closed form
  # F = Phi(a) + (1-2mu) phi(a) m(beta) in 140-digit decimal arithmetic from the
  # exact binary inputs, kept out of this package (Python decimal); log T is
  # about -1e-9 or -1e-20.
  # q, mu and phi are hexadecimal doubles on purpose: R's decimal parser is not
  # correctly rounded for 17 digits where long double is double (macOS arm64), and
  # for mu = 0.9999 one ulp of q moves log T by 7e-10 relative.
  ref <- rbind(
    c(0x1.bd537ae851f25p-14, 0x1.a36e2eb1c432dp-14, 0x1.0000000000000p+0, 1, -9.999999999999589e-10),
    c(0x1.cc1716b30ab02p-14, 0x1.a36e2eb1c432dp-14, 0x1.0000000000000p+0, 1, -1.000000000000050e-20),
    c(0x1.7d17b4931b363p-1, 0x1.47ae147ae147bp-7, 0x1.f400000000000p+9, 1, -9.999999999999893e-10),
    c(0x1.c623f2355c662p-1, 0x1.47ae147ae147bp-7, 0x1.f400000000000p+9, 1, -9.999999999999959e-21),
    c(0x1.8f4ec35e4e4fcp-1, 0x1.999999999999ap-4, 0x1.4000000000000p+3, 1, -1.000000000000001e-09),
    c(0x1.ca1befce55cd6p-1, 0x1.999999999999ap-4, 0x1.4000000000000p+3, 1, -9.999999999999887e-21),
    c(0x1.992ea30a81a7fp-1, 0x1.3333333333333p-2, 0x1.0000000000000p+0, 1, -1.000000000000003e-09),
    c(0x1.eb054cabbe43fp-1, 0x1.0000000000000p-1, 0x1.0000000000000p+0, 1, -9.999999999999882e-21),
    c(0x1.b08ae133dcdf0p-1, 0x1.6666666666666p-1, 0x1.999999999999ap-4, 1, -1.000000000000011e-09),
    c(0x1.ff4ab86664d17p-1, 0x1.ccccccccccccdp-1, 0x1.4000000000000p+3, 1, -9.999999999999056e-21),
    c(0x1.fff3a7fe15415p-1, 0x1.fff2e48e8a71ep-1, 0x1.0000000000000p+0, 1, -1.000000000335342e-09),
    c(0x1.7e534b47bc4f2p-14, 0x1.a36e2eb1c432dp-14, 0x1.0000000000000p+0, 0, -9.999999999999820e-21),
    c(0x1.6a8f33365d20dp-10, 0x1.999999999999ap-4, 0x1.4000000000000p+3, 0, -9.999999999999944e-21),
    c(0x1.5f49229e9e4d2p-5, 0x1.3333333333333p-2, 0x1.0000000000000p+0, 0, -1.000000000000001e-09),
    c(0x1.582717d6727c7p-4, 0x1.0000000000000p-1, 0x1.0000000000000p+0, 0, -1.000000000000002e-09),
    c(0x1.9930ebc74708dp-2, 0x1.6666666666666p-1, 0x1.999999999999ap-4, 0, -1.000000000000017e-20),
    c(0x1.af20818d51955p-4, 0x1.ccccccccccccdp-1, 0x1.4000000000000p+3, 0, -9.999999999999966e-21),
    c(0x1.fff19f474a67bp-1, 0x1.fff2e48e8a71ep-1, 0x1.0000000000000p+0, 0, -1.000000000218489e-20)
  )
  for (k in seq_len(nrow(ref))) {
    got <- psimplex(ref[k, 1L], ref[k, 2L], ref[k, 3L], lower.tail = ref[k, 4L] == 1,
                    log.p = TRUE)
    expect_lt(rel_err(got, ref[k, 5L]), 1e-11,
              label = paste("relative error at", paste(ref[k, 1:4], collapse = ", ")))
  }
  # The case that was wrong by a factor 4.6.
  expect_lt(rel_err(psimplex(0.9, 0.1, 10, log.p = TRUE), -7.96788938016e-22), 1e-10)
})

test_that("log F and log S down to |log T| = 1e-300 are minus the opposite tail", {
  # Pick q from the small opposite tail p, so that log T = log1p(-p) = -p to 1e-16
  # (the exact complement), and compare with the log.p value of the large tail.
  p <- 10^c(-8, -10, -12, -16, -20, -30, -100, -300)
  for (par in list(c(0.1, 1), c(0.5, 1), c(0.9, 10), c(1e-4, 1), c(0.3, 0.01))) {
    for (lower in c(TRUE, FALSE)) {
      q <- qsimplex(p, par[1L], par[2L], lower.tail = !lower)
      small <- psimplex(q, par[1L], par[2L], lower.tail = !lower)   # linear scale, accurate
      big <- psimplex(q, par[1L], par[2L], lower.tail = lower, log.p = TRUE)
      expect_true(all(big <= 0 & is.finite(big)), info = paste(par, lower))
      expect_lt(max(rel_err(big, log1p(-small))), 1e-11,
                label = paste("relative error at", paste(par, lower, collapse = ", ")))
    }
  }
})

test_that("F and S are monotone in q on dense grids, mu = 0.1 and 0.9", {
  # Near F = 1 the old log-of-a-sum made both tails move in 1-ulp steps the wrong way.
  for (mu in c(0.1, 0.9)) {
    for (phi in c(1e-3, 0.1, 1)) {
      qs <- list(plogis(seq(qlogis(1e-6), qlogis(1 - 1e-6), length.out = 20001L)))
      qs[[2L]] <- mu + seq(-8, 8, length.out = 20001L) * sd_simplex(mu, phi)
      med <- qsimplex(0.5, mu, phi)
      qs[[3L]] <- med * (1 + seq(-2000, 2000) * 1e-13)       # micro grid across F = 1/2
      for (q in qs) {
        q <- sort(unique(q[q > 0 & q < 1]))
        info <- paste(mu, phi)
        expect_true(all(diff(psimplex(q, mu, phi, log.p = TRUE)) >= 0), info = info)
        expect_true(all(diff(psimplex(q, mu, phi, lower.tail = FALSE, log.p = TRUE)) <= 0),
                    info = info)
        expect_true(all(diff(psimplex(q, mu, phi)) >= 0), info = info)
        expect_true(all(diff(psimplex(q, mu, phi, lower.tail = FALSE)) <= 0), info = info)
      }
    }
  }
})

test_that("qsimplex round-trips log p = -1e-30 and -1e-12 on both tails", {
  # psimplex(qsimplex(-1e-30, 0.1, 1, log.p = TRUE), log.p = TRUE) was -3.2e-30.
  for (par in list(c(0.1, 1), c(0.5, 1), c(0.9, 1), c(0.3, 10))) {
    for (lower in c(TRUE, FALSE)) {
      for (lp in c(-1e-12, -1e-30)) {
        q <- qsimplex(lp, par[1L], par[2L], lower.tail = lower, log.p = TRUE)
        expect_true(q > 0 && q < 1, info = paste(par, lower, lp))
        expect_lt(rel_err(psimplex(q, par[1L], par[2L], lower.tail = lower, log.p = TRUE), lp),
                  1e-9, label = paste("relative error at", paste(par, lower, lp, collapse = ", ")))
      }
    }
  }
})
