# Parallelism is a headline feature of this package, so it is tested rather than
# assumed: results must not depend on the thread count beyond floating-point
# reassociation, and long calls must be interruptible.

test_that("the distribution functions are bit-identical across thread counts", {
  set.seed(4)
  n <- 2e5L
  q <- runif(n, 0.01, 0.99)
  mu <- runif(n, 0.05, 0.95)
  phi <- exp(runif(n, -2, 2))

  # Pure element-wise maps: no reduction, so there is nothing to reassociate.
  for (k in c(2L, 4L, 8L)) {
    expect_identical(dsimplex(q, mu, phi, n_threads = k),
                     dsimplex(q, mu, phi, n_threads = 1L))
    expect_identical(psimplex(q, mu, phi, n_threads = k),
                     psimplex(q, mu, phi, n_threads = 1L))
    expect_identical(qsimplex(q, mu, phi, n_threads = k),
                     qsimplex(q, mu, phi, n_threads = 1L))
  }
})

test_that("the objective agrees across thread counts to floating-point noise", {
  set.seed(5)
  n <- 5e4L
  p <- 6L
  X <- cbind(1, matrix(rnorm(n * (p - 1L)), n, p - 1L))
  Z <- cbind(1, rnorm(n))
  b <- c(0.2, rep(0.15, p - 1L))
  y <- rsimplex(n, simplex_linkinv(as.numeric(X %*% b), "logit"), 1)
  th <- c(b, -0.3, 0.2)

  ref <- fastsimplexreg:::simplex_eval_cpp(th, y, X, Z, 1L, 1L)
  for (link in 1:4) {
    r1 <- fastsimplexreg:::simplex_eval_cpp(th, y, X, Z, link, 1L)
    for (k in c(2L, 4L, 8L)) {
      rk <- fastsimplexreg:::simplex_eval_cpp(th, y, X, Z, link, k)
      # The per-thread accumulators are summed in thread order, so the total is
      # reassociated -- relative differences must stay at the rounding level.
      expect_lt(abs(rk$value - r1$value) / max(abs(r1$value), 1), 1e-12)
      expect_lt(max(abs(rk$gradient - r1$gradient)) /
                  max(max(abs(r1$gradient)), 1), 1e-12)
    }
  }
  expect_true(is.finite(ref$value))
})

test_that("a fit gives the same answer at any thread count", {
  set.seed(6)
  n <- 3e4L
  d <- data.frame(x1 = rnorm(n), x2 = rbinom(n, 1L, 0.4), z1 = rnorm(n))
  d$y <- rsimplex(n, simplex_linkinv(0.2 + 0.7 * d$x1 - 0.4 * d$x2, "logit"),
                  exp(-0.5 + 0.4 * d$z1))
  f1 <- fastsimplexreg(y ~ x1 + x2 | z1, data = d, n_threads = 1L)
  for (k in c(2L, 4L)) {
    fk <- fastsimplexreg(y ~ x1 + x2 | z1, data = d, n_threads = k)
    expect_equal(fk$logLik, f1$logLik, tolerance = 1e-8)
    expect_equal(coef(fk), coef(f1), tolerance = 1e-6)
    expect_equal(fk$standard_errors, f1$standard_errors, tolerance = 1e-5)
  }
})

test_that("a long call can be interrupted", {
  skip_on_cran()
  skip_on_os("windows")
  skip_if_not(nzchar(Sys.which("bash")), "bash is needed to deliver SIGINT")

  set.seed(1)
  n <- 2e6L
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  d$y <- rsimplex(n, simplex_linkinv(0.2 + 0.4 * d$x1 - 0.3 * d$x2, "logit"), 1)

  pid <- Sys.getpid()
  system2("bash", c("-c", shQuote(sprintf("sleep 2; kill -INT %d", pid))),
          wait = FALSE)
  t0 <- Sys.time()
  res <- tryCatch({
    fastsimplexreg(y ~ x1 + x2, data = d, n_threads = 1L, inference = FALSE)
    "completed"
  }, interrupt = function(e) "interrupted")
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  # The check sits at the top of each BFGS iteration, so the latency is one
  # objective evaluation. Without it the whole fit ran to completion (~10 s).
  expect_identical(res, "interrupted")
  expect_lt(elapsed, 6)
})
