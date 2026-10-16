# Extractor methods of the mixed fit: coef(model = "all"), deviance, model.matrix,
# terms and simulate.
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

test_that("coef(model = 'all') has unique, (phi)_-prefixed names", {
  cf <- coef(fr, "all")
  expect_identical(names(cf),
                   c("(Intercept)", "x1", "(phi)_(Intercept)", "(phi)_z1"))
  expect_false(anyDuplicated(names(cf)) > 0L)
  # Fixed effects only: the variance components are not part of coef().
  expect_length(cf, length(coef(fr, "mean")) + length(coef(fr, "dispersion")))
  expect_equal(unname(cf), unname(c(coef(fr, "mean"), coef(fr, "dispersion"))))
  # The per-submodel views keep the bare names.
  expect_identical(names(coef(fr, "dispersion")), c("(Intercept)", "z1"))
  # An aliased fit keeps its NA slots, under the same names.
  fa <- suppressWarnings(fit_mixed(y ~ x1 + x2 | z1 + z2, dat_a))
  expect_identical(names(coef(fa, "all")), names(fa$aliased))
  expect_true(is.na(coef(fa, "all")[["x2"]]))
  expect_true(is.na(coef(fa, "all")[["(phi)_z2"]]))
})

test_that("deviance() uses the unit deviance at the conditional fitted means", {
  y <- dat_a$y
  mu <- fitted(fr)
  phi <- fitted(fr, "dispersion")
  d <- (y - mu)^2 / (y * (1 - y) * mu^2 * (1 - mu)^2)
  expect_false(is.null(deviance(fr)))
  expect_length(deviance(fr), 1L)
  expect_equal(deviance(fr), sum(d))
  expect_equal(deviance(fr, type = "unscaled"), sum(d))
  expect_equal(deviance(fr, type = "scaled"), sum(d / phi))
  # Equals the sum of squared deviance residuals.
  expect_equal(deviance(fr, "scaled"), sum(residuals(fr, "deviance")^2))
  expect_error(deviance(fr, type = "other"))
  # Falls back to the stored residuals when the response was not kept.
  fn <- fit_mixed(y ~ x1 | z1, dat_a, y = FALSE)
  expect_equal(deviance(fn), deviance(fr), tolerance = 1e-6)
})

# Random intercept and slope, a factor and an offset in the mean: exercises
# contrasts, xlevels and the three design parts.
dat_q <- local({
  d <- sim_mixed_alias(J = 20L, nj = 6L, seed = 7L)
  set.seed(8L)
  d$f <- factor(sample(c("a", "b", "c"), nrow(d), replace = TRUE))
  d$o <- runif(nrow(d), 0, 0.2)
  d
})

test_that("model.matrix() returns the three design parts", {
  fx <- fastsimplexregmixed(y ~ x1 + f + offset(o) | z1, random = ~ 1 + x1 | g,
                            data = dat_q, nAGQ = 7L, n_threads = 1L, x = TRUE)
  fm <- fastsimplexregmixed(y ~ x1 + f + offset(o) | z1, random = ~ 1 + x1 | g,
                            data = dat_q, nAGQ = 7L, n_threads = 1L)
  n <- nrow(dat_q)

  mm <- model.matrix(fx)                       # default is the mean part
  expect_identical(dim(mm), c(n, 4L))
  expect_identical(colnames(mm), c("(Intercept)", "x1", "fb", "fc"))
  expect_identical(colnames(model.matrix(fx, "dispersion")), c("(Intercept)", "z1"))
  expect_identical(colnames(model.matrix(fx, "random")), c("(Intercept)", "x1"))

  # Stored (x = TRUE) and rebuilt-from-model-frame designs coincide.
  expect_null(fm$x)
  for (part in c("mean", "dispersion", "random")) {
    expect_equal(model.matrix(fm, part), model.matrix(fx, part), info = part)
  }
  expect_equal(unname(model.matrix(fm, "random")[, 2L]), dat_q$x1)

  # With neither the design nor the model frame, say so.
  fn <- fastsimplexregmixed(y ~ x1 | z1, random = ~ 1 | g, data = dat_q,
                            nAGQ = 7L, n_threads = 1L, model = FALSE)
  expect_error(model.matrix(fn), "Refit with x = TRUE or model = TRUE")
  expect_error(model.matrix(fx, "other"))

  # Without a dispersion part the dispersion design is the intercept.
  f1 <- fit_mixed(y ~ x1, dat_q)
  expect_identical(dim(model.matrix(f1, "dispersion")), c(n, 1L))
  expect_identical(colnames(model.matrix(f1, "dispersion")), "(Intercept)")
})

test_that("terms() returns terms objects for the three parts", {
  expect_s3_class(terms(fr), "terms")
  expect_s3_class(terms(fr, "dispersion"), "terms")
  expect_s3_class(terms(fr, "random"), "terms")
  expect_identical(attr(terms(fr), "term.labels"), "x1")        # default: mean
  expect_identical(attr(terms(fr, "dispersion"), "term.labels"), "z1")
  expect_identical(attr(terms(fr, "random"), "term.labels"), character(0))
  expect_error(terms(fr, "other"))
})

test_that("simulate() draws from the conditional fit and handles the seed", {
  n <- nrow(dat_a)
  s1 <- simulate(fr, nsim = 3L, seed = 42L)
  s2 <- simulate(fr, nsim = 3L, seed = 42L)
  expect_s3_class(s1, "data.frame")
  expect_identical(dim(s1), c(n, 3L))
  expect_identical(names(s1), paste0("sim_", 1:3))
  expect_true(all(as.matrix(s1) > 0 & as.matrix(s1) < 1))
  expect_identical(s1, s2)                              # reproducible
  expect_equal(c(attr(s1, "seed")), 42)
  expect_false(isTRUE(all.equal(s1$sim_1, s1$sim_2)))   # columns differ
  expect_false(isTRUE(all.equal(s1, simulate(fr, nsim = 3L, seed = 43L),
                                check.attributes = FALSE)))
  expect_identical(rownames(s1), rownames(fr$model))

  # A supplied seed leaves the caller's RNG stream untouched.
  set.seed(1L)
  before <- .Random.seed
  invisible(simulate(fr, nsim = 2L, seed = 5L))
  expect_identical(.Random.seed, before)
  # Without a seed the stream is used (and advanced), as in stats::simulate().
  set.seed(9L)
  start_state <- .Random.seed
  a <- simulate(fr, nsim = 1L)
  expect_false(identical(.Random.seed, start_state))
  set.seed(9L)
  expect_identical(simulate(fr, nsim = 1L), a)

  # Conditional on the estimated random effects: draws center on the
  # conditional fitted means, not on the population-level ones.
  big <- simulate(fr, nsim = 300L, seed = 1L)
  m <- rowMeans(big)
  expect_lt(mean(abs(m - fitted(fr))), 0.02)
  expect_gt(mean(abs(m - predict(fr, re.form = NA))), mean(abs(m - fitted(fr))))

  expect_error(simulate(fr, nsim = 0L), "positive integer")
  expect_error(simulate(fr, nsim = c(1L, 2L)), "positive integer")
})
