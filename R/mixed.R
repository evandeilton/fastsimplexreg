# mixed.R
# R interface for the two-level simplex mixed model with variable dispersion.

# Internal: parse a random-effects specification `~ terms | group` into a
# one-sided random-effects formula and the grouping variable name. Enforces the
# v1 contract of a single grouping factor (two-level nesting).
.parse_random <- function(random) {
  if (!inherits(random, "formula") || length(random) != 2L) {
    stop("'random' must be a one-sided formula such as ~ 1 + x1 | group.", call. = FALSE)
  }
  bar <- random[[2L]]
  if (!(is.call(bar) && identical(bar[[1L]], as.name("|")))) {
    stop("'random' must contain a single grouping bar, e.g. ~ 1 + x1 | group.", call. = FALSE)
  }
  re_expr <- bar[[2L]]
  grp_expr <- bar[[3L]]
  if (is.call(re_expr) && identical(re_expr[[1L]], as.name("|"))) {
    stop("Only one grouping factor is supported in this version (2-level nesting).", call. = FALSE)
  }
  if (!is.name(grp_expr)) {
    stop("The grouping factor must be a single variable, e.g. ~ 1 + x1 | group.", call. = FALSE)
  }
  list(
    re_formula = stats::reformulate(deparse(re_expr, width.cutoff = 500L)),
    group = as.character(grp_expr)
  )
}


# Internal: assemble the response, the three design matrices (mean X,
# dispersion W, random Z) and the grouping factor from ONE shared model frame,
# so that NA handling, factor levels and contrasts are aligned across all parts.
.build_simplex_mixed_matrices <- function(formula, random, data,
                                          subset = NULL, na.action = stats::na.omit) {
  if (!requireNamespace("Formula", quietly = TRUE)) {
    stop("Package 'Formula' is required. Install it with install.packages('Formula').", call. = FALSE)
  }
  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula such as y ~ x1 + x2 | z1.", call. = FALSE)
  }

  Fu <- Formula::Formula(formula)
  dims <- length(Fu)
  if (dims[1L] != 1L) {
    stop("The model must contain exactly one response component.", call. = FALSE)
  }
  if (dims[2L] < 1L || dims[2L] > 2L) {
    stop("Use one or two RHS components: y ~ mean_terms | dispersion_terms.", call. = FALSE)
  }
  has_dispersion <- dims[2L] == 2L

  re <- .parse_random(random)

  # `subset` arrives here ALREADY RESOLVED to a plain index/logical vector by
  # .simplex_eval_subset(); see R/simplexreg.R. Applied here rather than
  # forwarded to model.frame(), whose own NSE would resolve the bare symbol in
  # the formula's environment instead of using the value supplied here.
  if (!is.null(subset)) {
    data <- data[subset, , drop = FALSE]
  }

  # One combined Formula: y ~ mean | dispersion | random | group. Building a
  # single model frame guarantees identical row dropping across all four parts.
  mean_f <- stats::formula(Fu, lhs = 1L, rhs = 1L)
  disp_f <- if (has_dispersion) stats::formula(Fu, lhs = 0L, rhs = 2L) else ~1
  grp_f <- stats::reformulate(re$group)
  Fc <- Formula::as.Formula(mean_f, disp_f, re$re_formula, grp_f)

  mf <- stats::model.frame(Fc, data = data, na.action = na.action, drop.unused.levels = TRUE)

  response <- Formula::model.part(Fc, data = mf, lhs = 1L, drop = TRUE)
  X <- stats::model.matrix(Fc, data = mf, rhs = 1L)
  W <- stats::model.matrix(Fc, data = mf, rhs = 2L)
  Z <- stats::model.matrix(Fc, data = mf, rhs = 3L)
  grp <- Formula::model.part(Fc, data = mf, rhs = 4L, drop = TRUE)

  if (!is.numeric(response)) {
    stop("The response must be numeric.", call. = FALSE)
  }
  response <- as.numeric(response)
  if (any(!is.finite(response)) || any(response <= 0 | response >= 1)) {
    stop("All response values must be finite and strictly inside (0, 1).", call. = FALSE)
  }
  if (ncol(Z) < 1L) {
    stop("The random-effects design must contain at least one term.", call. = FALSE)
  }

  storage.mode(X) <- "double"
  storage.mode(W) <- "double"
  storage.mode(Z) <- "double"

  n_obs <- nrow(X)
  tm <- stats::terms(Fc, rhs = 1L)
  td <- stats::terms(Fc, rhs = 2L)
  offset_mu <- .simplex_check_offset(.simplex_offset(tm, mf), n_obs, "mean")
  offset_phi <- .simplex_check_offset(.simplex_offset(td, mf), n_obs, "dispersion")

  group <- droplevels(as.factor(grp))

  list(
    formula = Fu,
    random = random,
    combined = Fc,
    model = mf,
    y = response,
    X = X, W = W, Z = Z,
    offset_mu = offset_mu,
    offset_phi = offset_phi,
    group = group,
    group_name = re$group,
    terms_mean = stats::terms(Fc, rhs = 1L),
    terms_dispersion = stats::terms(Fc, rhs = 2L),
    terms_random = stats::terms(Fc, rhs = 3L),
    xlevels_mean = stats::.getXlevels(stats::terms(Fc, rhs = 1L), mf),
    xlevels_dispersion = stats::.getXlevels(stats::terms(Fc, rhs = 2L), mf),
    xlevels_random = stats::.getXlevels(stats::terms(Fc, rhs = 3L), mf),
    contrasts_mean = attr(X, "contrasts"),
    contrasts_dispersion = attr(W, "contrasts"),
    contrasts_random = attr(Z, "contrasts"),
    has_dispersion_formula = has_dispersion
  )
}


# Internal: starting values c(beta, gamma, omega). Fixed effects reuse the
# marginal (random-effect-ignoring) starting values; omega starts from a small
# non-degenerate covariance Sigma0 = 0.1 * I.
.simplex_mixed_start <- function(y, X, W, Z, link) {
  fixed <- .simplex_start(y, X, W, link = link)
  q <- ncol(Z)
  Sigma0 <- diag(0.1, q, q)
  D0 <- t(chol(Sigma0))  # lower-triangular Cholesky
  omega0 <- simplex_mixed_omega_from_D_cpp(D0)
  c(fixed, omega0)
}


# Internal: labels for the packed omega (variance-component) parameters, in the
# same column-major lower-triangular order the C++ backend uses.
.omega_labels <- function(re_names) {
  q <- length(re_names)
  labs <- character(0)
  for (j in seq_len(q)) {
    # The diagonal entry is log D[j, j], where Sigma = D D'. That is the
    # marginal standard deviation ONLY for j = 1; for j >= 2 it is the
    # CONDITIONAL standard deviation of random effect j given effects 1..j-1.
    # Labelling it "logsd." understated a random slope's SD by a factor of
    # nearly two in a measured q = 2 fit (0.389 against a true 0.758), read
    # straight off fit$omega. "logchol." names what the number actually is;
    # VarCorr() remains the only supported route to marginal SDs.
    labs <- c(labs, if (q == 1L) paste0("logsd.", re_names[j])
                    else paste0("logchol.", re_names[j], ".", re_names[j]))
    for (r in seq_len(q)[-seq_len(j)]) {
      labs <- c(labs, paste0("chol.", re_names[r], ".", re_names[j]))
    }
  }
  labs
}


#' Fast Simplex Mixed-Effects Regression with Variable Dispersion
#'
#' @description
#' Fits a two-level (nested) simplex mixed model for continuous proportions in
#' the open interval \eqn{(0, 1)} by maximum marginal likelihood, using adaptive
#' Gauss-Hermite quadrature (AGHQ). Conditional on a cluster random effect
#' \eqn{b_j \sim N_q(0, \Sigma)},
#' \deqn{y_{ij} \mid b_j \sim \mathrm{Simplex}(\mu_{ij}, \phi_{ij}), \quad
#'       g(\mu_{ij}) = x_{ij}^\top\beta + z_{ij}^\top b_j, \quad
#'       \log\phi_{ij} = w_{ij}^\top\gamma.}
#' The fixed-effects mean and dispersion submodels use the same multi-part
#' \pkg{Formula} interface as [fastsimplexreg()]; the random effects and the
#' grouping factor are given through `random`.
#'
#' @details
#' The marginal likelihood integrates the cluster random effects out with AGHQ
#' (`nAGQ` points per dimension). `nAGQ = 1` is the Laplace approximation and
#' is accepted only with a warning: the analytic score is the score of the
#' exact marginal likelihood, not of the `nAGQ`-point quadrature, so at
#' `nAGQ = 1` the two disagree by about 66% and the resulting Wald intervals
#' cover 57% rather than 95%. Use `nAGQ >= 5`, and `nAGQ >= 11` when
#' reporting inference.
#' The per-cluster inner mode-finding, the quadrature and the analytic score are
#' implemented in C++ (RcppArmadillo, BLAS) and parallelised over clusters with
#' OpenMP, so the fit scales to large nested data sets. The random-effect
#' covariance \eqn{\Sigma = D D^\top} is estimated on an unconstrained
#' log-Cholesky scale, guaranteeing a positive-definite estimate.
#'
#' This version supports a single grouping factor (two-level nesting), Gaussian
#' random effects in the mean submodel, and fixed-effect (variable) dispersion.
#'
#' @param formula A multi-part formula `y ~ mean_terms | dispersion_terms`. When
#'   the dispersion part is omitted the dispersion is constant. Each part may
#'   carry its own `offset()` term, added to the linear predictor of that
#'   submodel on its own link scale, kept separate per submodel and rebuilt from `newdata`
#'   in [predict.simplex_fast_mixed()].
#' @param data A `data.frame` containing the response, covariates and grouping
#'   factor.
#' @param random A one-sided formula giving the random-effects design and the
#'   grouping factor, `~ z1 + z2 | group` (lme4-style bar). `~ 1 | group` is a
#'   random intercept; `~ 1 + x1 | group` a random intercept and slope.
#' @param link Mean link: one of `"logit"`, `"probit"`, `"cloglog"` or
#'   `"neglog"`. The dispersion uses a log link.
#' @param nAGQ Number of adaptive Gauss-Hermite quadrature points per random-
#'   effect dimension. Values below 5 are accepted but warn: the standard
#'   errors are then unreliable (see Details). `nAGQ = 1` is the Laplace
#'   approximation.
#' @param start Optional starting vector `c(beta, gamma, omega)`. When `NULL`,
#'   fast link-specific values are used.
#' @param maxit Maximum number of BFGS iterations.
#' @param rel_tol Relative objective tolerance.
#' @param grad_tol Infinity-norm gradient tolerance.
#' @param n_threads Number of OpenMP threads; the loop over clusters is
#'   parallelised. Zero uses all threads available to the backend. What
#'   parallelism buys here is governed by the work per CLUSTER, not by the number
#'   of clusters: measured on this design at 64000 observations, speed-up at 4
#'   threads was 1.2x with clusters of 4 observations, 3.0x with 32 and 3.6x with
#'   128. Eight threads is a loss for clusters smaller than about 32. If
#'   \pkg{RhpcBLASctl} is installed it is used to pin the BLAS to one thread for
#'   the duration of the fit and restore it afterwards, which was worth a further
#'   34\% at four threads; installing it is optional.
#'
#'   Results are reproducible for a fixed `n_threads`, but the per-thread
#'   accumulators are summed in thread order, so different thread counts differ
#'   by floating-point reassociation -- around `1e-16` relative on the
#'   log-likelihood.
#' @param inference Logical; compute the Hessian, covariance matrix and standard
#'   errors.
#' @param hessian_rel_step Relative step for the finite-difference Hessian.
#' @param inner_maxit Maximum iterations of the per-cluster inner solver.
#'   Must be at least 10: the AGHQ expansion is taken at the posterior mode,
#'   so a truncated inner solve expands around the wrong point.
#' @param inner_tol Convergence tolerance of the inner solver.
#' @param trace Logical; print optimiser progress.
#' @param subset Optional expression selecting a subset of observations,
#'   evaluated inside `data` as in [stats::lm()] -- for example
#'   `subset = x1 > 0`. A plain index, logical or row-name vector also works.
#' @param na.action Missing-data handler.
#' @param model,x,y Logical; store the model frame, the design matrices, and the
#'   response in the fitted object.
#'
#' @return An object of S3 class `"simplex_fast_mixed"`: a list whose main
#'   components are `coefficients` (a list with the `mean` and `dispersion`
#'   fixed-effect estimates), `par` (the full vector `c(beta, gamma, omega)`),
#'   `omega` (the packed log-Cholesky parameters), `D` (the Cholesky factor) and
#'   `Sigma` (the estimated random-effect covariance), `ranef` (the
#'   empirical-Bayes modes, a groups-by-`q` matrix) with `ranef.postvar`,
#'   `standard_errors`, `vcov`, `fitted.values` (means conditional on the
#'   random-effect modes), `dispersion.values`, `linear.predictors`, `residuals`
#'   (response residuals), `logLik` (the marginal log-likelihood), `AIC`, `BIC`,
#'   `nobs`, `ngrps`, `groups`, `nAGQ`, `q`, `df`, `df.residual`, `convergence`,
#'   `message`, `iterations`, the `offset` actually applied to each submodel and
#'   the stored `terms`/`design` metadata used for prediction. The same
#'   inference diagnostics as [fastsimplexreg()] are stored -- `vcov_rank`,
#'   `vcov_pseudo`, `vcov_eigenvalues`, `vcov_condition` and `n_saturated`.
#'
#' @references
#' Barndorff-Nielsen, O. E. and Jorgensen, B. (1991). Some parametric models on
#' the simplex. *Journal of Multivariate Analysis*, **39**(1), 106--116.
#' \doi{10.1016/0047-259X(91)90008-P}
#'
#' Pinheiro, J. C. and Bates, D. M. (1995). Approximations to the log-likelihood
#' function in the nonlinear mixed-effects model. *Journal of Computational and
#' Graphical Statistics*, **4**(1), 12--35.
#' \doi{10.1080/10618600.1995.10474663}
#'
#' @seealso [fastsimplexreg()], [ranef()], [VarCorr()]
#'
#' @examples
#' set.seed(1)
#' J <- 60; nj <- 8; n <- J * nj
#' dat <- data.frame(
#'   g  = factor(rep(seq_len(J), each = nj)),
#'   x1 = rnorm(n),
#'   z1 = rnorm(n)
#' )
#' b <- rnorm(J, 0, 0.7)[dat$g]
#' mu <- simplex_linkinv(0.3 - 0.6 * dat$x1 + b, "logit")
#' dat$y <- rsimplex(n, mu, exp(-0.4 + 0.3 * dat$z1))
#' fit <- fastsimplexregmixed(y ~ x1 | z1, random = ~ 1 | g, data = dat,
#'                            nAGQ = 7, n_threads = 1)
#' summary(fit)
#' VarCorr(fit)
#'
#' # Real data: gasoline yield with a random intercept per crude-oil batch.
#' if (requireNamespace("betareg", quietly = TRUE)) {
#'   data("GasolineYield", package = "betareg")
#'   gy <- fastsimplexregmixed(yield ~ temp, random = ~ 1 | batch,
#'                             data = GasolineYield, link = "logit", nAGQ = 15)
#'   summary(gy)
#' }
#'
#' @export
fastsimplexregmixed <- function(
    formula,
    data,
    random,
    link = c("logit", "probit", "cloglog", "neglog"),
    nAGQ = 11L,
    start = NULL,
    maxit = 300L,
    rel_tol = 1e-9,
    grad_tol = 1e-6,
    n_threads = 1L,
    inference = TRUE,
    hessian_rel_step = 1e-5,
    inner_maxit = 50L,
    inner_tol = 1e-8,
    trace = FALSE,
    subset = NULL,
    na.action = stats::na.omit,
    model = TRUE,
    x = FALSE,
    y = TRUE) {

  if (missing(random)) {
    stop("'random' must be supplied, e.g. random = ~ 1 | group.", call. = FALSE)
  }
  nAGQ_in <- nAGQ
  nAGQ <- as.integer(nAGQ)
  if (length(nAGQ) != 1L || is.na(nAGQ) || nAGQ < 1L) {
    stop("'nAGQ' must be a single positive integer.", call. = FALSE)
  }
  if (length(nAGQ_in) == 1L && is.numeric(nAGQ_in) && !is.na(nAGQ_in) &&
      nAGQ_in != nAGQ) {
    warning("'nAGQ' was truncated from ", format(nAGQ_in), " to ", nAGQ, ".",
            call. = FALSE)
  }
  # Guard against a small nAGQ. The analytic gradient is the score of the TRUE
  # marginal likelihood (Fisher's identity), not of its nAGQ-point quadrature
  # approximation, so the two only agree as nAGQ grows. Measured against
  # numDeriv on the exported objective, the relative gradient error is 65.7% at
  # nAGQ = 1, 5.5e-3 at 5, 5.0e-7 at 11 and 4.3e-9 at 21. At nAGQ = 1 that
  # propagated into a 8.9% error in the standard errors, 32% non-convergence
  # and 57.3% empirical coverage of a nominal 95% Wald interval (against 93.3%
  # from nAGQ = 5 upward). nAGQ = 1 is therefore not a supported configuration
  # for inference, only for a quick exploratory fit.
  if (nAGQ < 5L) {
    warning("nAGQ = ", nAGQ, " is below the supported minimum of 5. The ",
            "analytic gradient approximates the score of the exact marginal ",
            "likelihood, not of the ", nAGQ, "-point quadrature, so the ",
            "optimiser may stop early and the standard errors are not ",
            "reliable (measured coverage at nAGQ = 1 is 57%, not 95%). Use ",
            "nAGQ >= 5, and nAGQ >= 11 when reporting inference.",
            call. = FALSE)
  }
  inner_maxit <- as.integer(inner_maxit)
  if (length(inner_maxit) != 1L || is.na(inner_maxit) || inner_maxit < 1L) {
    stop("'inner_maxit' must be a single positive integer.", call. = FALSE)
  }
  # The per-cluster inner solver must be allowed to actually reach the
  # posterior mode: the AGHQ expansion is taken AT that mode, so a truncated
  # inner solve silently expands around the wrong point. Measured on a stress
  # design, the marginal negative log-likelihood came out as 5442.85, 3048.23
  # and 487.01 for inner_maxit 1, 2 and 3 against a correct value of -260.00,
  # with a singular Hessian in all three cases -- and no diagnostic.
  if (inner_maxit < 10L) {
    stop("'inner_maxit' must be at least 10. With fewer iterations the inner ",
         "solver need not reach the posterior mode, and the AGHQ expansion is ",
         "then taken around the wrong point, silently returning a wrong ",
         "marginal likelihood.", call. = FALSE)
  }

  link_spec <- .normalize_simplex_link(link)
  # `subset` is non-standard-evaluated, like lm()/glm(). See .simplex_eval_subset().
  subset_idx <- .simplex_eval_subset(substitute(subset), data, parent.frame())
  design <- .build_simplex_mixed_matrices(formula, random, data,
                                          subset = subset_idx, na.action = na.action)

  response <- design$y
  X <- design$X
  W <- design$W
  Z <- design$Z
  group <- design$group
  p <- ncol(X); r <- ncol(W); q <- ncol(Z)
  J <- nlevels(group)
  m <- q * (q + 1L) / 2L

  # Guard the tensor-product quadrature against a node-count explosion: the grid
  # has nAGQ^q points per cluster, which grows very fast with q.
  n_nodes <- nAGQ^q
  if (n_nodes > 1e5) {
    stop("The adaptive quadrature grid would have nAGQ^q = ", format(n_nodes),
         " nodes per cluster (q = ", q, " random effects, nAGQ = ", nAGQ,
         "). Reduce 'nAGQ' or the number of random-effect terms.", call. = FALSE)
  }
  if (q > 3L) {
    warning("q = ", q, " random-effect terms: adaptive Gauss-Hermite quadrature ",
            "is intended for q <= 3. Estimation may be slow and less accurate.",
            call. = FALSE)
  }

  # Degenerate cluster structures. Without these guards a single-level grouping
  # factor, or a design in which every cluster is a singleton, fits happily and
  # reports convergence = 0 while the variance component collapses to ~1e-8
  # with a standard error in the thousands. lme4::glmer refuses both.
  if (J < 2L) {
    stop("The grouping factor '", design$group_name, "' has ", J,
         " level(s). A random effect needs at least 2 sampled levels.",
         call. = FALSE)
  }
  grp_sizes <- tabulate(as.integer(group), nbins = J)
  if (all(grp_sizes <= 1L)) {
    warning("Every cluster of '", design$group_name, "' contains a single ",
            "observation, so the random effect is not separable from the ",
            "residual variation. The variance component will collapse towards ",
            "zero and its standard error is meaningless.", call. = FALSE)
  } else if (min(grp_sizes) < q) {
    warning(sum(grp_sizes < q), " of ", J, " cluster(s) of '",
            design$group_name, "' have fewer than q = ", q, " observations, ",
            "so their random effects are not identified by their own data ",
            "alone and rest entirely on the shrinkage towards Sigma.",
            call. = FALSE)
  }

  # Group-contiguous ordering + CSR offsets.
  gi <- as.integer(group)
  ord <- order(gi)
  inv_ord <- order(ord)
  starts <- as.integer(c(0L, cumsum(tabulate(gi, nbins = J))))

  y_ord <- response[ord]
  # Offsets are per-observation, so they follow the group-contiguous reordering
  # exactly as y and the design matrices do.
  off_mu_ord <- if (is.null(design$offset_mu)) NULL else design$offset_mu[ord]
  off_phi_ord <- if (is.null(design$offset_phi)) NULL else design$offset_phi[ord]
  X_ord <- X[ord, , drop = FALSE]
  W_ord <- W[ord, , drop = FALSE]
  Z_ord <- Z[ord, , drop = FALSE]
  storage.mode(X_ord) <- storage.mode(W_ord) <- storage.mode(Z_ord) <- "double"

  if (is.null(start)) {
    start <- .simplex_mixed_start(response, X, W, Z, link = link_spec$name)
  } else {
    start <- as.numeric(start)
  }
  if (length(start) != p + r + m || any(!is.finite(start))) {
    stop("'start' must be a finite numeric vector of length ncol(X) + ncol(W) + q(q+1)/2.",
         call. = FALSE)
  }

  # BLAS thread pinning, for the duration of this call only.
  #
  # The mixed backend parallelises over CLUSTERS and each cluster does several
  # small Armadillo products and decompositions. A threaded BLAS then layers its
  # own team under ours -- 8 OpenMP threads times up to 24 BLAS threads on a
  # 24-core box -- and the machine spends its time on scheduling. Measured on a
  # J = 2000, nj = 8, nAGQ = 11 fit with identical data (log-likelihood equal to
  # 12 digits in every run):
  #
  #   n_threads   BLAS at 24   BLAS at 1   gain
  #           1       6.10 s      5.81 s    +5%
  #           4       4.02 s      2.66 s   +34%
  #           8       7.01 s      5.40 s   +23%
  #
  # The gap widens with the OpenMP team, which is the signature of
  # oversubscription. Restored on exit, so the user's global BLAS setting is
  # never left changed -- and skipped entirely when RhpcBLASctl is absent, which
  # only costs the speed-up.
  if (n_threads != 1L && requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    .blas_old <- try(RhpcBLASctl::blas_get_num_procs(), silent = TRUE)
    if (!inherits(.blas_old, "try-error") && is.numeric(.blas_old)) {
      on.exit(try(RhpcBLASctl::blas_set_num_threads(.blas_old), silent = TRUE),
              add = TRUE)
      try(RhpcBLASctl::blas_set_num_threads(1L), silent = TRUE)
    }
  }

  opt <- simplex_mixed_bfgs_cpp(
    start = start, y = y_ord, X = X_ord, Z = Z_ord, W = W_ord,
    starts = starts, q = as.integer(q), mean_link = link_spec$id,
    nAGQ = nAGQ, maxit = as.integer(maxit), rel_tol = as.numeric(rel_tol),
    grad_tol = as.numeric(grad_tol), n_threads = as.integer(n_threads),
    inner_maxit = inner_maxit, inner_tol = as.numeric(inner_tol),
    trace = isTRUE(trace), off_mu_ = off_mu_ord, off_phi_ = off_phi_ord
  )

  theta <- as.numeric(opt$par)
  beta <- theta[seq_len(p)]
  gamma <- theta[p + seq_len(r)]
  omega <- theta[p + r + seq_len(m)]
  names(beta) <- colnames(X)
  names(gamma) <- colnames(W)

  re_names <- colnames(Z)
  D <- simplex_mixed_D_from_omega_cpp(omega, as.integer(q))
  Sigma <- D %*% t(D)
  dimnames(Sigma) <- list(re_names, re_names)

  # Random-effect predictions (empirical Bayes modes + posterior covariances).
  re <- simplex_mixed_ranef_cpp(theta, y_ord, X_ord, Z_ord, W_ord, starts,
                                as.integer(q), link_spec$id, as.integer(n_threads),
                                inner_maxit, as.numeric(inner_tol),
                                off_mu_ = off_mu_ord, off_phi_ = off_phi_ord)
  ranef_mat <- re$b
  dimnames(ranef_mat) <- list(levels(group), re_names)

  # Conditional fitted values (include random effects), mapped to original order.
  pred <- simplex_mixed_predict_cpp(theta, X_ord, Z_ord, W_ord, starts,
                                    as.integer(q), re$b, link_spec$id, TRUE,
                                    off_mu_ = off_mu_ord, off_phi_ = off_phi_ord)
  mu_ord <- as.numeric(pred$mu)
  phi_ord <- as.numeric(pred$phi)
  eta_mu_ord <- as.numeric(pred$eta_mu)
  eta_phi_ord <- as.numeric(pred$eta_phi)
  mu <- mu_ord[inv_ord]
  phi <- phi_ord[inv_ord]
  eta_mu <- eta_mu_ord[inv_ord]
  eta_phi <- eta_phi_ord[inv_ord]

  logLik_value <- -as.numeric(opt$value)
  k <- length(theta)
  n <- length(response)

  par_names <- c(.simplex_par_names(colnames(X), colnames(W)),
                 .omega_labels(re_names))
  names(theta) <- par_names

  converged <- as.integer(opt$convergence) == 0L
  if (!converged) {
    warning("fastsimplexregmixed() did not converge (code ", opt$convergence, ": ",
            opt$message, "). Estimates and standard errors are unreliable; try a ",
            "larger 'nAGQ' or different starting values.", call. = FALSE)
  }

  vc <- NULL
  se <- stats::setNames(rep(NA_real_, k), par_names)
  hessian <- NULL
  vcov_rank <- NA_integer_
  vcov_pseudo <- NA
  vcov_eigenvalues <- NULL
  vcov_condition <- NA_real_
  # Standard errors only at a converged fit (see fastsimplexreg()).
  if (isTRUE(inference) && converged) {
    hessian <- simplex_mixed_hessian_fd_cpp(
      theta = theta, y = y_ord, X = X_ord, Z = Z_ord, W = W_ord,
      starts = starts, q = as.integer(q), mean_link = link_spec$id, nAGQ = nAGQ,
      rel_step = as.numeric(hessian_rel_step), n_threads = as.integer(n_threads),
      inner_maxit = inner_maxit, inner_tol = as.numeric(inner_tol),
      off_mu_ = off_mu_ord, off_phi_ = off_phi_ord
    )
    dimnames(hessian) <- list(par_names, par_names)

    # Same fail-safe path as the fixed-effects fit (see R/inference.R).
    inf <- .simplex_vcov(hessian, par_names, what = "fastsimplexregmixed()")
    vc <- inf$vcov
    se <- inf$se
    vcov_rank <- inf$rank
    vcov_pseudo <- inf$pseudo
    vcov_eigenvalues <- inf$eigenvalues
    vcov_condition <- inf$condition
  }

  .warn_saturated(opt$n_saturated, n, what = "fastsimplexregmixed()")

  obs_names <- rownames(design$model)
  na_act <- attr(design$model, "na.action")

  out <- list(
    call = match.call(),
    formula = design$formula,
    random = random,
    link = list(mean = link_spec$name, dispersion = "log"),
    coefficients = list(
      mean = stats::setNames(beta, colnames(X)),
      dispersion = stats::setNames(gamma, colnames(W))
    ),
    par = theta,
    omega = stats::setNames(omega, .omega_labels(re_names)),
    D = D,
    Sigma = Sigma,
    ranef = ranef_mat,
    ranef.postvar = re$postvar,
    standard_errors = stats::setNames(se, par_names),
    vcov = vc,
    vcov_rank = vcov_rank,
    vcov_pseudo = vcov_pseudo,
    vcov_eigenvalues = vcov_eigenvalues,
    vcov_condition = vcov_condition,
    n_saturated = as.integer(opt$n_saturated),
    hessian = hessian,
    na.action = na_act,
    fitted.values = stats::setNames(mu, obs_names),
    dispersion.values = stats::setNames(phi, obs_names),
    linear.predictors = list(mean = stats::setNames(eta_mu, obs_names),
                             dispersion = stats::setNames(eta_phi, obs_names)),
    residuals = stats::setNames(response - mu, obs_names),
    logLik = logLik_value,
    AIC = -2 * logLik_value + 2 * k,
    BIC = -2 * logLik_value + log(n) * k,
    nobs = n,
    ngrps = J,
    groups = group,
    group_name = design$group_name,
    df = k,
    df.residual = n - k,
    nAGQ = nAGQ,
    q = q,
    convergence = as.integer(opt$convergence),
    message = as.character(opt$message),
    iterations = as.integer(opt$iterations),
    function_evaluations = as.integer(opt$function_evaluations),
    gradient_evaluations = as.integer(opt$gradient_evaluations),
    gradient = as.numeric(opt$gradient),
    offset = list(mean = design$offset_mu, dispersion = design$offset_phi),
    terms = list(mean = design$terms_mean, dispersion = design$terms_dispersion,
                 random = design$terms_random),
    design = design[c(
      "terms_mean", "terms_dispersion", "terms_random",
      "xlevels_mean", "xlevels_dispersion", "xlevels_random",
      "contrasts_mean", "contrasts_dispersion", "contrasts_random",
      "group_name", "has_dispersion_formula"
    )],
    order = list(ord = ord, inv_ord = inv_ord, starts = starts),
    n_threads = as.integer(n_threads)
  )

  if (isTRUE(model)) out$model <- design$model
  if (isTRUE(x)) out$x <- list(mean = X, dispersion = W, random = Z)
  if (isTRUE(y)) out$y <- response

  structure(out, class = "simplex_fast_mixed")
}
