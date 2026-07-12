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

  # Apply the subset explicitly (the NSE-safe way; never forward to model.frame).
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

  group <- droplevels(as.factor(grp))

  list(
    formula = Fu,
    random = random,
    combined = Fc,
    model = mf,
    y = response,
    X = X, W = W, Z = Z,
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
  for (c in seq_len(q)) {
    labs <- c(labs, paste0("logsd.", re_names[c]))
    for (r in seq_len(q)[-seq_len(c)]) {
      labs <- c(labs, paste0("chol.", re_names[r], ".", re_names[c]))
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
#' (`nAGQ` points per dimension; `nAGQ = 1` gives the Laplace approximation).
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
#'   the dispersion part is omitted the dispersion is constant.
#' @param data A `data.frame` containing the response, covariates and grouping
#'   factor.
#' @param random A one-sided formula giving the random-effects design and the
#'   grouping factor, `~ z1 + z2 | group` (lme4-style bar). `~ 1 | group` is a
#'   random intercept; `~ 1 + x1 | group` a random intercept and slope.
#' @param link Mean link: one of `"logit"`, `"probit"`, `"cloglog"` or
#'   `"neglog"`. The dispersion uses a log link.
#' @param nAGQ Number of adaptive Gauss-Hermite quadrature points per random-
#'   effect dimension. `nAGQ = 1` is the Laplace approximation.
#' @param start Optional starting vector `c(beta, gamma, omega)`. When `NULL`,
#'   fast link-specific values are used.
#' @param maxit Maximum number of BFGS iterations.
#' @param rel_tol Relative objective tolerance.
#' @param grad_tol Infinity-norm gradient tolerance.
#' @param n_threads Number of OpenMP threads (the cluster loop is parallelised).
#'   Zero uses all available threads. Parallelism helps most when the per-cluster
#'   work is substantial (two or more random effects, or larger clusters); for
#'   many tiny clusters a small `n_threads` (or `1`) can be faster, because a
#'   multi-threaded BLAS may otherwise oversubscribe the cores. Results can differ
#'   by a negligible amount (around `1e-13`) between thread counts.
#' @param inference Logical; compute the Hessian, covariance matrix and standard
#'   errors.
#' @param hessian_rel_step Relative step for the finite-difference Hessian.
#' @param inner_maxit Maximum iterations of the per-cluster inner solver.
#' @param inner_tol Convergence tolerance of the inner solver.
#' @param trace Logical; print optimiser progress.
#' @param subset Optional index vector selecting observations.
#' @param na.action Missing-data handler.
#' @param model,x,y Logical; store the model frame, the design matrices, and the
#'   response in the fitted object.
#'
#' @return An object of class `"simplex_fast_mixed"`.
#'
#' @references
#' Barndorff-Nielsen, O. E. and Jorgensen, B. (1991). Some parametric models on
#' the simplex. *Journal of Multivariate Analysis*, **39**(1), 106-116.
#'
#' Pinheiro, J. C. and Bates, D. M. (1995). Approximations to the log-likelihood
#' function in the nonlinear mixed-effects model. *Journal of Computational and
#' Graphical Statistics*, **4**(1), 12-35.
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
  nAGQ <- as.integer(nAGQ)
  if (length(nAGQ) != 1L || is.na(nAGQ) || nAGQ < 1L) {
    stop("'nAGQ' must be a single positive integer.", call. = FALSE)
  }

  link_spec <- .normalize_simplex_link(link)
  design <- .build_simplex_mixed_matrices(formula, random, data,
                                          subset = subset, na.action = na.action)

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

  # Group-contiguous ordering + CSR offsets.
  gi <- as.integer(group)
  ord <- order(gi)
  inv_ord <- order(ord)
  starts <- as.integer(c(0L, cumsum(tabulate(gi, nbins = J))))

  y_ord <- response[ord]
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

  opt <- simplex_mixed_bfgs_cpp(
    start = start, y = y_ord, X = X_ord, Z = Z_ord, W = W_ord,
    starts = starts, q = as.integer(q), mean_link = link_spec$id,
    nAGQ = nAGQ, maxit = as.integer(maxit), rel_tol = as.numeric(rel_tol),
    grad_tol = as.numeric(grad_tol), n_threads = as.integer(n_threads),
    inner_maxit = as.integer(inner_maxit), inner_tol = as.numeric(inner_tol),
    trace = isTRUE(trace)
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
                                as.integer(inner_maxit), as.numeric(inner_tol))
  ranef_mat <- re$b
  dimnames(ranef_mat) <- list(levels(group), re_names)

  # Conditional fitted values (include random effects), mapped to original order.
  pred <- simplex_mixed_predict_cpp(theta, X_ord, Z_ord, W_ord, starts,
                                    as.integer(q), re$b, link_spec$id, TRUE)
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

  par_names <- c(colnames(X), colnames(W), .omega_labels(re_names))
  names(theta) <- par_names

  converged <- as.integer(opt$convergence) == 0L
  if (!converged) {
    warning("fastsimplexregmixed() did not converge (code ", opt$convergence, ": ",
            opt$message, "). Estimates and standard errors are unreliable; try a ",
            "larger 'nAGQ' or different starting values.", call. = FALSE)
  }

  vc <- NULL
  se <- rep(NA_real_, k)
  hessian <- NULL
  # Standard errors only at a converged fit (see fastsimplexreg()).
  if (isTRUE(inference) && converged) {
    hessian <- simplex_mixed_hessian_fd_cpp(
      theta = theta, y = y_ord, X = X_ord, Z = Z_ord, W = W_ord,
      starts = starts, q = as.integer(q), mean_link = link_spec$id, nAGQ = nAGQ,
      rel_step = as.numeric(hessian_rel_step), n_threads = as.integer(n_threads),
      inner_maxit = as.integer(inner_maxit), inner_tol = as.numeric(inner_tol)
    )
    vc <- tryCatch(
      solve(hessian),
      error = function(e) qr.solve(hessian, diag(nrow(hessian)), tol = 1e-10)
    )
    vc <- 0.5 * (vc + t(vc))
    dimnames(vc) <- list(par_names, par_names)
    se <- sqrt(pmax(diag(vc), 0))
  }

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
    hessian = hessian,
    fitted.values = mu,
    dispersion.values = phi,
    linear.predictors = list(mean = eta_mu, dispersion = eta_phi),
    residuals = response - mu,
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
