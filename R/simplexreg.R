# Internal helper: build one shared model frame so that subset handling, NA
# handling, factors and contrasts are perfectly aligned across the mean and
# dispersion components. Uses the Formula package for multi-part formulas.
.build_simplex_matrices <- function(formula, data, subset = NULL, na.action = stats::na.omit) {
  if (!requireNamespace("Formula", quietly = TRUE)) {
    stop("Package 'Formula' is required. Install it with install.packages('Formula').", call. = FALSE)
  }
  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula such as y ~ x1 + x2 | z1 + z2.", call. = FALSE)
  }

  F <- Formula::Formula(formula)
  dims <- length(F)

  if (dims[1L] != 1L) {
    stop("The model must contain exactly one response component.", call. = FALSE)
  }
  if (dims[2L] < 1L || dims[2L] > 2L) {
    stop("Use one or two RHS components: y ~ mean_terms | dispersion_terms.", call. = FALSE)
  }

  # Apply the observation subset explicitly. Forwarding 'subset' straight to
  # model.frame() would trigger its non-standard evaluation, which resolves the
  # bare symbol in the formula's environment (picking up base::subset) instead
  # of using the value supplied here. 'subset' is documented as an index vector.
  if (!is.null(subset)) {
    data <- data[subset, , drop = FALSE]
  }

  mf <- stats::model.frame(
    F,
    data = data,
    na.action = na.action,
    drop.unused.levels = TRUE
  )

  response <- Formula::model.part(F, data = mf, lhs = 1L, drop = TRUE)
  X <- stats::model.matrix(F, data = mf, rhs = 1L)

  if (dims[2L] == 2L) {
    Z <- stats::model.matrix(F, data = mf, rhs = 2L)
    terms_dispersion <- stats::terms(F, rhs = 2L)
  } else {
    Z <- matrix(
      1.0,
      nrow = nrow(X),
      ncol = 1L,
      dimnames = list(rownames(X), "(Intercept)")
    )
    terms_dispersion <- stats::terms(~ 1)
  }

  if (!is.numeric(response)) {
    stop("The response must be numeric.", call. = FALSE)
  }
  response <- as.numeric(response)
  if (any(!is.finite(response)) || any(response <= 0 | response >= 1)) {
    stop("All response values must be finite and strictly inside (0, 1).", call. = FALSE)
  }

  storage.mode(X) <- "double"
  storage.mode(Z) <- "double"

  terms_mean <- stats::terms(F, rhs = 1L)

  list(
    formula = F,
    model = mf,
    y = response,
    X = X,
    Z = Z,
    terms_mean = terms_mean,
    terms_dispersion = terms_dispersion,
    xlevels_mean = stats::.getXlevels(terms_mean, mf),
    xlevels_dispersion = if (dims[2L] == 2L) {
      stats::.getXlevels(terms_dispersion, mf)
    } else {
      list()
    },
    contrasts_mean = attr(X, "contrasts"),
    contrasts_dispersion = attr(Z, "contrasts"),
    has_dispersion_formula = dims[2L] == 2L
  )
}


# Internal helper: stable, link-specific starting values c(beta, gamma).
.simplex_start <- function(y, X, Z, link) {
  p <- ncol(X)
  q <- ncol(Z)
  beta <- numeric(p)
  gamma <- numeric(q)

  mu0 <- min(1 - 1e-6, max(1e-6, mean(y)))
  eta0 <- switch(
    link,
    logit = stats::qlogis(mu0),
    probit = stats::qnorm(mu0),
    cloglog = log(-log1p(-mu0)),
    neglog = -log(-log(mu0)),
    stop("Unsupported link.", call. = FALSE)
  )

  intercept_x <- which(colnames(X) == "(Intercept)")
  if (length(intercept_x)) {
    beta[intercept_x[1L]] <- eta0
  }

  qmu <- mu0 * (1 - mu0)
  dev0 <- (y - mu0)^2 / (y * (1 - y) * qmu^2)
  phi0 <- max(mean(dev0), 1e-6)
  intercept_z <- which(colnames(Z) == "(Intercept)")
  if (length(intercept_z)) {
    gamma[intercept_z[1L]] <- log(phi0)
  }

  c(beta, gamma)
}


#' @title Fit a Fast Simplex Regression with Variable Dispersion
#'
#' @description
#' Fits, by maximum likelihood, a simplex regression model with separate
#' submodels for the mean and the dispersion. The interface uses the multi-part
#' formulas of the \pkg{Formula} package:
#'
#' `y ~ x1 + x2 | z1 + z2`
#'
#' The first right-hand side component models the mean \eqn{\mu}; the second
#' component models the dispersion \eqn{\phi}. When the second component is
#' omitted, as in `y ~ x1 + x2`, the dispersion is constant (equivalent to
#' `| 1`).
#'
#' The mean supports the `logit`, `probit`, `cloglog` and `neglog` links; the
#' dispersion uses a log link. The log-likelihood, the analytic score, the link
#' inverses and the BFGS optimiser run entirely in C++. Matrix-vector products
#' use Armadillo/BLAS and the per-observation loop may use OpenMP.
#'
#' @param formula A multi-part formula, for example `y ~ x1 + x2 | z1 + z2`.
#' @param data A `data.frame` containing the response and covariates.
#' @param link Character string selecting the mean link: `"logit"`, `"probit"`,
#'   `"cloglog"` or `"neglog"`.
#' @param start Optional numeric starting vector `c(beta, gamma)`. When `NULL`,
#'   fast link-specific starting values are used.
#' @param maxit Integer; the maximum number of BFGS iterations.
#' @param rel_tol Numeric; relative tolerance on the objective function.
#' @param grad_tol Numeric; tolerance on the infinity norm of the gradient.
#' @param n_threads Integer number of OpenMP threads. Use `0` to request all
#'   threads available to the backend.
#' @param inference Logical; if `TRUE`, computes the Hessian, the
#'   variance-covariance matrix and the standard errors.
#' @param hessian_rel_step Numeric; the initial relative step for the Hessian,
#'   obtained by central differences of the analytic gradient.
#' @param trace Logical; if `TRUE`, prints optimiser progress.
#' @param subset Optional vector specifying a subset of observations.
#' @param na.action A function indicating how to handle missing values.
#' @param model Logical; if `TRUE`, stores the model frame in the fitted object.
#' @param x Logical; if `TRUE`, stores the design matrices `X` and `Z`.
#' @param y Logical; if `TRUE`, stores the response in the fitted object.
#'
#' @return An object of S3 class `"simplex_fast"`: a list whose main components
#'   are `coefficients` (a list with `mean` and `dispersion` estimates), `par`
#'   (the full coefficient vector), `standard_errors`, `vcov`, `fitted.values`
#'   (fitted means), `dispersion.values` (fitted dispersions),
#'   `linear.predictors`, `residuals` (response residuals), `logLik`, `AIC`,
#'   `BIC`, `nobs`, `df.residual`, `convergence`, `message`, `iterations` and
#'   the stored `terms`/`design` metadata used for prediction.
#'
#' @references
#' Barndorff-Nielsen, O. E. and Jorgensen, B. (1991).
#' Some parametric models on the simplex.
#' *Journal of Multivariate Analysis*, **39**(1), 106--116.
#'
#' Zhang, P., Qiu, Z. and Shi, C. (2016).
#' simplexreg: An R Package for Regression Analysis of Proportional Data Using
#' the Simplex Distribution.
#' *Journal of Statistical Software*, **71**(11), 1--21.
#'
#' @seealso [dsimplex()], [rsimplex()], [simplex_linkinv()],
#'   [predict.simplex_fast()], [summary.simplex_fast()]
#'
#' @examples
#' # Simulated data with variable dispersion.
#' set.seed(123)
#' n <- 500
#' dat <- data.frame(x1 = rnorm(n), x2 = rbinom(n, 1, 0.4), z1 = rnorm(n))
#' mu <- simplex_linkinv(-0.4 + 0.8 * dat$x1 - 0.5 * dat$x2, link = "logit")
#' phi <- exp(-1 + 0.6 * dat$z1)
#' dat$y <- rsimplex(n, mu, phi)
#'
#' fit <- fastsimplexreg(y ~ x1 + x2 | z1, data = dat, link = "logit",
#'                    n_threads = 1L)
#' summary(fit)
#' coef(fit)
#' head(predict(fit, type = "both"))
#'
#' # Real data: reading accuracy from the 'betareg' package.
#' if (requireNamespace("betareg", quietly = TRUE)) {
#'   data("ReadingSkills", package = "betareg")
#'   rs <- fastsimplexreg(accuracy ~ dyslexia + iq | dyslexia,
#'                        data = ReadingSkills, link = "logit")
#'   summary(rs)
#' }
#'
#' @export
fastsimplexreg <- function(
    formula,
    data,
    link = c("logit", "probit", "cloglog", "neglog"),
    start = NULL,
    maxit = 300L,
    rel_tol = 1e-9,
    grad_tol = 1e-6,
    n_threads = 1L,
    inference = TRUE,
    hessian_rel_step = 1e-5,
    trace = FALSE,
    subset = NULL,
    na.action = stats::na.omit,
    model = TRUE,
    x = FALSE,
    y = TRUE) {

  link_spec <- .normalize_simplex_link(link)
  design <- .build_simplex_matrices(
    formula = formula,
    data = data,
    subset = subset,
    na.action = na.action
  )

  response <- design$y
  X <- design$X
  Z <- design$Z
  p <- ncol(X)
  q <- ncol(Z)

  if (is.null(start)) {
    start <- .simplex_start(response, X, Z, link = link_spec$name)
  } else {
    start <- as.numeric(start)
  }
  if (length(start) != p + q || any(!is.finite(start))) {
    stop("'start' must be a finite numeric vector of length ncol(X) + ncol(Z).", call. = FALSE)
  }

  opt <- simplex_bfgs_cpp(
    start = start,
    y = response,
    X = X,
    Z = Z,
    mean_link = link_spec$id,
    maxit = as.integer(maxit),
    rel_tol = as.numeric(rel_tol),
    grad_tol = as.numeric(grad_tol),
    n_threads = as.integer(n_threads),
    trace = isTRUE(trace)
  )

  theta <- as.numeric(opt$par)
  # Bare coefficient names, matching the convention of other simplex/beta
  # regression packages. The mean and dispersion submodels are distinguished by
  # position (the first p entries are the mean coefficients, the remaining q are
  # the dispersion coefficients) rather than by a name prefix.
  names(theta) <- c(colnames(X), colnames(Z))

  pred <- simplex_predict_cpp(theta, X, Z, mean_link = link_spec$id)
  logLik_value <- -as.numeric(opt$value)
  k <- length(theta)
  n <- length(response)

  converged <- as.integer(opt$convergence) == 0L
  if (!converged) {
    warning("fastsimplexreg() did not converge (code ", opt$convergence, ": ",
            opt$message, "). Estimates and standard errors are unreliable.",
            call. = FALSE)
  }

  vc <- NULL
  se <- rep(NA_real_, k)
  hessian <- NULL
  # Standard errors are only computed at a converged (stationary) fit; at a
  # non-converged point the Hessian is meaningless, so leave them NA.
  if (isTRUE(inference) && converged) {
    hessian <- simplex_hessian_fd_cpp(
      theta = theta,
      y = response,
      X = X,
      Z = Z,
      mean_link = link_spec$id,
      rel_step = as.numeric(hessian_rel_step),
      n_threads = as.integer(n_threads)
    )

    vc <- tryCatch(
      solve(hessian),
      error = function(e) qr.solve(hessian, diag(nrow(hessian)), tol = 1e-10)
    )
    vc <- 0.5 * (vc + t(vc))
    dimnames(vc) <- list(names(theta), names(theta))
    se <- sqrt(pmax(diag(vc), 0))
  }

  coefficients <- list(
    mean = stats::setNames(theta[seq_len(p)], colnames(X)),
    dispersion = stats::setNames(theta[p + seq_len(q)], colnames(Z))
  )

  out <- list(
    call = match.call(),
    formula = design$formula,
    link = list(mean = link_spec$name, dispersion = "log"),
    coefficients = coefficients,
    par = theta,
    standard_errors = stats::setNames(se, names(theta)),
    vcov = vc,
    hessian = hessian,
    fitted.values = as.numeric(pred$mu),
    dispersion.values = as.numeric(pred$phi),
    linear.predictors = list(
      mean = as.numeric(pred$eta_mu),
      dispersion = as.numeric(pred$eta_phi)
    ),
    residuals = response - as.numeric(pred$mu),
    logLik = logLik_value,
    AIC = -2 * logLik_value + 2 * k,
    BIC = -2 * logLik_value + log(n) * k,
    nobs = n,
    df.residual = n - k,
    convergence = as.integer(opt$convergence),
    message = as.character(opt$message),
    iterations = as.integer(opt$iterations),
    function_evaluations = as.integer(opt$function_evaluations),
    gradient_evaluations = as.integer(opt$gradient_evaluations),
    gradient = as.numeric(opt$gradient),
    terms = list(mean = design$terms_mean, dispersion = design$terms_dispersion),
    design = design[c(
      "terms_mean",
      "terms_dispersion",
      "xlevels_mean",
      "xlevels_dispersion",
      "contrasts_mean",
      "contrasts_dispersion",
      "has_dispersion_formula"
    )],
    n_threads = as.integer(n_threads)
  )

  if (isTRUE(model)) out$model <- design$model
  if (isTRUE(x)) out$x <- list(mean = X, dispersion = Z)
  if (isTRUE(y)) out$y <- response

  structure(out, class = "simplex_fast")
}
