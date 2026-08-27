# Internal: resolve the `subset` argument the way stats::lm/glm do -- as an
# EXPRESSION evaluated inside `data` first and only then in the caller.
#
# Until 0.2.4 `subset` was documented like glm's ("a vector specifying a subset
# of observations") but implemented as a plain index vector: `data[subset, ]`.
# The difference is silent and dangerous. `subset = x1 > 0` errored loudly when
# no `x1` existed in the caller -- but when an unrelated object of that name did
# exist, which is the common case inside a function, the fit used THAT vector
# and quietly estimated the model on the wrong rows.
#
# `expr` is the result of substitute(subset) in the calling fitter; `envir` is
# its parent.frame(). Every call form that worked before still works: a plain
# index/logical/character vector evaluates to itself.
.simplex_eval_subset <- function(expr, data, envir) {
  if (is.null(expr)) return(NULL)
  r <- eval(expr, data, envir)
  if (is.null(r)) return(NULL)
  if (is.logical(r)) {
    if (length(r) != nrow(data)) {
      stop("A logical 'subset' must have length nrow(data) (", nrow(data),
           "), not ", length(r), ".", call. = FALSE)
    }
    # Match stats::model.frame: an NA in a logical subset drops the row.
    return(r & !is.na(r))
  }
  r
}


# Internal: evaluate the offset(s) of ONE model part.
#
# stats::model.offset() cannot be used here: with a multi-part Formula it sums
# the offsets of every part, so `y ~ x + offset(a) | z + offset(b)` would give
# a + b for both the mean and the dispersion. Each part's own `terms` object
# carries an "offset" attribute indexing into its "variables", which is what
# keeps the two apart.
#
# Works from either a model frame (where the offset already exists as a column
# named e.g. "offset(a)") or a raw `newdata` (where the expression is
# evaluated). Returns NULL when the part has no offset, else a numeric vector.
.simplex_offset <- function(tt, data) {
  idx <- attr(tt, "offset")
  if (is.null(idx) || !length(idx)) return(NULL)
  vars <- attr(tt, "variables")
  env <- environment(tt)
  if (is.null(env)) env <- parent.frame()
  vals <- lapply(idx, function(i) {
    call_i <- vars[[i + 1L]]
    nm <- deparse(call_i, width.cutoff = 500L)
    v <- if (!is.null(data[[nm]])) data[[nm]] else eval(call_i, data, env)
    as.numeric(v)
  })
  Reduce(`+`, vals)
}


# Internal: validate an offset against the number of rows it must align with.
.simplex_check_offset <- function(off, n, what) {
  if (is.null(off)) return(NULL)
  if (length(off) == 1L) off <- rep(off, n)
  if (length(off) != n) {
    stop("The ", what, " offset has length ", length(off), " but the model has ",
         n, " observation(s).", call. = FALSE)
  }
  if (any(!is.finite(off))) {
    stop("The ", what, " offset contains non-finite values.", call. = FALSE)
  }
  off
}


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

  fml <- Formula::Formula(formula)
  dims <- length(fml)

  if (dims[1L] != 1L) {
    stop("The model must contain exactly one response component.", call. = FALSE)
  }
  if (dims[2L] < 1L || dims[2L] > 2L) {
    stop("Use one or two RHS components: y ~ mean_terms | dispersion_terms.", call. = FALSE)
  }

  # `subset` arrives here ALREADY RESOLVED to a plain index/logical vector by
  # .simplex_eval_subset(), which the fitter calls with substitute(subset) so
  # that an expression like `x1 > 0` is evaluated inside `data`. It is applied
  # here rather than forwarded to model.frame(), whose own non-standard
  # evaluation would resolve the bare symbol in the formula's environment
  # (picking up base::subset) instead of using the value supplied here.
  if (!is.null(subset)) {
    data <- data[subset, , drop = FALSE]
  }

  mf <- stats::model.frame(
    fml,
    data = data,
    na.action = na.action,
    drop.unused.levels = TRUE
  )

  response <- Formula::model.part(fml, data = mf, lhs = 1L, drop = TRUE)
  X <- stats::model.matrix(fml, data = mf, rhs = 1L)

  if (dims[2L] == 2L) {
    Z <- stats::model.matrix(fml, data = mf, rhs = 2L)
    terms_dispersion <- stats::terms(fml, rhs = 2L)
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

  terms_mean <- stats::terms(fml, rhs = 1L)

  n_obs <- nrow(X)
  offset_mu <- .simplex_check_offset(.simplex_offset(terms_mean, mf), n_obs, "mean")
  offset_phi <- .simplex_check_offset(.simplex_offset(terms_dispersion, mf),
                                      n_obs, "dispersion")

  list(
    formula = fml,
    model = mf,
    y = response,
    X = X,
    Z = Z,
    offset_mu = offset_mu,
    offset_phi = offset_phi,
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


# Internal helper: build the prediction design matrices for `newdata` from the
# stored terms/contrasts/xlevels of a fit.
#
# Two things this guarantees that a bare model.matrix() call does not:
#   1. ALL parts (mean, dispersion, and for the mixed fit the random design)
#      are built from the SAME set of complete rows, so they can never end up
#      with different row counts when the missing values sit in different
#      columns;
#   2. the rows dropped for missingness are recorded, so the caller can re-expand
#      the predictions back to nrow(newdata) with NA in the dropped positions --
#      the convention of stats::predict.lm(), instead of silently returning a
#      shorter, unaligned vector.
# `extra_vars` names columns that must be PRESENT but are not part of the
# completeness test (the grouping factor of a mixed fit: a missing level means
# "no known cluster", which the unseen-level path handles, not a missing
# covariate).
# Returns a list with `matrices` (in the order of `terms_list`), `keep` (logical
# over the rows of newdata) and `n` (nrow(newdata)).
.simplex_predict_design <- function(newdata, terms_list, contrasts_list,
                                    xlev_list, extra_vars = character(0)) {
  if (!is.data.frame(newdata)) newdata <- as.data.frame(newdata)
  tms <- lapply(terms_list, stats::delete.response)
  model_vars <- unique(unlist(lapply(tms, all.vars)))
  absent <- setdiff(unique(c(model_vars, extra_vars)), names(newdata))
  if (length(absent)) {
    stop("Variable(s) required by the model are missing from 'newdata': ",
         paste(absent, collapse = ", "), ".", call. = FALSE)
  }

  n <- nrow(newdata)
  keep <- if (length(model_vars)) {
    stats::complete.cases(newdata[, model_vars, drop = FALSE])
  } else {
    rep(TRUE, n)
  }
  nd <- newdata[keep, , drop = FALSE]

  matrices <- lapply(seq_along(tms), function(i) {
    mm <- stats::model.matrix(tms[[i]], data = nd,
                              contrasts.arg = contrasts_list[[i]],
                              xlev = xlev_list[[i]])
    storage.mode(mm) <- "double"
    mm
  })
  names(matrices) <- names(terms_list)
  list(matrices = matrices, keep = keep, n = n)
}


# Internal helper: re-expand a prediction computed on the complete rows back to
# the full length of `newdata`, with NA where a row was dropped.
.simplex_expand <- function(values, keep) {
  if (all(keep)) return(values)
  out <- rep(NA_real_, length(keep))
  out[keep] <- values
  out
}


# Internal: detect columns of a design matrix that are aliased (linearly
# dependent on earlier columns), using the same pivoted QR that lm() uses.
#
# Without this the optimiser happily converges on a rank-deficient design and
# reports a finite estimate for EVERY column, splitting one identified effect
# arbitrarily across the collinear group. With x2 = 2 * x1 it returned
# x1 = 0.105 and x2 = 0.210, two numbers that mean nothing individually --
# only x1 + 2*x2 = 0.524 is identified, which is exactly the coefficient the
# reduced model gives. lm()/glm() report NA for the aliased column instead.
#
# Returns a logical vector over the columns, TRUE where aliased.
.simplex_aliased <- function(M) {
  if (ncol(M) == 0L) return(logical(0))
  qrM <- qr(M, tol = 1e-7, LAPACK = FALSE)
  aliased <- rep(TRUE, ncol(M))
  if (qrM$rank > 0L) aliased[qrM$pivot[seq_len(qrM$rank)]] <- FALSE
  aliased
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
#' @param inference Logical; if `TRUE`, computes the information matrix, the
#'   variance-covariance matrix and the standard errors.
#' @param information Character; which information matrix to invert for the
#'   standard errors. `"observed"` (default) uses the observed information, the
#'   Hessian of the negative log-likelihood obtained by central differences of
#'   the analytic score. `"expected"` uses the exact Fisher information, which
#'   for the simplex is available in closed form and is block diagonal in
#'   \eqn{(\beta, \gamma)}: it needs no finite differencing, is positive
#'   definite by construction, and is roughly twenty times cheaper. The two
#'   agree asymptotically and, at \eqn{n = 4000}, to within 0.4\%. The default
#'   stays `"observed"` because Efron and Hinkley (1978) argue it is the better
#'   variance estimator for conditional inference; `"expected"` is the more
#'   robust choice when the observed information is ill-conditioned.
#' @param hessian_rel_step Numeric; the initial relative step for the Hessian,
#'   obtained by central differences of the analytic gradient.
#' @param trace Logical; if `TRUE`, prints optimiser progress.
#' @param subset Optional expression selecting a subset of observations,
#'   evaluated inside `data` as in [stats::lm()] -- for example
#'   `subset = x1 > 0`. A plain index, logical or row-name vector also
#'   works. An `NA` in a logical subset drops that row.
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
#'   the stored `terms`/`design` metadata used for prediction. Inference
#'   diagnostics are also stored: `vcov_rank` (rank of the observed information
#'   matrix), `vcov_pseudo` (`TRUE` when a Moore-Penrose pseudo-inverse was
#'   required because the matrix was rank deficient or indefinite),
#'   `vcov_eigenvalues`, and `n_saturated` (observations whose fitted mean hit
#'   the numerical boundary of the likelihood path). Standard errors of
#'   parameters that the data do not identify are `NA`, never `0`.
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
#' Efron, B. and Hinkley, D. V. (1978).
#' Assessing the accuracy of the maximum likelihood estimator: observed versus
#' expected Fisher information.
#' *Biometrika*, **65**(3), 457--483.
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
    information = c("observed", "expected"),
    hessian_rel_step = 1e-5,
    trace = FALSE,
    subset = NULL,
    na.action = stats::na.omit,
    model = TRUE,
    x = FALSE,
    y = TRUE) {

  link_spec <- .normalize_simplex_link(link)
  information <- match.arg(information)
  # `subset` is non-standard-evaluated, like lm()/glm(): resolved inside `data`
  # first, then in the caller. See .simplex_eval_subset().
  subset_idx <- .simplex_eval_subset(substitute(subset), data, parent.frame())
  design <- .build_simplex_matrices(
    formula = formula,
    data = data,
    subset = subset_idx,
    na.action = na.action
  )

  response <- design$y
  X <- design$X
  Z <- design$Z

  # Drop aliased columns before fitting and put NA back afterwards, as lm()
  # does. Fitting on the full rank-deficient design would return an arbitrary
  # split of one identified effect across the collinear group.
  alias_x <- .simplex_aliased(X)
  alias_z <- .simplex_aliased(Z)
  if (any(alias_x) || any(alias_z)) {
    warning("Design is rank deficient: ",
            paste(c(colnames(X)[alias_x], colnames(Z)[alias_z]),
                  collapse = ", "),
            " ", if (sum(alias_x, alias_z) == 1L) "is" else "are",
            " a linear combination of the other columns and cannot be ",
            "estimated. Coefficient(s) reported as NA.", call. = FALSE)
  }
  X_full <- X
  Z_full <- Z
  X <- X[, !alias_x, drop = FALSE]
  Z <- Z[, !alias_z, drop = FALSE]
  if (ncol(X) == 0L || ncol(Z) == 0L) {
    stop("Every column of the ", if (ncol(X) == 0L) "mean" else "dispersion",
         " design is aliased; the model has no estimable parameters.",
         call. = FALSE)
  }
  p <- ncol(X)
  q <- ncol(Z)
  # numeric(0) is the backend's "no offset" sentinel.
  off_mu <- if (is.null(design$offset_mu)) numeric(0) else design$offset_mu
  off_phi <- if (is.null(design$offset_phi)) numeric(0) else design$offset_phi

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
    trace = isTRUE(trace),
    off_mu_ = off_mu,
    off_phi_ = off_phi
  )

  theta <- as.numeric(opt$par)
  # Dispersion coefficients carry a "(phi)_" prefix in the FULL parameter
  # vector, as betareg does. Distinguishing the two submodels by position alone
  # is not enough: with an intercept in both, names(theta) repeated
  # "(Intercept)", so vcov(fit)["(Intercept)", "(Intercept)"] silently returned
  # the MEAN intercept's variance whatever the user meant, and
  # confint(fit, parm = "(Intercept)") returned two indistinguishable rows.
  # The per-submodel tables in summary()/coef(model=) keep the bare names.
  names(theta) <- .simplex_par_names(colnames(X), colnames(Z))

  # `theta` currently spans the ESTIMABLE parameters only. Keep that vector for
  # everything numerical (prediction, the information matrix); the NA-padded
  # full-design version is built after inference, for reporting.
  theta_est <- theta
  aliased <- c(alias_x, alias_z)
  names(aliased) <- .simplex_par_names(colnames(X_full), colnames(Z_full))
  p_full <- ncol(X_full)
  q_full <- ncol(Z_full)

  pred <- simplex_predict_cpp(theta_est, X, Z, mean_link = link_spec$id,
                              off_mu_ = off_mu, off_phi_ = off_phi)
  logLik_value <- -as.numeric(opt$value)
  k <- length(theta_est)      # estimable parameters only (aliased ones dropped)
  n <- length(response)

  converged <- as.integer(opt$convergence) == 0L
  if (!converged) {
    warning("fastsimplexreg() did not converge (code ", opt$convergence, ": ",
            opt$message, "). Estimates and standard errors are unreliable.",
            call. = FALSE)
  }

  vc <- NULL
  se <- stats::setNames(rep(NA_real_, k), names(theta_est))
  hessian <- NULL
  vcov_rank <- NA_integer_
  vcov_pseudo <- NA
  vcov_eigenvalues <- NULL
  vcov_condition <- NA_real_
  # Standard errors are only computed at a converged (stationary) fit; at a
  # non-converged point the Hessian is meaningless, so leave them NA.
  if (isTRUE(inference) && converged) {
    hessian <- if (information == "expected") {
      # Exact, block diagonal and positive definite by construction; see
      # .simplex_expected_info(). No finite differencing, hence no noise floor.
      .simplex_expected_info(X, Z, as.numeric(pred$eta_mu),
                             as.numeric(pred$eta_phi), link_spec$name)
    } else {
      simplex_hessian_fd_cpp(
        theta = theta_est,
        y = response,
        X = X,
        Z = Z,
        mean_link = link_spec$id,
        rel_step = as.numeric(hessian_rel_step),
        n_threads = as.integer(n_threads),
        off_mu_ = off_mu,
        off_phi_ = off_phi
      )
    }
    est_names <- names(theta)[!aliased]
    dimnames(hessian) <- list(est_names, est_names)

    # Fail-safe inversion: a rank-deficient or indefinite Hessian yields NA
    # standard errors and a warning, never a confident zero. See R/inference.R.
    inf <- .simplex_vcov(hessian, est_names, what = "fastsimplexreg()")
    vc <- inf$vcov
    se <- inf$se
    vcov_rank <- inf$rank
    vcov_pseudo <- inf$pseudo
    vcov_eigenvalues <- inf$eigenvalues
    vcov_condition <- inf$condition
  }

  # Re-expand to the FULL design: aliased coefficients are reported as NA, like
  # lm(), never as an arbitrary share of an identified effect. vcov and the
  # information matrix stay over the ESTIMABLE parameters only -- there is no
  # curvature in an aliased direction to report.
  if (any(aliased)) {
    theta <- stats::setNames(rep(NA_real_, length(aliased)), names(aliased))
    theta[!aliased] <- theta_est
    se_full <- stats::setNames(rep(NA_real_, length(aliased)), names(aliased))
    se_full[!aliased] <- se
    se <- se_full
  }

  # Saturation of the mean link is reported, not applied silently.
  .warn_saturated(opt$n_saturated, n, what = "fastsimplexreg()")

  # Observation labels, so fitted()/residuals()/predict() can be aligned back to
  # the source rows without assuming row order; and the na.action object, so
  # na.exclude() can actually do what it promises (see .simplex_pad()).
  obs_names <- rownames(design$model)
  na_act <- attr(design$model, "na.action")

  coefficients <- list(
    mean = stats::setNames(theta[seq_len(p_full)], colnames(X_full)),
    dispersion = stats::setNames(theta[p_full + seq_len(q_full)], colnames(Z_full))
  )

  out <- list(
    call = match.call(),
    formula = design$formula,
    link = list(mean = link_spec$name, dispersion = "log"),
    coefficients = coefficients,
    par = theta,
    aliased = aliased,
    standard_errors = stats::setNames(se, names(theta)),
    vcov = vc,
    vcov_rank = vcov_rank,
    vcov_pseudo = vcov_pseudo,
    vcov_eigenvalues = vcov_eigenvalues,
    vcov_condition = vcov_condition,
    n_saturated = as.integer(opt$n_saturated),
    hessian = hessian,
    na.action = na_act,
    fitted.values = stats::setNames(as.numeric(pred$mu), obs_names),
    dispersion.values = stats::setNames(as.numeric(pred$phi), obs_names),
    linear.predictors = list(
      mean = stats::setNames(as.numeric(pred$eta_mu), obs_names),
      dispersion = stats::setNames(as.numeric(pred$eta_phi), obs_names)
    ),
    residuals = stats::setNames(response - as.numeric(pred$mu), obs_names),
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
    offset = list(mean = design$offset_mu, dispersion = design$offset_phi),
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
  if (isTRUE(x)) out$x <- list(mean = X_full, dispersion = Z_full)
  if (isTRUE(y)) out$y <- response

  structure(out, class = "simplex_fast")
}
