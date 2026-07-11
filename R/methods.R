#' @title Extractor Methods for Simplex Regression Fits
#'
#' @description
#' Standard extractor methods for objects of class `"simplex_fast"` produced by
#' [fastsimplexreg()]. They mirror the conventions of the corresponding generics
#' for linear and generalised linear models.
#'
#' \describe{
#'   \item{`coef`}{Returns the estimated coefficients. `model = "all"` returns
#'     the full vector, `"mean"` the mean submodel and `"dispersion"` the
#'     dispersion submodel.}
#'   \item{`vcov`}{Returns the variance-covariance matrix of the estimates.
#'     Errors if the model was fitted with `inference = FALSE`.}
#'   \item{`logLik`}{Returns the maximised log-likelihood, with attributes `df`
#'     (number of estimated parameters) and `nobs`, and class `"logLik"`.}
#'   \item{`nobs`}{Returns the number of observations used in the fit.}
#'   \item{`fitted`}{Returns the fitted means (`model = "mean"`) or fitted
#'     dispersions (`model = "dispersion"`).}
#'   \item{`residuals`}{Returns residuals. `type = "response"` gives
#'     \eqn{y - \hat\mu}; `type = "pearson"` gives
#'     \eqn{(y - \hat\mu) / \sqrt{\hat\phi\, V(\hat\mu)}} with the simplex unit
#'     variance function \eqn{V(\mu) = \{\mu(1-\mu)\}^3}, i.e. the first-order
#'     dispersion-model approximation \eqn{\mathrm{Var}(Y) \approx \phi\,V(\mu)};
#'     `type = "deviance"` gives the signed deviance residual
#'     \eqn{\mathrm{sign}(y-\hat\mu)\sqrt{d(y;\hat\mu)/\hat\phi}} based on the
#'     simplex unit deviance
#'     \eqn{d(y;\mu) = (y-\mu)^2 / \{y(1-y)\mu^2(1-\mu)^2\}}.}
#'   \item{`confint`}{Returns Wald confidence intervals
#'     \eqn{\hat\theta \pm z_{1-\alpha/2}\,\mathrm{SE}(\hat\theta)}, using the
#'     covariance matrix from `vcov`.}
#'   \item{`deviance`}{Returns the scaled deviance
#'     \eqn{\sum_i d(y_i;\hat\mu_i)/\hat\phi_i}, equal to the sum of squared
#'     deviance residuals.}
#'   \item{`model.matrix`}{Returns the mean (`model = "mean"`) or dispersion
#'     (`model = "dispersion"`) design matrix. Requires the fit to have stored
#'     the design (`x = TRUE`) or the model frame (`model = TRUE`).}
#'   \item{`terms`}{Returns the `"terms"` object for the mean or dispersion
#'     submodel.}
#'   \item{`formula`}{Returns the (multi-part) model formula.}
#'   \item{`model.frame`}{Returns the stored model frame (requires
#'     `model = TRUE` at fitting time).}
#'   \item{`update`}{Refits the model with a modified formula and/or arguments.
#'     The multi-part (`|`) structure is updated correctly through the
#'     \pkg{Formula} package.}
#' }
#'
#' @param object,x,formula An object of class `"simplex_fast"` (the argument is
#'   named `object`, `x` or `formula` to match the corresponding generic).
#' @param formula. For `update`, a change to the model formula, following the
#'   conventions of [stats::update.formula()]; the two-part mean/dispersion
#'   structure is handled via the \pkg{Formula} package.
#' @param evaluate For `update`, logical; if `TRUE` (default) the updated call
#'   is evaluated and the refitted model returned, otherwise the updated call is
#'   returned unevaluated.
#' @param model For `coef`, one of `"all"`, `"mean"` or `"dispersion"`; for
#'   `fitted`, `model.matrix` and `terms`, one of `"mean"` or `"dispersion"`.
#' @param type For `residuals`, one of `"response"`, `"pearson"` or
#'   `"deviance"`.
#' @param parm For `confint`, a specification of which parameters to report,
#'   either a vector of numeric indices or of names. Defaults to all.
#' @param level For `confint`, the confidence level.
#' @param ... Additional arguments, currently ignored.
#'
#' @return `coef`, `fitted`, `residuals` and `nobs` return numeric vectors;
#'   `vcov`, `confint` and `model.matrix` return matrices; `logLik` returns a
#'   `"logLik"` object; `deviance` returns a single number; `terms`, `formula`
#'   and `model.frame` return the corresponding model-description objects.
#'
#' @seealso [fastsimplexreg()], [summary.simplex_fast()], [predict.simplex_fast()],
#'   [plot.simplex_fast()], [simulate.simplex_fast()]
#'
#' @examples
#' set.seed(1)
#' n <- 300
#' dat <- data.frame(x1 = rnorm(n), z1 = rnorm(n))
#' mu <- simplex_linkinv(0.2 + 0.7 * dat$x1, link = "logit")
#' dat$y <- rsimplex(n, mu, exp(-0.5 + 0.4 * dat$z1))
#' fit <- fastsimplexreg(y ~ x1 | z1, data = dat, n_threads = 1L, x = TRUE)
#' coef(fit)
#' coef(fit, model = "mean")
#' vcov(fit)
#' logLik(fit)
#' AIC(fit)
#' nobs(fit)
#' deviance(fit)
#' head(fitted(fit))
#' head(residuals(fit, type = "deviance"))
#' confint(fit)
#' head(model.matrix(fit, model = "mean"))
#' formula(fit)
#'
#' @name simplex_fast-methods
#' @rdname simplex_fast-methods
#' @export
coef.simplex_fast <- function(object, model = c("all", "mean", "dispersion"), ...) {
  model <- match.arg(model)
  switch(
    model,
    all = object$par,
    mean = object$coefficients$mean,
    dispersion = object$coefficients$dispersion
  )
}


#' @rdname simplex_fast-methods
#' @export
vcov.simplex_fast <- function(object, ...) {
  if (is.null(object$vcov)) {
    stop("Covariance matrix was not computed. Refit with inference = TRUE.", call. = FALSE)
  }
  object$vcov
}


#' @rdname simplex_fast-methods
#' @export
logLik.simplex_fast <- function(object, ...) {
  out <- object$logLik
  attr(out, "df") <- length(object$par)
  attr(out, "nobs") <- object$nobs
  class(out) <- "logLik"
  out
}


#' @rdname simplex_fast-methods
#' @export
nobs.simplex_fast <- function(object, ...) object$nobs


#' @rdname simplex_fast-methods
#' @export
fitted.simplex_fast <- function(object, model = c("mean", "dispersion"), ...) {
  model <- match.arg(model)
  if (model == "mean") object$fitted.values else object$dispersion.values
}


#' @rdname simplex_fast-methods
#' @export
residuals.simplex_fast <- function(object, type = c("response", "pearson", "deviance"), ...) {
  type <- match.arg(type)
  mu <- object$fitted.values
  y <- mu + object$residuals  # reconstruct the response as mu + (y - mu)
  phi <- object$dispersion.values
  switch(
    type,
    response = y - mu,
    # Pearson residuals use the simplex unit variance function
    # V(mu) = {mu (1 - mu)}^3 scaled by the dispersion phi, i.e. the first-order
    # dispersion-model approximation Var(Y) ~ phi * V(mu).
    pearson = (y - mu) / sqrt(phi * (mu * (1 - mu))^3),
    # Signed deviance residuals from the simplex unit deviance
    # d(y; mu) = (y - mu)^2 / {y (1 - y) mu^2 (1 - mu)^2}.
    deviance = {
      d <- (y - mu)^2 / (y * (1 - y) * mu^2 * (1 - mu)^2)
      sign(y - mu) * sqrt(d / phi)
    }
  )
}


#' @rdname simplex_fast-methods
#' @export
deviance.simplex_fast <- function(object, ...) {
  mu <- object$fitted.values
  y <- mu + object$residuals
  phi <- object$dispersion.values
  d <- (y - mu)^2 / (y * (1 - y) * mu^2 * (1 - mu)^2)
  sum(d / phi)
}


#' @rdname simplex_fast-methods
#' @export
model.matrix.simplex_fast <- function(object, model = c("mean", "dispersion"), ...) {
  model <- match.arg(model)
  key <- if (model == "mean") "mean" else "dispersion"

  if (!is.null(object$x)) {
    return(object$x[[key]])
  }
  if (is.null(object$model)) {
    stop(
      "Design matrix not stored. Refit with x = TRUE or model = TRUE.",
      call. = FALSE
    )
  }

  if (model == "mean") {
    mm <- stats::model.matrix(
      object$terms$mean,
      data = object$model,
      contrasts.arg = object$design$contrasts_mean
    )
  } else if (isTRUE(object$design$has_dispersion_formula)) {
    mm <- stats::model.matrix(
      object$terms$dispersion,
      data = object$model,
      contrasts.arg = object$design$contrasts_dispersion
    )
  } else {
    mm <- matrix(
      1.0,
      nrow = object$nobs,
      ncol = 1L,
      dimnames = list(rownames(object$model), "(Intercept)")
    )
  }
  mm
}


#' @rdname simplex_fast-methods
#' @export
terms.simplex_fast <- function(x, model = c("mean", "dispersion"), ...) {
  model <- match.arg(model)
  x$terms[[model]]
}


#' @rdname simplex_fast-methods
#' @export
formula.simplex_fast <- function(x, ...) {
  stats::formula(x$formula)
}


#' @rdname simplex_fast-methods
#' @export
model.frame.simplex_fast <- function(formula, ...) {
  if (is.null(formula$model)) {
    stop("Model frame not stored. Refit with model = TRUE.", call. = FALSE)
  }
  formula$model
}


#' @rdname simplex_fast-methods
#' @export
update.simplex_fast <- function(object, formula., ..., evaluate = TRUE) {
  call <- object$call
  if (is.null(call)) {
    stop("The fitted object does not contain a call to update.", call. = FALSE)
  }
  extras <- match.call(expand.dots = FALSE)$`...`

  if (!missing(formula.)) {
    # Update through the stored Formula object so that the two-part
    # mean | dispersion structure is preserved; stats::update.formula would
    # mishandle the '|' operator.
    call$formula <- stats::update(object$formula, formula.)
  }

  if (length(extras)) {
    existing <- !is.na(match(names(extras), names(call)))
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
    if (any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }

  if (evaluate) eval(call, parent.frame()) else call
}


#' @rdname simplex_fast-methods
#' @export
confint.simplex_fast <- function(object, parm, level = 0.95, ...) {
  if (is.null(object$vcov)) {
    stop("Covariance matrix was not computed. Refit with inference = TRUE.", call. = FALSE)
  }
  est <- object$par
  se <- object$standard_errors
  pnames <- names(est)

  # Select parameters by position so that duplicated coefficient names (e.g. a
  # "(Intercept)" in both the mean and dispersion submodels) are never confused.
  if (missing(parm) || is.null(parm)) {
    idx <- seq_along(est)
  } else if (is.numeric(parm)) {
    idx <- as.integer(parm)
  } else {
    idx <- which(pnames %in% parm)
  }
  idx <- idx[idx >= 1L & idx <= length(est)]
  if (!length(idx)) {
    stop("No valid parameters selected in 'parm'.", call. = FALSE)
  }

  a <- (1 - level) / 2
  z <- stats::qnorm(1 - a)
  ci <- cbind(est[idx] - z * se[idx], est[idx] + z * se[idx])
  colnames(ci) <- paste0(format(100 * c(a, 1 - a), trim = TRUE, digits = 3), " %")
  rownames(ci) <- pnames[idx]
  ci
}


#' @title Predictions from a Simplex Regression Fit
#'
#' @description
#' Computes predictions from a fitted `"simplex_fast"` model, either on the
#' data used for fitting or on new data. When `newdata` is supplied, the stored
#' `terms`, `xlevels` and `contrasts` are reused so that the design matrices are
#' built consistently with the fit.
#'
#' @param object An object of class `"simplex_fast"`.
#' @param newdata Optional `data.frame` of new observations. When `NULL`, the
#'   fitted values are returned.
#' @param type Type of prediction: `"response"` or `"mean"` (fitted mean
#'   \eqn{\mu}), `"dispersion"` (fitted \eqn{\phi}), `"link"` (a list with the
#'   linear predictors `mean` and `dispersion`) or `"both"` (a `data.frame`
#'   with columns `mu` and `phi`).
#' @param ... Additional arguments, currently ignored.
#'
#' @return A numeric vector, a list or a `data.frame`, depending on `type`.
#'
#' @seealso [fastsimplexreg()]
#'
#' @examples
#' set.seed(2)
#' n <- 300
#' dat <- data.frame(x1 = rnorm(n), z1 = rnorm(n))
#' mu <- simplex_linkinv(0.1 + 0.6 * dat$x1, link = "logit")
#' dat$y <- rsimplex(n, mu, exp(-0.5 + 0.3 * dat$z1))
#' fit <- fastsimplexreg(y ~ x1 | z1, data = dat, n_threads = 1L)
#' head(predict(fit, type = "response"))
#' head(predict(fit, newdata = dat[1:5, ], type = "both"))
#'
#' @export
predict.simplex_fast <- function(
    object,
    newdata = NULL,
    type = c("response", "mean", "dispersion", "link", "both"),
    ...) {
  type <- match.arg(type)

  if (is.null(newdata)) {
    mu <- object$fitted.values
    phi <- object$dispersion.values
    eta_mu <- object$linear.predictors$mean
    eta_phi <- object$linear.predictors$dispersion
  } else {
    X <- stats::model.matrix(
      stats::delete.response(object$design$terms_mean),
      data = newdata,
      contrasts.arg = object$design$contrasts_mean,
      xlev = object$design$xlevels_mean
    )

    if (isTRUE(object$design$has_dispersion_formula)) {
      Z <- stats::model.matrix(
        stats::delete.response(object$design$terms_dispersion),
        data = newdata,
        contrasts.arg = object$design$contrasts_dispersion,
        xlev = object$design$xlevels_dispersion
      )
    } else {
      Z <- matrix(
        1.0,
        nrow = nrow(X),
        ncol = 1L,
        dimnames = list(rownames(X), "(Intercept)")
      )
    }

    storage.mode(X) <- "double"
    storage.mode(Z) <- "double"
    pred <- simplex_predict_cpp(
      object$par,
      X,
      Z,
      mean_link = unname(.simplex_links[[object$link$mean]])
    )
    mu <- as.numeric(pred$mu)
    phi <- as.numeric(pred$phi)
    eta_mu <- as.numeric(pred$eta_mu)
    eta_phi <- as.numeric(pred$eta_phi)
  }

  switch(
    type,
    response = mu,
    mean = mu,
    dispersion = phi,
    link = list(mean = eta_mu, dispersion = eta_phi),
    both = data.frame(mu = mu, phi = phi)
  )
}


#' @title Print a Simplex Regression Fit
#'
#' @description
#' Compactly prints a fitted `"simplex_fast"` object: the formula, links,
#' number of observations, fit statistics and the mean and dispersion
#' coefficients.
#'
#' @param x An object of class `"simplex_fast"`.
#' @param digits Integer; the number of significant digits to display.
#' @param ... Additional arguments, currently ignored.
#'
#' @return The object `x`, invisibly.
#'
#' @seealso [fastsimplexreg()], [summary.simplex_fast()]
#'
#' @examples
#' set.seed(3)
#' dat <- data.frame(x1 = rnorm(200))
#' dat$y <- rsimplex(200, simplex_linkinv(0.3 + 0.5 * dat$x1, "logit"), 1)
#' fit <- fastsimplexreg(y ~ x1, data = dat, n_threads = 1L)
#' print(fit)
#'
#' @export
print.simplex_fast <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nFast simplex regression with variable dispersion\n")
  cat("Formula: ")
  print(x$formula)
  cat("Mean link:", x$link$mean, "| Dispersion link:", x$link$dispersion, "\n")
  cat("Observations:", x$nobs, "\n")
  cat("Log-likelihood:", formatC(x$logLik, digits = digits, format = "fg"), "\n")
  cat("AIC:", formatC(x$AIC, digits = digits, format = "fg"), "\n")
  cat("BIC:", formatC(x$BIC, digits = digits, format = "fg"), "\n")
  cat("Convergence:", x$convergence, "-", x$message, "\n\n")

  cat("Mean coefficients [", x$link$mean, " link]:\n", sep = "")
  print(round(x$coefficients$mean, digits))
  cat("\nDispersion coefficients [log link]:\n")
  print(round(x$coefficients$dispersion, digits))
  invisible(x)
}


#' @title Summarise a Simplex Regression Fit
#'
#' @description
#' Produces a summary of a fitted `"simplex_fast"` object, including coefficient
#' tables with standard errors, Wald z-statistics and p-values for both the mean
#' and dispersion submodels.
#'
#' @param object An object of class `"simplex_fast"`.
#' @param x An object of class `"summary.simplex_fast"`.
#' @param digits Integer; the number of significant digits to display.
#' @param ... Additional arguments, currently ignored.
#'
#' @return An object of class `"summary.simplex_fast"`, a list whose main
#'   component `coefficients` is itself a list with the `mean` and `dispersion`
#'   coefficient tables (each with columns `Estimate`, `Std. Error`, `z value`
#'   and `Pr(>|z|)`), together with the Pearson residuals, the links, fit
#'   statistics (log-likelihood, AIC, BIC, deviance) and optimiser diagnostics.
#'   The `print` method returns its argument invisibly.
#'
#' @seealso [fastsimplexreg()], [print.simplex_fast()]
#'
#' @examples
#' set.seed(4)
#' dat <- data.frame(x1 = rnorm(300), z1 = rnorm(300))
#' mu <- simplex_linkinv(0.2 + 0.6 * dat$x1, "logit")
#' dat$y <- rsimplex(300, mu, exp(-0.4 + 0.3 * dat$z1))
#' fit <- fastsimplexreg(y ~ x1 | z1, data = dat, n_threads = 1L)
#' summary(fit)
#'
#' @export
summary.simplex_fast <- function(object, ...) {
  p <- length(object$coefficients$mean)
  q <- length(object$coefficients$dispersion)
  est <- object$par
  se <- object$standard_errors
  z <- est / se
  pval <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)

  full <- cbind(
    Estimate = est,
    `Std. Error` = se,
    `z value` = z,
    `Pr(>|z|)` = pval
  )

  mean_tab <- full[seq_len(p), , drop = FALSE]
  disp_tab <- full[p + seq_len(q), , drop = FALSE]
  rownames(mean_tab) <- names(object$coefficients$mean)
  rownames(disp_tab) <- names(object$coefficients$dispersion)

  structure(
    list(
      call = object$call,
      formula = object$formula,
      link = object$link,
      # A list with separate mean and dispersion coefficient tables, matching
      # the layout used by other simplex/beta regression packages.
      coefficients = list(mean = mean_tab, dispersion = disp_tab),
      pearson.residuals = stats::residuals(object, type = "pearson"),
      logLik = object$logLik,
      AIC = object$AIC,
      BIC = object$BIC,
      deviance = stats::deviance(object),
      nobs = object$nobs,
      convergence = object$convergence,
      message = object$message,
      iterations = object$iterations,
      function_evaluations = object$function_evaluations,
      gradient_evaluations = object$gradient_evaluations
    ),
    class = "summary.simplex_fast"
  )
}


#' @rdname summary.simplex_fast
#' @export
print.summary.simplex_fast <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n")
  print(x$call)

  cat("\nPearson residuals:\n")
  res_q <- stats::quantile(x$pearson.residuals, c(0, 0.25, 0.5, 0.75, 1), names = FALSE)
  names(res_q) <- c("Min", "1Q", "Median", "3Q", "Max")
  print(round(res_q, digits + 1L))

  cat("\nCoefficients (mean model with ", x$link$mean, " link):\n", sep = "")
  printCoefmat(x$coefficients$mean, digits = digits, signif.legend = FALSE)

  cat("\nCoefficients (dispersion model with ", x$link$dispersion, " link):\n", sep = "")
  printCoefmat(x$coefficients$dispersion, digits = digits)

  cat("\nLog-likelihood:", formatC(x$logLik, digits = digits, format = "fg"))
  cat(" | AIC:", formatC(x$AIC, digits = digits, format = "fg"))
  cat(" | BIC:", formatC(x$BIC, digits = digits, format = "fg"), "\n")
  cat("Deviance:", formatC(x$deviance, digits = digits, format = "fg"))
  cat(" | Observations:", x$nobs, "| Iterations:", x$iterations, "\n")
  cat("Convergence:", x$convergence, "-", x$message, "\n")
  invisible(x)
}


#' Simulate Responses from a Simplex Regression Fit
#'
#' Simulates new response vectors from a fitted `"simplex_fast"` model by drawing
#' from the simplex distribution at the fitted means \eqn{\hat\mu_i} and
#' dispersions \eqn{\hat\phi_i}, using [rsimplex()].
#'
#' @param object An object of class `"simplex_fast"`.
#' @param nsim Number of response vectors to simulate.
#' @param seed Optional seed for the random number generator. Handled following
#'   the convention of [stats::simulate()]: when supplied, the current RNG state
#'   is restored on exit and the seed is recorded in the `"seed"` attribute of
#'   the result.
#' @param ... Additional arguments, currently ignored.
#'
#' @return A `data.frame` with `nsim` columns (`sim_1`, `sim_2`, ...), each a
#'   simulated response of length `nobs(object)`.
#'
#' @seealso [rsimplex()], [fastsimplexreg()]
#'
#' @examples
#' set.seed(5)
#' dat <- data.frame(x1 = rnorm(200))
#' dat$y <- rsimplex(200, simplex_linkinv(0.3 + 0.5 * dat$x1, "logit"), 1)
#' fit <- fastsimplexreg(y ~ x1, data = dat, n_threads = 1L)
#' sims <- simulate(fit, nsim = 3, seed = 42)
#' str(sims)
#'
#' @export
simulate.simplex_fast <- function(object, nsim = 1, seed = NULL, ...) {
  nsim <- as.integer(nsim)
  if (length(nsim) != 1L || is.na(nsim) || nsim < 1L) {
    stop("'nsim' must be a single positive integer.", call. = FALSE)
  }

  # RNG-state bookkeeping, mirroring stats:::simulate.lm.
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    stats::runif(1)
  }
  if (is.null(seed)) {
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  } else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }

  n <- object$nobs
  mu <- object$fitted.values
  phi <- object$dispersion.values

  val <- as.data.frame(
    replicate(nsim, rsimplex(n, mu = mu, phi = phi), simplify = "matrix")
  )
  names(val) <- paste0("sim_", seq_len(nsim))
  if (!is.null(object$y)) {
    row.names(val) <- names(object$y)
  }
  attr(val, "seed") <- RNGstate
  val
}
