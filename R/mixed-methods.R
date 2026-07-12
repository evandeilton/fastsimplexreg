# mixed-methods.R
# S3 methods for objects of class "simplex_fast_mixed" produced by
# fastsimplexregmixed().

# Re-export the nlme generics so that ranef() and VarCorr() are usable directly
# from fastsimplexreg without attaching nlme, while still dispatching to nlme's
# own methods for lme/nlme objects.

#' @importFrom nlme ranef
#' @export
nlme::ranef

#' @importFrom nlme VarCorr
#' @export
nlme::VarCorr

#' Number of Groups in a Mixed-Model Fit
#'
#' Generic and method returning the number of groups (clusters) of the single
#' grouping factor in a fitted `"simplex_fast_mixed"` model.
#'
#' @param object A fitted model object.
#' @param ... Additional arguments, currently ignored.
#' @return An integer, the number of groups.
#' @export
ngrps <- function(object, ...) UseMethod("ngrps")

#' @rdname ngrps
#' @export
ngrps.simplex_fast_mixed <- function(object, ...) object$ngrps


#' Extractor Methods for Simplex Mixed-Model Fits
#'
#' @description
#' Standard extractor methods for objects of class `"simplex_fast_mixed"`
#' produced by [fastsimplexregmixed()].
#'
#' \describe{
#'   \item{`coef`}{Fixed-effect coefficients (`model = "all"`, `"mean"` or
#'     `"dispersion"`). Random effects are obtained with `ranef()`.}
#'   \item{`vcov`}{Covariance matrix of the estimated parameters
#'     `c(beta, gamma, omega)`.}
#'   \item{`logLik`}{Maximised marginal log-likelihood, with `df = p + r +
#'     q(q+1)/2` and `nobs`.}
#'   \item{`nobs`}{Number of observations.}
#'   \item{`ngrps`}{Number of groups.}
#'   \item{`fitted`}{Fitted means (conditional on the empirical-Bayes random
#'     effects) or fitted dispersions.}
#'   \item{`residuals`}{Response, Pearson or deviance residuals, conditional on
#'     the empirical-Bayes random effects.}
#'   \item{`ranef`}{Empirical-Bayes random-effect modes (a groups-by-`q`
#'     matrix); with `postVar = TRUE`, the posterior covariances are attached as
#'     the `"postVar"` attribute.}
#'   \item{`VarCorr`}{The estimated random-effect covariance matrix
#'     \eqn{\Sigma}, with standard deviations and correlations.}
#' }
#'
#' @param object,x A fitted `"simplex_fast_mixed"` object.
#' @param model For `coef`, one of `"all"`, `"mean"` or `"dispersion"`; for
#'   `fitted`, one of `"mean"` or `"dispersion"`.
#' @param type For `residuals`, one of `"response"`, `"pearson"` or
#'   `"deviance"`.
#' @param postVar For `ranef`, logical; attach posterior covariances.
#' @param sigma For `VarCorr`, an optional scale multiplier (kept for
#'   compatibility with the generic; defaults to 1).
#' @param digits For the `VarCorr` print method, the number of significant
#'   digits to display.
#' @param ... Additional arguments, currently ignored.
#'
#' @return `coef`, `fitted` and `residuals` return numeric vectors; `vcov`
#'   returns a matrix; `ranef` returns a matrix; `VarCorr` returns the
#'   covariance matrix with `stddev`/`correlation` attributes; `logLik` returns
#'   a `"logLik"` object.
#'
#' @seealso [fastsimplexregmixed()]
#'
#' @name simplex_fast_mixed-methods
#' @rdname simplex_fast_mixed-methods
#' @export
coef.simplex_fast_mixed <- function(object, model = c("all", "mean", "dispersion"), ...) {
  model <- match.arg(model)
  switch(
    model,
    all = c(object$coefficients$mean, object$coefficients$dispersion),
    mean = object$coefficients$mean,
    dispersion = object$coefficients$dispersion
  )
}

#' @rdname simplex_fast_mixed-methods
#' @export
vcov.simplex_fast_mixed <- function(object, ...) {
  if (is.null(object$vcov)) {
    stop("Covariance matrix was not computed. Refit with inference = TRUE.", call. = FALSE)
  }
  object$vcov
}

#' @rdname simplex_fast_mixed-methods
#' @export
logLik.simplex_fast_mixed <- function(object, ...) {
  out <- object$logLik
  attr(out, "df") <- object$df
  attr(out, "nobs") <- object$nobs
  class(out) <- "logLik"
  out
}

#' @rdname simplex_fast_mixed-methods
#' @export
nobs.simplex_fast_mixed <- function(object, ...) object$nobs

#' @rdname simplex_fast_mixed-methods
#' @export
fitted.simplex_fast_mixed <- function(object, model = c("mean", "dispersion"), ...) {
  model <- match.arg(model)
  if (model == "mean") object$fitted.values else object$dispersion.values
}

#' @rdname simplex_fast_mixed-methods
#' @export
residuals.simplex_fast_mixed <- function(object, type = c("response", "pearson", "deviance"), ...) {
  type <- match.arg(type)
  mu <- object$fitted.values
  y <- mu + object$residuals
  phi <- object$dispersion.values
  switch(
    type,
    response = y - mu,
    pearson = (y - mu) / sqrt(phi * (mu * (1 - mu))^3),
    deviance = {
      d <- (y - mu)^2 / (y * (1 - y) * mu^2 * (1 - mu)^2)
      sign(y - mu) * sqrt(d / phi)
    }
  )
}

#' @rdname simplex_fast_mixed-methods
#' @importFrom nlme ranef
#' @export
ranef.simplex_fast_mixed <- function(object, postVar = FALSE, ...) {
  out <- object$ranef
  if (isTRUE(postVar)) {
    attr(out, "postVar") <- object$ranef.postvar
  }
  out
}

#' @rdname simplex_fast_mixed-methods
#' @importFrom nlme VarCorr
#' @export
VarCorr.simplex_fast_mixed <- function(x, sigma = 1, ...) {
  Sigma <- x$Sigma
  sd <- sqrt(diag(Sigma))
  corr <- suppressWarnings(stats::cov2cor(Sigma))
  structure(Sigma, stddev = sd, correlation = corr, group = x$group_name,
            class = "VarCorr.simplex_fast_mixed")
}

#' @rdname simplex_fast_mixed-methods
#' @export
print.VarCorr.simplex_fast_mixed <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  sd <- attr(x, "stddev")
  corr <- attr(x, "correlation")
  cat("Random effects covariance (group: ", attr(x, "group"), ")\n", sep = "")
  tab <- cbind(Variance = diag(unclass(x)), `Std.Dev.` = sd)
  print(round(tab, digits))
  if (length(sd) > 1L) {
    cat("\nCorrelations:\n")
    print(round(corr, digits))
  }
  invisible(x)
}


#' Predictions from a Simplex Mixed-Model Fit
#'
#' @param object A fitted `"simplex_fast_mixed"` object.
#' @param newdata Optional new data. When `NULL`, in-sample predictions are
#'   returned.
#' @param type Type of prediction: `"response"`/`"mean"`, `"dispersion"`,
#'   `"link"` or `"both"`.
#' @param re.form Controls the random effects. `NULL` (default) includes the
#'   estimated random effects for groups seen in the fit; `NA` (or `~0`) gives
#'   population-level predictions (random effects set to zero).
#' @param ... Additional arguments, currently ignored.
#'
#' @return A numeric vector, list or `data.frame`, depending on `type`.
#'
#' @seealso [fastsimplexregmixed()]
#' @export
predict.simplex_fast_mixed <- function(object, newdata = NULL,
                                       type = c("response", "mean", "dispersion", "link", "both"),
                                       re.form = NULL, ...) {
  type <- match.arg(type)
  population <- (length(re.form) == 1L && is.na(re.form)) ||
    (inherits(re.form, "formula") && identical(all.vars(re.form), character(0)))

  if (is.null(newdata)) {
    if (population) {
      d <- object$design
      X <- stats::model.matrix(stats::delete.response(d$terms_mean),
                               data = object$model, contrasts.arg = d$contrasts_mean,
                               xlev = d$xlevels_mean)
      W <- stats::model.matrix(stats::delete.response(d$terms_dispersion),
                               data = object$model, contrasts.arg = d$contrasts_dispersion,
                               xlev = d$xlevels_dispersion)
      beta <- object$coefficients$mean
      gamma <- object$coefficients$dispersion
      eta_mu <- as.numeric(X %*% beta)
      eta_phi <- as.numeric(W %*% gamma)
      mu <- simplex_linkinv(eta_mu, object$link$mean)
      phi <- exp(eta_phi)
    } else {
      mu <- object$fitted.values
      phi <- object$dispersion.values
      eta_mu <- object$linear.predictors$mean
      eta_phi <- object$linear.predictors$dispersion
    }
  } else {
    d <- object$design
    X <- stats::model.matrix(stats::delete.response(d$terms_mean), data = newdata,
                             contrasts.arg = d$contrasts_mean, xlev = d$xlevels_mean)
    W <- stats::model.matrix(stats::delete.response(d$terms_dispersion), data = newdata,
                             contrasts.arg = d$contrasts_dispersion, xlev = d$xlevels_dispersion)
    Z <- stats::model.matrix(stats::delete.response(d$terms_random), data = newdata,
                             contrasts.arg = d$contrasts_random, xlev = d$xlevels_random)
    beta <- object$coefficients$mean
    gamma <- object$coefficients$dispersion
    eta_mu <- as.numeric(X %*% beta)
    eta_phi <- as.numeric(W %*% gamma)
    if (!population) {
      # Add random effects for groups present in the fit; zero for unseen groups.
      grp <- as.character(newdata[[d$group_name]])
      known <- match(grp, rownames(object$ranef))
      b <- matrix(0, nrow = nrow(X), ncol = ncol(object$ranef))
      seen <- !is.na(known)
      if (any(seen)) b[seen, ] <- object$ranef[known[seen], , drop = FALSE]
      eta_mu <- eta_mu + rowSums(Z * b)
    }
    mu <- simplex_linkinv(eta_mu, object$link$mean)
    phi <- exp(eta_phi)
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


#' Diagnostic Plots for a Simplex Mixed-Model Fit
#'
#' Diagnostic plots built with \pkg{ggplot2}, sharing the panels of
#' [plot.simplex_fast()] but based on the mixed-model fit (residuals are
#' conditional on the empirical-Bayes random effects).
#'
#' @param x A fitted `"simplex_fast_mixed"` object.
#' @param which Integer subset of `1:4` selecting panels.
#' @param type Residual type used in the panels.
#' @param smooth Logical; add a LOESS smoother.
#' @param ... Additional arguments, currently ignored.
#'
#' @return Invisibly, a `ggplot`/\pkg{patchwork} object or a list of `ggplot`s.
#' @seealso [fastsimplexregmixed()], [plot.simplex_fast()]
#' @export
plot.simplex_fast_mixed <- function(x, which = 1:4,
                                    type = c("deviance", "pearson", "response"),
                                    smooth = TRUE, ...) {
  .simplex_diag_plot(x, which = which, type = match.arg(type), smooth = smooth)
}


#' Print a Simplex Mixed-Model Fit
#'
#' @param x A fitted `"simplex_fast_mixed"` object.
#' @param digits Number of significant digits.
#' @param ... Additional arguments, currently ignored.
#' @return The object `x`, invisibly.
#' @seealso [fastsimplexregmixed()], [summary.simplex_fast_mixed()]
#' @export
print.simplex_fast_mixed <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nFast simplex mixed model with variable dispersion\n")
  cat("Formula: "); print(x$formula)
  cat("Random:  "); print(x$random)
  cat("Mean link:", x$link$mean, "| Dispersion link:", x$link$dispersion, "\n")
  cat("Observations:", x$nobs, "| Groups:", x$ngrps, "| nAGQ:", x$nAGQ, "\n")
  cat("Log-likelihood:", formatC(x$logLik, digits = digits, format = "fg"),
      "| AIC:", formatC(x$AIC, digits = digits, format = "fg"),
      "| BIC:", formatC(x$BIC, digits = digits, format = "fg"), "\n\n")

  cat("Mean coefficients [", x$link$mean, " link]:\n", sep = "")
  print(round(x$coefficients$mean, digits))
  cat("\nDispersion coefficients [log link]:\n")
  print(round(x$coefficients$dispersion, digits))
  cat("\nRandom-effect covariance (group: ", x$group_name, "):\n", sep = "")
  print(round(x$Sigma, digits))
  invisible(x)
}


#' Summarise a Simplex Mixed-Model Fit
#'
#' @param object A fitted `"simplex_fast_mixed"` object.
#' @param x A `"summary.simplex_fast_mixed"` object.
#' @param digits Number of significant digits.
#' @param ... Additional arguments, currently ignored.
#' @return An object of class `"summary.simplex_fast_mixed"`.
#' @seealso [fastsimplexregmixed()]
#' @export
summary.simplex_fast_mixed <- function(object, ...) {
  p <- length(object$coefficients$mean)
  r <- length(object$coefficients$dispersion)
  est <- object$par
  se <- object$standard_errors
  z <- est / se
  pval <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
  full <- cbind(Estimate = est, `Std. Error` = se, `z value` = z, `Pr(>|z|)` = pval)

  mean_tab <- full[seq_len(p), , drop = FALSE]
  disp_tab <- full[p + seq_len(r), , drop = FALSE]
  rownames(mean_tab) <- names(object$coefficients$mean)
  rownames(disp_tab) <- names(object$coefficients$dispersion)

  structure(
    list(
      call = object$call, formula = object$formula, random = object$random,
      link = object$link,
      coefficients = list(mean = mean_tab, dispersion = disp_tab),
      varcorr = VarCorr.simplex_fast_mixed(object),
      pearson.residuals = stats::residuals(object, type = "pearson"),
      logLik = object$logLik, AIC = object$AIC, BIC = object$BIC,
      nobs = object$nobs, ngrps = object$ngrps, nAGQ = object$nAGQ,
      convergence = object$convergence, message = object$message,
      iterations = object$iterations
    ),
    class = "summary.simplex_fast_mixed"
  )
}

#' @rdname summary.simplex_fast_mixed
#' @export
print.summary.simplex_fast_mixed <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n"); print(x$call)

  if (!is.null(x$convergence) && x$convergence != 0L) {
    cat("\n*** MODEL DID NOT CONVERGE (code ", x$convergence, ": ", x$message,
        ") -- results below are UNRELIABLE. ***\n", sep = "")
  }

  cat("\nPearson residuals:\n")
  rq <- stats::quantile(x$pearson.residuals, c(0, 0.25, 0.5, 0.75, 1), names = FALSE)
  names(rq) <- c("Min", "1Q", "Median", "3Q", "Max")
  print(round(rq, digits + 1L))

  cat("\nCoefficients (mean model with ", x$link$mean, " link):\n", sep = "")
  stats::printCoefmat(x$coefficients$mean, digits = digits, signif.legend = FALSE)
  cat("\nCoefficients (dispersion model with ", x$link$dispersion, " link):\n", sep = "")
  stats::printCoefmat(x$coefficients$dispersion, digits = digits)

  cat("\nRandom effects:\n")
  print(x$varcorr, digits = digits)

  cat("\nLog-likelihood:", formatC(x$logLik, digits = digits, format = "fg"),
      "| AIC:", formatC(x$AIC, digits = digits, format = "fg"),
      "| BIC:", formatC(x$BIC, digits = digits, format = "fg"), "\n")
  cat("Observations:", x$nobs, "| Groups:", x$ngrps, "| nAGQ:", x$nAGQ,
      "| Iterations:", x$iterations, "\n")
  cat("Convergence:", x$convergence, "-", x$message, "\n")
  invisible(x)
}
