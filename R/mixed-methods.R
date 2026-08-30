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
#'
#' @examples
#' set.seed(1)
#' J <- 40; nj <- 8; n <- J * nj
#' dat <- data.frame(g = factor(rep(seq_len(J), each = nj)), x1 = rnorm(n))
#' b <- rnorm(J, 0, 0.7)[dat$g]
#' dat$y <- rsimplex(n, simplex_linkinv(0.3 - 0.6 * dat$x1 + b, "logit"), 1)
#' fit <- fastsimplexregmixed(y ~ x1, random = ~ 1 | g, data = dat,
#'                            nAGQ = 7, n_threads = 1)
#' ngrps(fit)
#' @export
ngrps <- function(object, ...) UseMethod("ngrps")

#' @rdname ngrps
#' @export
#' @rawNamespace export(ngrps.simplex_fast_mixed)
# The METHOD is exported, not only registered. `lme4` defines its own,
# independent `ngrps` generic, and its UseMethod() searches lme4's registration
# table -- which cannot contain ours. Attaching lme4 after fastsimplexreg
# therefore masked our generic and `ngrps(fit)` failed with "Cannot extract the
# number of groups from this object". With the method on the search path either
# generic finds it. (`ranef` and `VarCorr` were never affected: lme4 re-exports
# the same nlme generic objects, so there is only one generic in play.)
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
#'   \item{`residuals`}{Randomised quantile (the default), response, Pearson or
#'     deviance residuals, conditional on the empirical-Bayes random effects.}
#'   \item{`ranef`}{Empirical-Bayes random-effect modes (a groups-by-`q`
#'     matrix); with `postVar = TRUE`, the posterior covariances are attached as
#'     the `"postVar"` attribute.}
#'   \item{`VarCorr`}{The estimated random-effect covariance matrix
#'     \eqn{\Sigma}, with standard deviations and correlations.}
#'   \item{`confint`}{Wald confidence intervals over the full parameter vector
#'     `c(beta, gamma, omega)`. Note that the intervals for the variance
#'     components are on the unconstrained log-Cholesky scale, where a Wald
#'     interval is defensible; on the variance scale it would not be, because
#'     the null lies on the boundary.}
#' }
#'
#' @param object,x A fitted `"simplex_fast_mixed"` object.
#' @param model For `coef`, one of `"all"`, `"mean"` or `"dispersion"`; for
#'   `fitted`, one of `"mean"` or `"dispersion"`.
#' @param type For `residuals`, one of `"quantile"` (the default; randomised
#'   quantile residuals in the sense of Dunn and Smyth, 1996, which are exactly
#'   standard normal under a correctly specified model), `"response"`,
#'   `"pearson"` or `"deviance"`.
#' @param postVar For `ranef`, logical; attach posterior covariances.
#' @param parm For `confint`, which parameters to report: numeric positions or
#'   names, over the full vector `c(beta, gamma, omega)`. Defaults to all.
#' @param level For `confint`, the confidence level.
#' @param sigma For `VarCorr`, present only to match the signature of
#'   [nlme::VarCorr()]. A simplex mixed model has no residual scale parameter,
#'   so the argument rescales nothing; supplying anything other than `1` raises
#'   a warning and is ignored.
#' @param digits For the `VarCorr` print method, the number of significant
#'   digits to display.
#' @param ... Additional arguments, currently ignored.
#'
#' @return `coef`, `fitted` and `residuals` return numeric vectors; `vcov`
#'   returns a matrix; `ranef` returns a matrix; `VarCorr` returns the
#'   covariance matrix with `stddev`/`correlation` attributes; `logLik` returns
#'   a `"logLik"` object.
#'
#' @references
#' Dunn, P. K. and Smyth, G. K. (1996). Randomized quantile residuals.
#' *Journal of Computational and Graphical Statistics*, **5**(3), 236--244.
#' \doi{10.1080/10618600.1996.10474708}
#'
#' @seealso [fastsimplexregmixed()]
#'
#' @examples
#' set.seed(1)
#' J <- 40; nj <- 8; n <- J * nj
#' dat <- data.frame(g = factor(rep(seq_len(J), each = nj)), x1 = rnorm(n))
#' b <- rnorm(J, 0, 0.7)[dat$g]
#' dat$y <- rsimplex(n, simplex_linkinv(0.3 - 0.6 * dat$x1 + b, "logit"), 1)
#' fit <- fastsimplexregmixed(y ~ x1, random = ~ 1 | g, data = dat,
#'                            nAGQ = 7, n_threads = 1)
#'
#' coef(fit)
#' VarCorr(fit)
#' head(ranef(fit))
#' confint(fit)
#' ngrps(fit)
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
confint.simplex_fast_mixed <- function(object, parm, level = 0.95, ...) {
  # A method is REQUIRED here, not optional. Without one, dispatch fell through
  # to stats::confint.default, which indexes vcov() by NAME -- and coef() of a
  # mixed fit repeats "(Intercept)" across the mean and dispersion submodels.
  # Both rows then resolved to the first match, so the dispersion intercept was
  # reported with the mean intercept's interval, one that need not even contain
  # its own estimate; and the variance components were dropped entirely, since
  # coef() is shorter than the parameter vector. Selection is by position, over
  # the FULL parameter vector c(beta, gamma, omega).
  if (is.null(object$vcov)) {
    stop("Covariance matrix was not computed. Refit with inference = TRUE.",
         call. = FALSE)
  }
  .simplex_confint(object$par, object$standard_errors, parm, level,
                   missing(parm))
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
  .simplex_pad(object,
               if (model == "mean") object$fitted.values
               else object$dispersion.values)
}

#' @rdname simplex_fast_mixed-methods
#' @export
residuals.simplex_fast_mixed <- function(object, type = c("quantile", "response", "pearson", "deviance"), ...) {
  type <- match.arg(type)
  .simplex_pad(object, .simplex_resid_raw(object, type))
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
  # `sigma` exists only to match the signature of nlme::VarCorr(). A simplex
  # mixed model has no residual scale parameter to factor out of Sigma, so there
  # is nothing for it to rescale. Rather than accepting and discarding a value
  # the user believes is doing something, anything other than 1 is refused.
  if (!isTRUE(all.equal(unname(sigma), 1))) {
    warning("VarCorr(): 'sigma' has no meaning for a simplex mixed model -- ",
            "the covariance Sigma is reported on its own scale and the supplied ",
            "value is ignored.", call. = FALSE)
  }
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
#'
#' @examples
#' set.seed(1)
#' J <- 40; nj <- 8; n <- J * nj
#' dat <- data.frame(g = factor(rep(seq_len(J), each = nj)), x1 = rnorm(n))
#' b <- rnorm(J, 0, 0.7)[dat$g]
#' dat$y <- rsimplex(n, simplex_linkinv(0.3 - 0.6 * dat$x1 + b, "logit"), 1)
#' fit <- fastsimplexregmixed(y ~ x1, random = ~ 1 | g, data = dat,
#'                            nAGQ = 7, n_threads = 1)
#' head(predict(fit))
#' head(predict(fit, re.form = NA))
#' @export
predict.simplex_fast_mixed <- function(object, newdata = NULL,
                                       type = c("response", "mean", "dispersion", "link", "both"),
                                       re.form = NULL, ...) {
  type <- match.arg(type)
  population <- (length(re.form) == 1L && is.na(re.form)) ||
    (inherits(re.form, "formula") && identical(all.vars(re.form), character(0)))
  d <- object$design

  if (is.null(newdata)) {
    if (population) {
      # Population-level in-sample prediction needs the original design back.
      # Use the stored matrices when available, otherwise rebuild them from the
      # stored model frame. With neither, model.matrix() would fall through to
      # the formula's environment and silently build the design from whatever
      # objects happen to be visible in the caller -- returning predictions of
      # the wrong length from the wrong data.
      if (!is.null(object$x)) {
        X <- object$x$mean
        W <- object$x$dispersion
      } else if (!is.null(object$model)) {
        X <- stats::model.matrix(stats::delete.response(d$terms_mean),
                                 data = object$model,
                                 contrasts.arg = d$contrasts_mean,
                                 xlev = d$xlevels_mean)
        W <- stats::model.matrix(stats::delete.response(d$terms_dispersion),
                                 data = object$model,
                                 contrasts.arg = d$contrasts_dispersion,
                                 xlev = d$xlevels_dispersion)
      } else {
        stop("Population-level predictions need the design. Refit with ",
             "model = TRUE or x = TRUE, or supply 'newdata'.", call. = FALSE)
      }
      beta <- object$coefficients$mean
      gamma <- object$coefficients$dispersion
      eta_mu <- as.numeric(X %*% beta)
      eta_phi <- as.numeric(W %*% gamma)
      if (!is.null(object$offset$mean)) eta_mu <- eta_mu + object$offset$mean
      if (!is.null(object$offset$dispersion)) {
        eta_phi <- eta_phi + object$offset$dispersion
      }
      mu <- simplex_linkinv(eta_mu, object$link$mean)
      phi <- exp(eta_phi)
    } else {
      mu <- object$fitted.values
      phi <- object$dispersion.values
      eta_mu <- object$linear.predictors$mean
      eta_phi <- object$linear.predictors$dispersion
    }
  } else {
    terms_list <- list(mean = d$terms_mean,
                       dispersion = d$terms_dispersion,
                       random = d$terms_random)
    contrasts_list <- list(d$contrasts_mean, d$contrasts_dispersion,
                           d$contrasts_random)
    xlev_list <- list(d$xlevels_mean, d$xlevels_dispersion, d$xlevels_random)
    # A conditional prediction needs the grouping factor: without it every row
    # would silently fall back to a zero random effect, i.e. a population-level
    # prediction returned under the label of a conditional one.
    extra <- if (population) character(0) else d$group_name

    des <- .simplex_predict_design(newdata, terms_list, contrasts_list,
                                   xlev_list, extra_vars = extra)
    X <- des$matrices$mean
    W <- des$matrices$dispersion
    Z <- des$matrices$random

    beta <- object$coefficients$mean
    gamma <- object$coefficients$dispersion
    eta_mu <- as.numeric(X %*% beta)
    eta_phi <- as.numeric(W %*% gamma)

    # Offsets belong to the linear predictor and must be rebuilt from newdata.
    nd_keep <- newdata[des$keep, , drop = FALSE]
    off_mu <- .simplex_check_offset(.simplex_offset(d$terms_mean, nd_keep),
                                    nrow(X), "mean")
    off_phi <- .simplex_check_offset(.simplex_offset(d$terms_dispersion, nd_keep),
                                     nrow(W), "dispersion")
    if (!is.null(off_mu)) eta_mu <- eta_mu + off_mu
    if (!is.null(off_phi)) eta_phi <- eta_phi + off_phi

    if (!population) {
      # Add random effects for groups present in the fit; zero for unseen groups.
      grp <- as.character(newdata[[d$group_name]])[des$keep]
      known <- match(grp, rownames(object$ranef))
      b <- matrix(0, nrow = nrow(X), ncol = ncol(object$ranef))
      seen <- !is.na(known)
      if (any(seen)) b[seen, ] <- object$ranef[known[seen], , drop = FALSE]
      if (any(!seen)) {
        warning("predict(): ", sum(!seen), " row(s) of 'newdata' belong to ",
                "group level(s) not seen in the fit; their random effect is ",
                "taken to be zero.", call. = FALSE)
      }
      eta_mu <- eta_mu + rowSums(Z * b)
    }

    mu <- simplex_linkinv(eta_mu, object$link$mean)
    phi <- exp(eta_phi)
    # Keep the result aligned with, and the same length as, newdata.
    mu <- .simplex_expand(mu, des$keep)
    phi <- .simplex_expand(phi, des$keep)
    eta_mu <- .simplex_expand(eta_mu, des$keep)
    eta_phi <- .simplex_expand(eta_phi, des$keep)
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
#' @param type Type of residual used in panels 1-3: `"quantile"` (the default;
#'   randomised quantile residuals), `"deviance"`, `"pearson"` or
#'   `"response"`. See [residuals.simplex_fast_mixed()].
#' @param smooth Logical; add a LOESS smoother.
#' @param ... Additional arguments, currently ignored.
#'
#' @return Invisibly, a `ggplot`/\pkg{patchwork} object or a list of `ggplot`s.
#' @seealso [fastsimplexregmixed()], [plot.simplex_fast()]
#'
#' @examples
#' set.seed(1)
#' J <- 40; nj <- 8; n <- J * nj
#' dat <- data.frame(g = factor(rep(seq_len(J), each = nj)), x1 = rnorm(n))
#' b <- rnorm(J, 0, 0.7)[dat$g]
#' dat$y <- rsimplex(n, simplex_linkinv(0.3 - 0.6 * dat$x1 + b, "logit"), 1)
#' fit <- fastsimplexregmixed(y ~ x1, random = ~ 1 | g, data = dat,
#'                            nAGQ = 7, n_threads = 1)
#' if (requireNamespace("ggplot2", quietly = TRUE)) p <- plot(fit, which = 1:2)
#' @export
plot.simplex_fast_mixed <- function(x, which = 1:4,
                                    type = c("quantile", "deviance", "pearson", "response"),
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
#'
#' @examples
#' set.seed(1)
#' J <- 40; nj <- 8; n <- J * nj
#' dat <- data.frame(g = factor(rep(seq_len(J), each = nj)), x1 = rnorm(n))
#' b <- rnorm(J, 0, 0.7)[dat$g]
#' dat$y <- rsimplex(n, simplex_linkinv(0.3 - 0.6 * dat$x1 + b, "logit"), 1)
#' fit <- fastsimplexregmixed(y ~ x1, random = ~ 1 | g, data = dat,
#'                            nAGQ = 7, n_threads = 1)
#' print(fit)
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
#'
#' @examples
#' set.seed(1)
#' J <- 40; nj <- 8; n <- J * nj
#' dat <- data.frame(g = factor(rep(seq_len(J), each = nj)), x1 = rnorm(n))
#' b <- rnorm(J, 0, 0.7)[dat$g]
#' dat$y <- rsimplex(n, simplex_linkinv(0.3 - 0.6 * dat$x1 + b, "logit"), 1)
#' fit <- fastsimplexregmixed(y ~ x1, random = ~ 1 | g, data = dat,
#'                            nAGQ = 7, n_threads = 1)
#' summary(fit)
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
      # Quantile residuals, not Pearson: Pearson residuals rejected up to
      # 100% of CORRECT models in a Shapiro-Wilk check, because
      # Var(Y) = phi V(mu) only holds to first order (the measured ratio
      # falls to 0.437 at phi = 10).
      quantile.residuals = .simplex_resid_raw(object, "quantile"),
      logLik = object$logLik, AIC = object$AIC, BIC = object$BIC,
      nobs = object$nobs, ngrps = object$ngrps, nAGQ = object$nAGQ,
      vcov_rank = object$vcov_rank,
      vcov_pseudo = object$vcov_pseudo,
      vcov_condition = object$vcov_condition,
      n_saturated = object$n_saturated,
      aliased = object$aliased,
      npar = length(object$par),
      no_inference = is.null(object$vcov),
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

  cat("\nQuantile residuals:\n")
  rq <- stats::quantile(x$quantile.residuals, c(0, 0.25, 0.5, 0.75, 1), names = FALSE)
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
  .simplex_print_diagnostics(x)
  cat("Convergence:", x$convergence, "-", x$message, "\n")
  invisible(x)
}
