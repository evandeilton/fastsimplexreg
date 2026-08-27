# inference.R
# Fail-safe construction of the covariance matrix and standard errors from an
# observed information (negative Hessian) matrix.
#
# Rationale. Inverting the Hessian with solve() and then flooring the diagonal
# with pmax(., 0) is unsafe in exactly the situation that matters: with
# collinear or otherwise non-identified covariates, solve() SUCCEEDS, returns
# large negative variances, and the floor turns them into a standard error of
# exactly 0 -- i.e. z = Inf and p = 0 for a parameter the data cannot identify.
# A missing standard error must never be reported as a confident one, so this
# helper works from the eigen-decomposition of the symmetrised Hessian and:
#
#   * inverts only the strictly positive part of the spectrum (a Moore-Penrose
#     pseudo-inverse restricted to the directions that carry real curvature);
#   * refuses to turn a direction of negative curvature into a variance -- there
#     the fit is not a maximum and 1/lambda is not a variance at all;
#   * returns NA_real_ (never 0) for every parameter that loads on a discarded
#     direction, or whose variance is not strictly positive;
#   * warns, naming the affected parameters.

# Internal: eigen-based fail-safe inverse of an observed information matrix.
# Returns a list with `vcov`, `se`, `rank`, `eigenvalues`, `condition` and
# `pseudo` (TRUE when the plain inverse was not usable and the pseudo-inverse
# was taken).
#
# The rank test is performed on the EQUILIBRATED matrix D^-1 H D^-1 with
# D = diag(sqrt(diag(H))), i.e. on the correlation scale. Without that step the
# test measures scaling rather than rank: with the neglog link the observed
# information of a perfectly ordinary fit spans eigenvalues from 6e2 to 2e15
# (condition 3e12) purely because the parameters live on very different scales,
# while on the correlation scale its condition number is 6e5 -- comfortably full
# rank. Testing the raw matrix declares such a fit unidentified, which is a false
# alarm as damaging as the silent zero it replaced.
.simplex_vcov <- function(hessian, par_names, what = "fastsimplexreg()") {
  d <- nrow(hessian)
  H <- 0.5 * (hessian + t(hessian))

  if (anyNA(H) || !all(is.finite(H))) {
    warning(what, ": the observed information matrix contains non-finite ",
            "entries; no standard errors can be computed.", call. = FALSE)
    vc <- matrix(NA_real_, d, d, dimnames = list(par_names, par_names))
    return(list(vcov = vc, se = stats::setNames(rep(NA_real_, d), par_names),
                rank = NA_integer_, eigenvalues = rep(NA_real_, d),
                condition = NA_real_, pseudo = TRUE))
  }

  # Equilibrate when the diagonal allows it; fall back to the raw matrix when it
  # does not (a non-positive diagonal entry is itself a sign of indefiniteness,
  # which the eigenvalue test below then catches).
  dg <- diag(H)
  scaled <- all(is.finite(dg)) && all(dg > 0)
  sc <- if (scaled) sqrt(dg) else rep(1, d)
  Hs <- H / outer(sc, sc)

  eg <- eigen(Hs, symmetric = TRUE)
  lambda <- eg$values
  V <- eg$vectors

  # Numerical-rank cut-off. The relevant scale is NOT machine epsilon: this
  # matrix is a central-difference approximation of the analytic gradient with
  # relative step h (1e-5 by default), so its entries carry a relative error of
  # order max(h^2, eps/h) ~ 1e-10 -- roughly five orders of magnitude above eps.
  # An eigenvalue below that floor is indistinguishable from zero GIVEN HOW THE
  # MATRIX WAS COMPUTED, and 1/lambda is then pure finite-difference noise
  # dressed up as a variance. sqrt(eps) ~ 1.5e-8 sits comfortably above the
  # noise floor and comfortably below the smallest genuine eigenvalue observed
  # in practice (5e-6 on the correlation scale for the worst-scaled link),
  # which is the separation this test needs.
  #
  # Using eps directly is not merely conservative, it is wrong: it makes the
  # verdict depend on which BLAS computed the Hessian. A design with exactly
  # collinear covariates was correctly flagged on one platform and passed as
  # full rank on another, because the degenerate eigenvalue landed at 1e-13
  # either side of a 1e-15 cut.
  tol <- sqrt(.Machine$double.eps) * max(abs(lambda), 1)

  keep <- lambda > tol
  negative <- lambda < -tol
  rank <- sum(keep)
  condition <- if (rank == d) max(lambda) / min(lambda) else Inf

  inv_lambda <- numeric(d)
  inv_lambda[keep] <- 1 / lambda[keep]
  vc <- V %*% (inv_lambda * t(V))
  vc <- vc / outer(sc, sc)            # undo the equilibration
  vc <- 0.5 * (vc + t(vc))
  dimnames(vc) <- list(par_names, par_names)

  # A parameter is unidentified if it has appreciable loading on a discarded
  # eigen-direction; its variance from the pseudo-inverse is a minimum-norm
  # artefact, not an estimate, so it is reported as NA rather than as a number.
  affected <- rep(FALSE, d)
  if (any(!keep)) {
    loading <- sqrt(rowSums(V[, !keep, drop = FALSE]^2))
    affected <- loading > 1e-6
  }

  variance <- diag(vc)
  se <- rep(NA_real_, d)
  ok <- is.finite(variance) & variance > 0 & !affected
  se[ok] <- sqrt(variance[ok])
  names(se) <- par_names

  if (any(negative)) {
    warning(what, ": the observed information matrix is not positive definite (",
            sum(negative), " negative eigenvalue(s) on the correlation scale, ",
            "smallest = ", format(min(lambda), digits = 3),
            "). The optimiser did not reach a maximum, so the covariance matrix ",
            "is not a valid variance estimate; standard errors for ",
            paste(unique(par_names[affected]), collapse = ", "),
            " are reported as NA.", call. = FALSE)
  } else if (rank < d) {
    warning(what, ": the observed information matrix is rank deficient (rank ",
            rank, " of ", d, "). Parameter(s) ",
            paste(unique(par_names[affected]), collapse = ", "),
            " are not identified by the data -- check for collinear covariates. ",
            "A Moore-Penrose pseudo-inverse was used and their standard errors ",
            "are NA.", call. = FALSE)
  } else if (is.finite(condition) && condition > 1 / sqrt(.Machine$double.eps)) {
    # Full rank, so every standard error is a real number -- but a very large
    # one, and the user should know why rather than be left to wonder.
    warning(what, ": the observed information matrix is ill-conditioned ",
            "(condition number ", format(condition, digits = 3),
            " on the correlation scale). The standard errors are valid but ",
            "large; the data identify some parameters only weakly.",
            call. = FALSE)
  }

  list(vcov = vc, se = se, rank = rank, eigenvalues = lambda,
       condition = condition, pseudo = rank < d)
}


# Internal: derivative dmu/deta of each mean link, on the reporting scale.
.simplex_dmu_deta <- function(eta, link) {
  switch(
    link,
    logit = { m <- stats::plogis(eta); m * (1 - m) },
    probit = stats::dnorm(eta),
    cloglog = exp(eta - exp(eta)),
    neglog = exp(-eta - exp(-eta)),
    stop("Unsupported link.", call. = FALSE)
  )
}


# Internal: EXACT expected (Fisher) information for c(beta, gamma).
#
# The simplex is a proper dispersion model: d(Y; mu)/phi is exactly chi-squared
# with one degree of freedom. Three consequences follow, all exact rather than
# asymptotic:
#
#   * the dispersion score is -1/2 + dev/(2 phi), whose variance is exactly 1/2,
#     so the gamma block is exactly (1/2) Z'Z -- free of the data and of phi;
#   * d l/d eta_mu is proportional to 1/phi and has mean zero, so the cross
#     block E[-d^2 l / d eta_mu d eta_phi] is exactly ZERO: beta and gamma are
#     orthogonal, and the information is block diagonal;
#   * the mean block uses the exact I(mu) = 1/(phi V(mu)) + 3/(mu(1-mu)) with
#     V(mu) = {mu(1-mu)}^3, times (dmu/deta)^2.
#
# Sanity check that fixes all three at once: under constant dispersion this
# gives SE(gamma_0) = sqrt(2/n) exactly, which matches the observed-information
# value to eight digits at n = 4000.
.simplex_expected_info <- function(X, Z, eta_mu, eta_phi, link) {
  mu <- simplex_linkinv(eta_mu, link)
  phi <- exp(eta_phi)
  dmu <- .simplex_dmu_deta(eta_mu, link)
  u <- mu * (1 - mu)
  w_mu <- dmu^2 * (1 / (phi * u^3) + 3 / u)

  p <- ncol(X)
  q <- ncol(Z)
  info <- matrix(0, p + q, p + q)
  info[seq_len(p), seq_len(p)] <- crossprod(X, X * w_mu)
  info[p + seq_len(q), p + seq_len(q)] <- 0.5 * crossprod(Z)
  # The off-diagonal blocks stay exactly zero: beta and gamma are orthogonal.
  0.5 * (info + t(info))
}


# Internal: residuals on the COMPLETE rows only. The residuals() methods wrap
# this in .simplex_pad(); internal consumers (plot, summary) use it directly so
# that they stay aligned with fitted.values, which is never padded.
.simplex_resid_raw <- function(object, type) {
  mu <- object$fitted.values
  y <- .simplex_response(object)
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


# Internal: re-expand a per-observation vector over the rows that na.action
# removed, so that na.exclude() means what stats says it means.
#
# `na.action = na.exclude` was accepted and had no effect: fitted() and
# residuals() came back with the number of COMPLETE rows, not the number of
# rows in the data, so nothing could be aligned back to the source without
# knowing which rows had been dropped. stats::naresid()/napredict() do the
# padding; they are no-ops under na.omit, which keeps the default unchanged.
.simplex_pad <- function(object, values) {
  na_act <- object$na.action
  if (is.null(na_act)) return(values)
  stats::naresid(na_act, values)
}


# Internal: names for the FULL parameter vector. The dispersion block is
# prefixed so that a coefficient appearing in both submodels (typically
# "(Intercept)") is addressable unambiguously by name in coef(), vcov(),
# confint() and every downstream tool that indexes by name. This follows
# betareg, whose coef() reads "(Intercept)", "x1", "(phi)_(Intercept)".
.simplex_par_names <- function(mean_names, disp_names) {
  c(mean_names, paste0("(phi)_", disp_names))
}


# Internal: shared Wald-interval builder for both fit classes.
#
# Selection is BY POSITION throughout. That matters because `parm` may name a
# coefficient that appears in more than one submodel; resolving such a name by
# `%in%` returns every match, and letting stats::confint.default index vcov()
# by name returns the FIRST match for all of them -- which is how a mixed fit
# used to report the mean intercept's interval for the dispersion intercept.
.simplex_confint <- function(est, se, parm, level, missing_parm) {
  pnames <- names(est)
  if (missing_parm || is.null(parm)) {
    idx <- seq_along(est)
  } else if (is.numeric(parm)) {
    idx <- as.integer(parm)
  } else {
    idx <- which(pnames %in% parm)
  }
  idx <- idx[!is.na(idx) & idx >= 1L & idx <= length(est)]
  if (!length(idx)) {
    stop("No valid parameters selected in 'parm'.", call. = FALSE)
  }
  if (length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) {
    stop("'level' must be a single number strictly between 0 and 1.", call. = FALSE)
  }

  a <- (1 - level) / 2
  z <- stats::qnorm(1 - a)
  ci <- cbind(est[idx] - z * se[idx], est[idx] + z * se[idx])
  colnames(ci) <- paste0(format(100 * c(a, 1 - a), trim = TRUE, digits = 3), " %")
  rownames(ci) <- pnames[idx]
  ci
}


# Internal: report saturation of the mean link instead of applying it silently.
# `n` counts the observations whose fitted mean hit the numerical floor of the
# likelihood path, where the score contribution is exactly zero by construction.
.warn_saturated <- function(n, nobs, what = "fastsimplexreg()") {
  n <- as.integer(n)
  if (is.na(n) || n <= 0L) return(invisible(FALSE))
  warning(what, ": the mean saturated at the numerical boundary for ", n,
          " of ", nobs, " observation(s) at the final parameter value. Those ",
          "observations contribute exactly zero to the mean score, so the ",
          "corresponding coefficients are only weakly identified -- check for ",
          "separation or extreme covariate values.", call. = FALSE)
  invisible(TRUE)
}


# Internal: the response used by a fit. Prefer the stored copy (y = TRUE, the
# default); fall back to reconstructing it from the fitted means and the stored
# response residuals when the fit was made with y = FALSE.
.simplex_response <- function(object) {
  if (!is.null(object$y)) return(as.numeric(object$y))
  object$fitted.values + object$residuals
}
