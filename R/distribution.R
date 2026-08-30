# distribution.R
# The simplex distribution family, following base R's d/p/q/r conventions.
#
# Conventions honoured throughout (see ?dnorm, ?dbeta for the reference
# behaviour these mirror):
#   * `x`/`q`/`p`, `mu` and `phi` are recycled to their common maximum length;
#   * NA propagates as NA and NaN as NaN;
#   * a parameter outside its domain (mu outside (0,1), phi <= 0, or non-finite)
#     yields NaN with a single "NaNs produced" warning -- for rsimplex(), NaN
#     with an "NAs produced" warning, exactly as rbeta() does;
#   * `log`, `log.p` and `lower.tail` behave as in base R.

# Internal: apply the base-R "NaNs/NAs produced" warning convention. The C++
# kernels attach the number of out-of-domain parameter combinations as the
# attribute "n_invalid_par"; this strips it and raises the canonical warning.
.simplex_dist_finish <- function(values, msg = "NaNs produced") {
  n_bad <- attr(values, "n_invalid_par")
  attr(values, "n_invalid_par") <- NULL
  if (!is.null(n_bad) && isTRUE(n_bad > 0L)) {
    warning(msg, call. = FALSE)
  }
  values
}


#' @title The Simplex Distribution
#'
#' @description
#' Density, distribution function, quantile function and random generation for
#' the simplex distribution of Barndorff-Nielsen and Jorgensen (1991), with mean
#' `mu` and dispersion `phi` (the parameter often written \eqn{\sigma^2}). The
#' density is
#' \deqn{f(x; \mu, \phi) = [2\pi\phi\,(x(1-x))^3]^{-1/2}
#'   \exp\!\left\{-\frac{1}{2\phi}\,
#'   \frac{(x-\mu)^2}{x(1-x)\,\mu^2(1-\mu)^2}\right\},
#'   \qquad 0 < x < 1.}
#'
#' @details
#' These functions follow the conventions of base R's distribution family:
#' `x`/`q`/`p`, `mu` and `phi` are recycled to their common length; `NA`
#' propagates as `NA` and `NaN` as `NaN`; and a parameter outside its domain
#' (`mu` outside \eqn{(0,1)}, `phi` not positive, or either non-finite) produces
#' `NaN` with a warning, rather than an error or a silent zero. Values of `x`
#' outside the open support \eqn{(0, 1)} have density `0` (`-Inf` on the log
#' scale), which is a genuine density value and is therefore not a warning.
#'
#' The distribution function is available in **closed form**. Mapping to the
#' odds scale \eqn{X = Y/(1-Y)} turns the simplex density into a mixture of an
#' inverse Gaussian and its size-biased version, both of which integrate
#' exactly, giving
#' \deqn{F(y; \mu, \phi) = \Phi(a) + (1 - 2\mu)\,e^{k}\,\Phi(b),}
#' with
#' \deqn{a = \frac{y - \mu}{\mu(1-\mu)\sqrt{\phi\, y (1-y)}}, \qquad
#'   b = \frac{-(y + \mu - 2 y \mu)}{\mu(1-\mu)\sqrt{\phi\, y (1-y)}},
#'   \qquad k = \frac{2}{\phi\, \mu (1-\mu)}.}
#' The product \eqn{e^{k}\Phi(b)} is formed on the log scale and never
#' evaluated directly, so it neither overflows nor loses `log.p`:
#' `psimplex(0.15, mu = 0.5, phi = 0.01, log.p = TRUE)` returns about `-773.2`
#' rather than `-Inf`. The upper tail uses the exact reflection
#' \eqn{P(Y > y \mid \mu) = F(1 - y \mid 1 - \mu)}, which follows from the unit
#' deviance satisfying \eqn{d(y; \mu) = d(1-y; 1-\mu)} and keeps full relative
#' accuracy in both tails.
#'
#' `qsimplex()` inverts `psimplex()` by safeguarded Newton-bisection. Both are
#' more expensive than `dsimplex()` and accept `n_threads` for that reason.
#'
#' `rsimplex()` uses the exact inverse-Gaussian-mixture representation: with
#' \eqn{\epsilon = \mu/(1-\mu)} and \eqn{\tau = \phi (1-\mu)^2}, a variate
#' \eqn{x} is built from an inverse-Gaussian draw plus, with probability
#' \eqn{\mu}, a chi-squared(1) term, and mapped back to \eqn{(0,1)} through
#' \eqn{x/(1+x)}.
#'
#' @param x,q Numeric vector of quantiles.
#' @param p Numeric vector of probabilities.
#' @param n Number of observations to generate. If `length(n) > 1`, the length
#'   is taken to be the number required (the base-R convention).
#' @param mu Numeric vector of means in \eqn{(0, 1)}.
#' @param phi Numeric vector of positive dispersion values.
#' @param log,log.p Logical; if `TRUE`, probabilities/densities are given as
#'   \eqn{\log(p)}.
#' @param lower.tail Logical; if `TRUE` (default), probabilities are
#'   \eqn{P(X \le x)}, otherwise \eqn{P(X > x)}.
#' @param n_threads Integer number of OpenMP threads. Use `0` to request all
#'   threads available to the backend. Defaults to `1L` (serial).
#'
#' @return `dsimplex()` gives the density, `psimplex()` the distribution
#'   function, `qsimplex()` the quantile function, and `rsimplex()` generates
#'   random deviates. The length of the result of `rsimplex()` is `n`; for the
#'   other functions it is the maximum of the lengths of the numeric arguments.
#'
#' @references
#' Barndorff-Nielsen, O. E. and Jorgensen, B. (1991).
#' Some parametric models on the simplex.
#' *Journal of Multivariate Analysis*, **39**(1), 106--116.
#' \doi{10.1016/0047-259X(91)90008-P}
#'
#' @seealso [fastsimplexreg()]
#'
#' @examples
#' dsimplex(c(0.2, 0.5, 0.8), mu = 0.5, phi = 1)
#' dsimplex(c(0.2, 0.5, 0.8), mu = 0.5, phi = 1, log = TRUE)
#'
#' # Integrates to one over the support.
#' psimplex(1, mu = 0.4, phi = 2)
#'
#' # q is the inverse of p.
#' psimplex(qsimplex(c(0.1, 0.5, 0.9), mu = 0.4, phi = 2), mu = 0.4, phi = 2)
#'
#' set.seed(123)
#' y <- rsimplex(1000, mu = 0.35, phi = 0.8)
#' summary(y)
#'
#' @name simplex-distribution
#' @rdname simplex-distribution
#' @export
dsimplex <- function(x, mu, phi, log = FALSE, n_threads = 1L) {
  .simplex_dist_finish(dsimplex_cpp(
    y = as.numeric(x),
    mu = as.numeric(mu),
    phi = as.numeric(phi),
    log = isTRUE(log),
    n_threads = as.integer(n_threads)
  ))
}


#' @rdname simplex-distribution
#' @export
psimplex <- function(q, mu, phi, lower.tail = TRUE, log.p = FALSE,
                     n_threads = 1L) {
  .simplex_dist_finish(psimplex_cpp(
    q = as.numeric(q),
    mu = as.numeric(mu),
    phi = as.numeric(phi),
    lower_tail = isTRUE(lower.tail),
    log_p = isTRUE(log.p),
    n_threads = as.integer(n_threads)
  ))
}


#' @rdname simplex-distribution
#' @export
qsimplex <- function(p, mu, phi, lower.tail = TRUE, log.p = FALSE,
                     n_threads = 1L) {
  .simplex_dist_finish(qsimplex_cpp(
    p = as.numeric(p),
    mu = as.numeric(mu),
    phi = as.numeric(phi),
    lower_tail = isTRUE(lower.tail),
    log_p = isTRUE(log.p),
    n_threads = as.integer(n_threads)
  ))
}


#' @rdname simplex-distribution
#' @export
rsimplex <- function(n, mu, phi) {
  # Base-R convention: rnorm(c(1, 2, 3)) generates 3 deviates.
  if (length(n) > 1L) n <- length(n)
  n <- as.numeric(n)
  if (length(n) != 1L || is.na(n) || n < 0) {
    stop("'n' must be a single non-negative number.", call. = FALSE)
  }
  .simplex_dist_finish(
    rsimplex_cpp(n = as.numeric(n), mu = as.numeric(mu), phi = as.numeric(phi)),
    msg = "NAs produced"
  )
}
