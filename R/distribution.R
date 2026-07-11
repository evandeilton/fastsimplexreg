#' @title Simplex Distribution Density
#'
#' @description
#' Evaluates the probability density function of the simplex distribution of
#' Barndorff-Nielsen and Jorgensen (1991) for a mean `mu` and a dispersion
#' `phi` (the parameter often written \eqn{\sigma^2}). The density is
#' \deqn{f(x; \mu, \phi) = [2\pi\phi\,(x(1-x))^3]^{-1/2}
#'   \exp\!\left\{-\frac{1}{2\phi}\,
#'   \frac{(x-\mu)^2}{x(1-x)\,\mu^2(1-\mu)^2}\right\},
#'   \qquad 0 < x < 1.}
#' Values of `x` outside the open interval \eqn{(0, 1)} return `0` (or `-Inf`
#' on the log scale). The arguments `mu` and `phi` are recycled against `x`.
#' The computation is carried out in C++ and may use OpenMP threads.
#'
#' @param x Numeric vector of observations. Values must lie strictly inside
#'   \eqn{(0, 1)} to receive positive density.
#' @param mu Numeric vector of means in \eqn{(0, 1)}, of length one or
#'   `length(x)`.
#' @param phi Numeric vector of positive dispersion values, of length one or
#'   `length(x)`.
#' @param log Logical; if `TRUE`, log-densities are returned.
#' @param n_threads Integer number of OpenMP threads. Use `0` to request all
#'   threads available to the backend. Defaults to `1L` (serial).
#'
#' @return A numeric vector of densities (or log-densities when `log = TRUE`)
#'   with the recycled length of `x`, `mu` and `phi`.
#'
#' @references
#' Barndorff-Nielsen, O. E. and Jorgensen, B. (1991).
#' Some parametric models on the simplex.
#' *Journal of Multivariate Analysis*, **39**(1), 106--116.
#'
#' @seealso [rsimplex()], [fastsimplexreg()]
#'
#' @examples
#' dsimplex(c(0.2, 0.5, 0.8), mu = 0.5, phi = 1)
#' dsimplex(c(0.2, 0.5, 0.8), mu = 0.5, phi = 1, log = TRUE)
#'
#' # Integrates to one over the support.
#' integrate(function(u) dsimplex(u, mu = 0.4, phi = 2), 0, 1)$value
#'
#' @export
dsimplex <- function(x, mu, phi, log = FALSE, n_threads = 1L) {
  dsimplex_cpp(
    y = as.numeric(x),
    mu = as.numeric(mu),
    phi = as.numeric(phi),
    log = isTRUE(log),
    n_threads = as.integer(n_threads)
  )
}


#' @title Simplex Distribution Random Generation
#'
#' @description
#' Generates random deviates from the simplex distribution. The sampler is
#' implemented in C++ using an exact transformation based on an
#' inverse-Gaussian mixture representation. The arguments `mu` and `phi` are
#' recycled to length `n`.
#'
#' @param n Integer number of observations to generate.
#' @param mu Numeric vector of means in \eqn{(0, 1)}, of length one or `n`.
#' @param phi Numeric vector of positive dispersion values, of length one or
#'   `n`.
#'
#' @return A numeric vector of length `n` with values in \eqn{(0, 1)}.
#'
#' @references
#' Barndorff-Nielsen, O. E. and Jorgensen, B. (1991).
#' Some parametric models on the simplex.
#' *Journal of Multivariate Analysis*, **39**(1), 106--116.
#'
#' @seealso [dsimplex()], [fastsimplexreg()]
#'
#' @examples
#' set.seed(123)
#' y <- rsimplex(1000, mu = 0.35, phi = 0.8)
#' summary(y)
#'
#' @export
rsimplex <- function(n, mu, phi) {
  rsimplex_cpp(
    n = as.integer(n),
    mu = as.numeric(mu),
    phi = as.numeric(phi)
  )
}
