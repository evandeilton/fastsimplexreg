# Internal mean-link registry. Integer identifiers avoid string matching inside
# the C++ hot path and are shared, verbatim, across the R/C++ boundary:
# logit = 1, probit = 2, cloglog = 3, neglog = 4.
.simplex_links <- c(logit = 1L, probit = 2L, cloglog = 3L, neglog = 4L)

# Internal helper: validate a mean-link name and return its name and integer id.
.normalize_simplex_link <- function(link) {
  link <- match.arg(link, names(.simplex_links))
  list(name = link, id = unname(.simplex_links[[link]]))
}


#' @title Inverse of the Simplex Mean Link
#'
#' @description
#' Computes the mean \eqn{\mu = g^{-1}(\eta)} from a linear predictor `eta`,
#' using the same C++ backend employed when fitting a model with
#' [fastsimplexreg()]. Four links are supported:
#' \describe{
#'   \item{`logit`}{\eqn{g^{-1}(\eta) = 1 / (1 + e^{-\eta})}.}
#'   \item{`probit`}{\eqn{g^{-1}(\eta) = \Phi(\eta)}.}
#'   \item{`cloglog`}{\eqn{g^{-1}(\eta) = 1 - \exp(-\exp(\eta))}.}
#'   \item{`neglog`}{\eqn{g^{-1}(\eta) = \exp(-\exp(-\eta))}, following the
#'     definition in Zhang et al. (2016).}
#' }
#'
#' @param eta Numeric vector; the linear predictor for the mean.
#' @param link Character string selecting the mean link. One of `"logit"`,
#'   `"probit"`, `"cloglog"` or `"neglog"`.
#'
#' @return A numeric vector of means in \eqn{(0, 1)}, of the same length as
#'   `eta`.
#'
#' @references
#' Zhang, P., Qiu, Z. and Shi, C. (2016).
#' simplexreg: An R Package for Regression Analysis of Proportional Data Using
#' the Simplex Distribution.
#' *Journal of Statistical Software*, **71**(11), 1--21.
#'
#' @seealso [fastsimplexreg()]
#'
#' @examples
#' eta <- seq(-3, 3, length.out = 7)
#' simplex_linkinv(eta, link = "logit")
#' simplex_linkinv(eta, link = "probit")
#' simplex_linkinv(eta, link = "cloglog")
#' simplex_linkinv(eta, link = "neglog")
#'
#' @export
simplex_linkinv <- function(eta, link = c("logit", "probit", "cloglog", "neglog")) {
  link_spec <- .normalize_simplex_link(link)
  simplex_linkinv_cpp(as.numeric(eta), mean_link = link_spec$id)
}
