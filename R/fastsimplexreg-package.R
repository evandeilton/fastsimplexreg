#' @title fastsimplexreg: Fast Simplex Regression with Variable Dispersion
#'
#' @description
#' High-performance maximum-likelihood estimation of simplex regression models
#' for continuous proportions in the open interval \eqn{(0, 1)}. The package
#' supports separate submodels for the mean and the dispersion through a
#' multi-part [Formula::Formula] interface (`y ~ x1 + x2 | z1 + z2`), four mean
#' links (`logit`, `probit`, `cloglog`, `neglog`) and a log link for the
#' dispersion. The entire numerical hot path -- log-likelihood, analytic score,
#' native BFGS optimiser, density, random generation, prediction and link
#' inverses -- is implemented in C++ with RcppArmadillo, BLAS/LAPACK and
#' optional OpenMP parallelism, so that models scale to large data sets.
#'
#' @references
#' Barndorff-Nielsen, O. E. and Jorgensen, B. (1991).
#' Some parametric models on the simplex.
#' *Journal of Multivariate Analysis*, **39**(1), 106--116.
#' \doi{10.1016/0047-259X(91)90008-P}
#'
#' Jorgensen, B. (1997). *The Theory of Dispersion Models*.
#' Chapman & Hall, London.
#'
#' Zhang, P., Qiu, Z. and Shi, C. (2016).
#' simplexreg: An R Package for Regression Analysis of Proportional Data Using
#' the Simplex Distribution.
#' *Journal of Statistical Software*, **71**(11), 1--21.
#' \doi{10.18637/jss.v071.i11}
#'
#' @seealso [fastsimplexreg()], [dsimplex()], [rsimplex()], [simplex_linkinv()]
#'
#' @keywords internal
#'
#' @useDynLib fastsimplexreg, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats pnorm qnorm qlogis model.frame model.matrix terms
#' @importFrom stats delete.response na.omit setNames printCoefmat
#' @importFrom stats .getXlevels ppoints quantile complete.cases
#' @importFrom stats coef confint fitted logLik nobs predict residuals vcov
#' @importFrom stats deviance formula simulate update
#' @importFrom rlang .data
"_PACKAGE"


# Register ngrps.simplex_fast_mixed with lme4's INDEPENDENT ngrps generic.
#
# lme4 defines its own ngrps generic (with an ngrps.default that stops), so
# attaching lme4 after fastsimplexreg masked ours and ngrps(fit) failed with
# "Cannot extract the number of groups from this object". Exporting the method
# is not enough: UseMethod() consults the S3 registration table of the namespace
# where the GENERIC is defined, which is lme4's, and that table cannot contain a
# method we never registered into it. So we register it -- immediately if lme4
# is already loaded, and through a load hook otherwise.
#
# ranef and VarCorr need none of this: lme4 re-exports the very same nlme
# generic objects we import, so there is only ever one generic in play.
.register_lme4_ngrps <- function() {
  if (!isNamespaceLoaded("lme4")) return(invisible(FALSE))
  try(registerS3method("ngrps", "simplex_fast_mixed",
                       ngrps.simplex_fast_mixed,
                       envir = asNamespace("lme4")), silent = TRUE)
  invisible(TRUE)
}

.onLoad <- function(libname, pkgname) {
  .register_lme4_ngrps()
  setHook(packageEvent("lme4", "onLoad"),
          function(...) .register_lme4_ngrps())
  invisible()
}
