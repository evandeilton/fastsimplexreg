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
#' @aliases fastsimplexreg-package
#'
#' @useDynLib fastsimplexreg, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats pnorm qnorm qlogis model.frame model.matrix terms
#'   delete.response na.omit setNames printCoefmat .getXlevels ppoints quantile
#' @importFrom stats coef confint fitted logLik nobs predict residuals vcov
#'   deviance formula simulate update
#' @importFrom rlang .data
"_PACKAGE"
