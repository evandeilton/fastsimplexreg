# fastsimplexreg 0.2.0

* New `fastsimplexregmixed()`: a two-level (nested) simplex mixed model with
  variable dispersion, estimated by adaptive Gauss-Hermite quadrature (AGHQ;
  `nAGQ = 1` gives the Laplace approximation). Gaussian random effects in the
  mean submodel are specified with an lme4-style `random = ~ terms | group`
  bar; the covariance is estimated on an unconstrained log-Cholesky scale. The
  per-cluster inner mode-finding, quadrature and analytic score run in C++ and
  are parallelised over clusters with OpenMP.
* S3 methods for class `"simplex_fast_mixed"`: `print`, `summary`, `coef`,
  `ranef`, `VarCorr` (re-exported from `nlme`), `vcov`, `logLik`, `nobs`,
  `ngrps`, `fitted`, `residuals`, `predict` (with `re.form`) and `plot`.

# fastsimplexreg 0.1.0

Initial release.

* `fastsimplexreg()` fits simplex regression models with separate submodels for
  the mean and the dispersion via a multi-part `Formula` interface
  (`y ~ x1 + x2 | z1 + z2`), with `logit`, `probit`, `cloglog` and `neglog`
  mean links and a log dispersion link.
* Maximum-likelihood estimation uses an analytic score and a native BFGS
  optimiser implemented in C++ with `RcppArmadillo`, `BLAS`/`LAPACK` and
  optional `OpenMP` parallelism.
* `dsimplex()` and `rsimplex()` provide the density and random generation for
  the simplex distribution; `simplex_linkinv()` exposes the mean link inverses.
* A full set of S3 methods for class `"simplex_fast"`, matching the conventions
  of `lm`/`glm`/`betareg` fits: `coef()`, `vcov()`, `confint()` (Wald),
  `logLik()` (so that `AIC()` and `BIC()` work), `nobs()`, `deviance()`,
  `fitted()`, `residuals()` (response, Pearson and deviance), `predict()`,
  `simulate()`, `model.matrix()`, `terms()`, `formula()`, `model.frame()`,
  `update()` (multi-part-formula aware), `print()` and `summary()`.
* `plot()` method producing `ggplot2` diagnostic panels (residuals vs fitted,
  normal Q-Q, scale-location, observed vs fitted), optionally combined with
  `patchwork` when installed.
