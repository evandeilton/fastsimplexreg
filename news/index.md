# Changelog

## fastsimplexreg 0.2.1

Stability and performance release, following a three-way audit.

Stability / correctness:

- Non-convergence is now signalled:
  [`fastsimplexreg()`](https://evandeilton.github.io/fastsimplexreg/reference/fastsimplexreg.md)
  and
  [`fastsimplexregmixed()`](https://evandeilton.github.io/fastsimplexreg/reference/fastsimplexregmixed.md)
  emit a warning and withhold (NA) standard errors when the optimiser
  does not converge, and
  [`summary()`](https://rdrr.io/r/base/summary.html) prints a prominent
  banner. Previously a non-converged fit could return a
  confident-looking coefficient table.
- The mean-link inverse now keeps the analytic score and Hessian
  consistent with the objective at saturation (zeroing the derivatives
  when the mean clamps), fixing near-boundary line-search failures.
- The native BFGS declares soft convergence when the objective is
  already stationary but the gradient cannot be pushed below the
  tolerance (e.g. at the adaptive-quadrature noise floor).
- Malformed cluster offsets now raise a clean R error instead of
  aborting the process; the mixed finite-difference Hessian guards on
  evaluation validity (no more silent zero-variance columns); the
  `nAGQ^q` quadrature grid is capped.

Performance (mixed model):

- Single-pass adaptive quadrature with cached scores and hoisted
  constant terms, and tensor-weight pruning for `q >= 2`. A random-slope
  (`q = 2`) fit that took ~2 minutes now runs in a few seconds.
- Exact symmetrisation of the curvature matrix removes spurious Cholesky
  failures and redundant recomputation.

## fastsimplexreg 0.2.0

- New
  [`fastsimplexregmixed()`](https://evandeilton.github.io/fastsimplexreg/reference/fastsimplexregmixed.md):
  a two-level (nested) simplex mixed model with variable dispersion,
  estimated by adaptive Gauss-Hermite quadrature (AGHQ; `nAGQ = 1` gives
  the Laplace approximation). Gaussian random effects in the mean
  submodel are specified with an lme4-style `random = ~ terms | group`
  bar; the covariance is estimated on an unconstrained log-Cholesky
  scale. The per-cluster inner mode-finding, quadrature and analytic
  score run in C++ and are parallelised over clusters with OpenMP.
- S3 methods for class `"simplex_fast_mixed"`: `print`, `summary`,
  `coef`, `ranef`, `VarCorr` (re-exported from `nlme`), `vcov`,
  `logLik`, `nobs`, `ngrps`, `fitted`, `residuals`, `predict` (with
  `re.form`) and `plot`.
- The vignettes now analyse real data: the `sdac` (CD34+ cell recovery)
  and `retinal` (longitudinal intraocular gas) data sets from
  `simplexreg`.
- A `pkgdown` documentation website, GitHub Actions workflows
  (R-CMD-check, pkgdown, test-coverage), and an MIT license.

## fastsimplexreg 0.1.0

Initial release.

- [`fastsimplexreg()`](https://evandeilton.github.io/fastsimplexreg/reference/fastsimplexreg.md)
  fits simplex regression models with separate submodels for the mean
  and the dispersion via a multi-part `Formula` interface
  (`y ~ x1 + x2 | z1 + z2`), with `logit`, `probit`, `cloglog` and
  `neglog` mean links and a log dispersion link.
- Maximum-likelihood estimation uses an analytic score and a native BFGS
  optimiser implemented in C++ with `RcppArmadillo`, `BLAS`/`LAPACK` and
  optional `OpenMP` parallelism.
- [`dsimplex()`](https://evandeilton.github.io/fastsimplexreg/reference/dsimplex.md)
  and
  [`rsimplex()`](https://evandeilton.github.io/fastsimplexreg/reference/rsimplex.md)
  provide the density and random generation for the simplex
  distribution;
  [`simplex_linkinv()`](https://evandeilton.github.io/fastsimplexreg/reference/simplex_linkinv.md)
  exposes the mean link inverses.
- A full set of S3 methods for class `"simplex_fast"`, matching the
  conventions of `lm`/`glm`/`betareg` fits:
  [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html) (Wald),
  [`logLik()`](https://rdrr.io/r/stats/logLik.html) (so that
  [`AIC()`](https://rdrr.io/r/stats/AIC.html) and
  [`BIC()`](https://rdrr.io/r/stats/AIC.html) work),
  [`nobs()`](https://rdrr.io/r/stats/nobs.html),
  [`deviance()`](https://rdrr.io/r/stats/deviance.html),
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
  [`residuals()`](https://rdrr.io/r/stats/residuals.html) (response,
  Pearson and deviance),
  [`predict()`](https://rdrr.io/r/stats/predict.html),
  [`simulate()`](https://rdrr.io/r/stats/simulate.html),
  [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html),
  [`terms()`](https://rdrr.io/r/stats/terms.html),
  [`formula()`](https://rdrr.io/r/stats/formula.html),
  [`model.frame()`](https://rdrr.io/r/stats/model.frame.html),
  [`update()`](https://rdrr.io/r/stats/update.html) (multi-part-formula
  aware), [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html).
- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method
  producing `ggplot2` diagnostic panels (residuals vs fitted, normal
  Q-Q, scale-location, observed vs fitted), optionally combined with
  `patchwork` when installed.
