# fastsimplexreg 0.2.3

Correctness release following a full surgical audit of `R/`, `src/`, `tests/`
and `vignettes/`. All estimates, log-likelihoods, deviances, fitted values,
residuals, densities and random draws are numerically unchanged from 0.2.2;
what changed is what the package does when something is *wrong*.

## Breaking changes

* **The distribution functions now recycle instead of erroring.** `dsimplex()`
  and `rsimplex()` used to require `mu`/`phi` of length 1 or `length(x)`; they
  now recycle all arguments to their common maximum length, which is the base-R
  convention (`dnorm(1:3, mean = 1:2)`).
* **`NA` and invalid parameters no longer collapse to `0`.** `dsimplex()`
  returned `0` (`-Inf` on the log scale) for `NA` input and for an out-of-domain
  `mu`/`phi`, turning "missing" into "impossible". `NA` now propagates as `NA`,
  `NaN` as `NaN`, and an out-of-domain parameter gives `NaN` with the canonical
  "NaNs produced" warning. `rsimplex()` gives `NaN` with "NAs produced" instead
  of aborting the call.

## Inference (critical)

* **Standard errors are never reported as a confident zero again.** With
  collinear or otherwise unidentified covariates the Hessian could be singular
  or indefinite, `solve()` would still succeed, and `sqrt(pmax(diag(vcov), 0))`
  turned the resulting negative variances into `SE = 0` -- so `summary()`
  printed `z = Inf` and `p = 0` for a parameter the data cannot identify, with
  no warning. Inversion now goes through an eigen-decomposition of the
  equilibrated (correlation-scale) information matrix: only the strictly
  positive spectrum is inverted (a Moore-Penrose pseudo-inverse), directions of
  negative curvature are refused rather than turned into variances, every
  parameter loading on a discarded direction gets `SE = NA`, and the affected
  parameters are named in a warning. A merely ill-conditioned but full-rank fit
  keeps its (large, honest) standard errors and gets a separate warning.
* New components on the fitted object: `vcov_rank`, `vcov_pseudo`,
  `vcov_condition`, `vcov_eigenvalues` and `n_saturated`.
* **Saturation of the mean link is reported.** The likelihood path floors the
  mean at `1e-12`, which zeroes those observations' contribution to the score.
  That now raises a warning naming how many observations are affected.

## New

* `psimplex()` and `qsimplex()` complete the `d`/`p`/`q`/`r` family, with
  `lower.tail` and `log.p`. The CDF uses adaptive Gauss-Legendre quadrature with
  panels seeded around the mean (so a sharply peaked density is always
  resolved) and agrees with `stats::integrate()` to ~1e-15; `qsimplex()` inverts
  it by safeguarded Newton-bisection.

## Prediction

* `predict()` now preserves the length of `newdata`: rows dropped for
  missingness come back as `NA` instead of silently yielding a shorter,
  unaligned vector. All model parts are built from the same set of complete
  rows, so they can no longer end up with different row counts.
* `predict()` on a mixed fit **errors** when `newdata` lacks the grouping
  column. It previously returned population-level predictions under the label
  of conditional ones, silently. Rows whose group level was not seen in the fit
  still fall back to a zero random effect, but that substitution is now
  announced.
* Population-level prediction on a mixed fit stored with `model = FALSE` errors
  instead of resolving covariates in the caller's environment, where it could
  return predictions of the wrong length built from unrelated objects.
* A variable required by the model but absent from `newdata` is an error.

## Numerical / C++

* The reporting path (`simplex_linkinv()`, `predict()`, `fitted()`) no longer
  applies the likelihood path's `1e-12` floor. It is clamped only at the
  representable boundary, so fitted means keep their full dynamic range
  (`simplex_linkinv(-40, "logit")` is `4.2e-18`, not `1e-12`) while remaining
  strictly inside the open support `(0, 1)` that the density and the residual
  formulas require.
* `simplex_mixed_ranef_cpp()` no longer reads uninitialised memory: the Fisher
  fallback for the posterior covariance consumed a buffer that the observed
  loop could leave partially filled, and the return value of the link map was
  discarded.
* The mixed model's inner mode solver no longer commits a step its line search
  rejected, and no longer compares against the objective at an inadmissible
  point. A warm-started mode inherited from an outer trial point that was later
  rejected is now discarded when the cold start beats it, removing a path
  dependence in the objective.
* `VarCorr()`'s `sigma` argument is inert for a simplex mixed model; supplying
  anything other than `1` now warns instead of being silently ignored.

## Testing and documentation

* New `test-inference.R` and `test-distribution-conventions.R`.
* The analytic score is now compared directly against `numDeriv` for all four
  mean links (previously validated only indirectly, and only for `logit`).
* The benchmark vignette's accuracy table is generated from the shipped results
  instead of being hard-coded, so it cannot drift from the figures.
* `simulate()` now carries the row names of the model frame (the previous code
  read `names(object$y)`, which is always `NULL`).

# fastsimplexreg 0.2.2

* New `benchmark` vignette comparing `fastsimplexreg` with the CRAN packages
  `simplexreg` and `betareg` on accuracy, speed (scaling to n = 5e5), the four
  mean links, real data, and the mixed model, with professional ggplot2 figures.
  Estimates are numerically identical to `simplexreg` while fitting ~5-13x faster.
  The full reproducible study ships in `inst/benchmark/`.

# fastsimplexreg 0.2.1

Stability and performance release, following a three-way audit.

Stability / correctness:

* Non-convergence is now signalled: `fastsimplexreg()` and
  `fastsimplexregmixed()` emit a warning and withhold (NA) standard errors when
  the optimiser does not converge, and `summary()` prints a prominent banner.
  Previously a non-converged fit could return a confident-looking coefficient
  table.
* The mean-link inverse now keeps the analytic score and Hessian consistent with
  the objective at saturation (zeroing the derivatives when the mean clamps),
  fixing near-boundary line-search failures.
* The native BFGS declares soft convergence when the objective is already
  stationary but the gradient cannot be pushed below the tolerance (e.g. at the
  adaptive-quadrature noise floor).
* Malformed cluster offsets now raise a clean R error instead of aborting the
  process; the mixed finite-difference Hessian guards on evaluation validity
  (no more silent zero-variance columns); the `nAGQ^q` quadrature grid is capped.

Performance (mixed model):

* Single-pass adaptive quadrature with cached scores and hoisted constant terms,
  and tensor-weight pruning for `q >= 2`. A random-slope (`q = 2`) fit that took
  ~2 minutes now runs in a few seconds.
* Exact symmetrisation of the curvature matrix removes spurious Cholesky
  failures and redundant recomputation.

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
* The vignettes now analyse real data: the `sdac` (CD34+ cell recovery) and
  `retinal` (longitudinal intraocular gas) data sets from `simplexreg`.
* A `pkgdown` documentation website, GitHub Actions workflows (R-CMD-check,
  pkgdown, test-coverage), and an MIT license.

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
