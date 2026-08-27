# fastsimplexreg 0.2.4

Correctness release following a four-way consensus audit of 0.2.3
(mathematical/statistical rigour, R API contracts, CRAN compliance). Every item
below carries a measured proof in `tests/testthat`.

## Silently wrong results (fixed)

* **`confint()` on a mixed fit reported another parameter's interval.** There
  was no `confint.simplex_fast_mixed`, so dispatch fell through to
  `stats::confint.default`, which indexes `vcov()` by *name* -- and `coef()` of
  a mixed fit repeated `"(Intercept)"` across the mean and dispersion
  submodels. Both rows resolved to the first match: the dispersion intercept
  was reported with the mean intercept's interval, one that need not even
  contain its own estimate, while the variance components vanished from the
  table. There is now a method, selecting by position over the full parameter
  vector.
* **`offset()` in the formula was silently dropped.** The `terms` object
  recorded it, the design matrix omitted it, and nothing warned, so the fitted
  model was not the model the user wrote. Offsets are now honoured in both
  fitters, kept separate per submodel, and rebuilt from `newdata` in
  `predict()`.
* **`subset` was an index vector, not an expression.** The documentation
  reproduced `glm`'s wording. `subset = x1 > 0` errored loudly when no `x1`
  existed in the caller, but silently used *that* vector when one did --
  fitting the wrong rows with no symptom. It is now evaluated inside `data`
  first, as in `stats::lm()`. Every previously working call form still works.
* **Adaptive quadrature discarded the mass that mattered.** The AGHQ node
  pruning tested the product Gauss-Hermite weight, but the adaptive transform
  undoes the `e^{-t^2}` factor, so a node's real multiplier is `logW + t2`.
  Pruning discarded 5.7% of the effective quadrature mass at `nAGQ = 11` and
  48.7% at `nAGQ = 21`, which made the AGHQ sequence *stop converging*: raising
  `nAGQ` moved the marginal log-likelihood away from its limit. Pruning is
  removed.
* **Rank-deficient designs returned an arbitrary split.** With `x2 = 2 * x1`
  the fit reported `x1 = 0.105` and `x2 = 0.210`, two numbers that mean nothing
  individually. Aliased columns are now detected before fitting by the same
  pivoted QR `lm()` uses, and reported as `NA`.
* **`deviance()` was identically `nobs`.** The score equation for the
  dispersion submodel forces `sum(dev_i/phi_i) = n` whenever it has an
  intercept, so the reported "Deviance" was the same number for every model.
  The default is now the unscaled deviance; `type = "scaled"` still gives the
  old quantity, documented.

## Statistical accuracy

* **Exact Fisher information.** `I(mu) = 1/(phi mu^3 (1-mu)^3) + 3/(mu(1-mu))`;
  the second term was missing, understating the information by up to 4.75x.
  New `information = c("observed", "expected")` argument on `fastsimplexreg()`
  offers the closed-form expected information, which is exactly block diagonal
  in `(beta, gamma)` and positive definite by construction. The default stays
  `"observed"`, following Efron and Hinkley (1978).
* **`psimplex()` is a closed form.** The adaptive Gauss-Legendre quadrature is
  replaced by an exact expression, agreeing with seeded numerical integration
  to 1.8e-13. `log.p` is now computed on the log scale throughout:
  `psimplex(0.15, 0.5, 0.01, log.p = TRUE)` returned `-Inf` where the true
  value is `-773.216`.
* **Randomized quantile residuals, and they are now the default.** Under a
  correct model, Pearson residuals rejected up to 100% of fits in a
  Shapiro-Wilk check -- and Pearson was what `summary()` printed as its
  residual summary. Quantile residuals hold their nominal rate.
* **The AGHQ scaling matrix is ridged, not switched.** Replacing the observed
  curvature outright by the Fisher information made the objective
  discontinuous in `theta`, by 0.04 to 2.0 nats.

## Contracts and guards

* `na.action = na.exclude` now pads `fitted()`, `residuals()` and `predict()`
  back to `nrow(data)`, as `glm` does; it was accepted and ignored.
* Dispersion coefficients are named `(phi)_*` in the full parameter vector, as
  in **betareg**, so `vcov()` and `confint()` can be indexed by name
  unambiguously. The per-submodel tables in `summary()` keep bare names.
* The packed `omega` diagonal is labelled `logchol.*` for `q >= 2`: it is the
  *conditional* standard deviation there, not the marginal one, and reading it
  as an SD understated a random slope by a factor of 1.95.
* `fitted()`, `residuals()` and `predict()` carry the observation labels.
* `nAGQ < 5` warns; `inner_maxit < 10` is refused; a single-level grouping
  factor is refused and an all-singleton design warns.
* `ngrps()` no longer breaks when **lme4** is attached afterwards.
* `summary()` now reports the diagnostics the fit object already stored: vcov
  rank, pseudo-inverse use, ill-conditioning, saturated observations, aliased
  coefficients, and whether inference was computed at all.

## Packaging

* `Depends: R (>= 4.0.0)` -- no R 4.1 feature is used.
* `simplexreg`, `microbenchmark`, `parallel`, `utils` and `MASS` declared in
  `Suggests`; the shipped benchmark script needs them.
* Added `inst/CITATION` and a test-coverage workflow; dropped the redundant
  `SystemRequirements: C++17`; fixed a 404 badge in the README.

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

## Mixed-model convergence

* **Roughly a quarter of mixed fits were reported as optimiser failures when
  they had in fact converged**, and because standard errors are only computed at
  a converged fit, those perfectly good fits silently lost their inference.
  Measured over a grid of 240 `(J, nj, nAGQ, seed)` combinations, 26.2% ended
  with code 2; restarting the optimiser from the reported stopping point gained
  a median of 1.6e-11 in log-likelihood -- it was an optimum, not a failure.
  The failure rate is now 0% over the same grid.

  The rate depends on `nAGQ`, not on the number of groups: 55% at `nAGQ = 3`,
  47% at 5, 3% at 7 and 0% at the default 11, while being flat in `J` (25-29%)
  and in cluster size. The cause is that the analytic score is the exact score
  of the *true* marginal likelihood (Fisher's identity), not of its `nAGQ`-point
  quadrature approximation, so `grad_tol` is unreachable when `nAGQ` is small
  and the run ends on a line-search failure instead.

  Two changes: when a line search fails, the inverse-Hessian approximation is
  reset and the iteration is retried once from a clean steepest-descent
  direction (it goes stale, and this recovers real progress -- it moved 59 of
  180 stopping points, every one of them to a *higher* log-likelihood); and the
  soft-convergence test, which required a relative change below `rel_tol`
  (1e-9), now uses `sqrt(rel_tol)`, the scale at which these runs actually
  flatten out. A genuine failure -- a run that cannot take a single step -- is
  still reported as code 2.

  Cost: the fixed-effects path is unaffected in results and essentially
  unaffected in speed (+3.2% function evaluations, in 3 of 24 fits, all
  `neglog`; wall clock within noise). Mixed fits that already converged use
  about 50% more function evaluations, which buys the extra accuracy above; fits
  that used to abort naturally cost more now that they run to completion.

* Known limitation: supplying `start` equal to the optimum itself gives the
  optimiser no iteration history from which to judge stationarity, so it
  conservatively reports code 2 rather than risk labelling a stall as success.

* The optimiser trace (`trace = TRUE`) now prints the objective at full double
  precision and reports the relative change per iteration. At the previous six
  significant digits, successive iterations near the optimum printed
  identically -- exactly the regime the trace exists to diagnose.

## Testing and documentation

* New `test-inference.R` and `test-distribution-conventions.R`, plus regression
  tests pinning the convergence-reporting behaviour in both directions.
* The rank cut-off used by the fail-safe covariance is derived from the accuracy
  of the finite-difference Hessian (`sqrt(eps)` on the correlation scale) rather
  than from machine epsilon. With the tighter cut the verdict depended on which
  BLAS computed the Hessian: an exactly collinear design was flagged on Linux
  and passed as full rank on Windows.
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
