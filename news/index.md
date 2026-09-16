# Changelog

## fastsimplexreg 0.2.4

Correctness release following a four-way consensus audit of 0.2.3
(mathematical/statistical rigour, R API contracts, CRAN compliance).
Every item below carries a measured proof in `tests/testthat`.

### Silently wrong results (fixed)

- **[`confint()`](https://rdrr.io/r/stats/confint.html) on a mixed fit
  reported another parameter’s interval.** There was no
  `confint.simplex_fast_mixed`, so dispatch fell through to
  [`stats::confint.default`](https://rdrr.io/r/stats/confint.html),
  which indexes [`vcov()`](https://rdrr.io/r/stats/vcov.html) by *name*
  – and [`coef()`](https://rdrr.io/r/stats/coef.html) of a mixed fit
  repeated `"(Intercept)"` across the mean and dispersion submodels.
  Both rows resolved to the first match: the dispersion intercept was
  reported with the mean intercept’s interval, one that need not even
  contain its own estimate, while the variance components vanished from
  the table. There is now a method, selecting by position over the full
  parameter vector.
- **[`offset()`](https://rdrr.io/r/stats/offset.html) in the formula was
  silently dropped.** The `terms` object recorded it, the design matrix
  omitted it, and nothing warned, so the fitted model was not the model
  the user wrote. Offsets are now honoured in both fitters, kept
  separate per submodel, and rebuilt from `newdata` in
  [`predict()`](https://rdrr.io/r/stats/predict.html).
- **`subset` was an index vector, not an expression.** The documentation
  reproduced `glm`’s wording. `subset = x1 > 0` errored loudly when no
  `x1` existed in the caller, but silently used *that* vector when one
  did – fitting the wrong rows with no symptom. It is now evaluated
  inside `data` first, as in
  [`stats::lm()`](https://rdrr.io/r/stats/lm.html). Every previously
  working call form still works.
- **Adaptive quadrature discarded the mass that mattered.** The AGHQ
  node pruning tested the product Gauss-Hermite weight, but the adaptive
  transform undoes the `e^{-t^2}` factor, so a node’s real multiplier is
  `logW + t2`. Pruning discarded 5.7% of the effective quadrature mass
  at `nAGQ = 11` and 48.7% at `nAGQ = 21`, which made the AGHQ sequence
  *stop converging*: raising `nAGQ` moved the marginal log-likelihood
  away from its limit. Pruning is removed.
- **Rank-deficient designs returned an arbitrary split.** With
  `x2 = 2 * x1` the fit reported `x1 = 0.105` and `x2 = 0.210`, two
  numbers that mean nothing individually. Aliased columns are now
  detected before fitting by the same pivoted QR
  [`lm()`](https://rdrr.io/r/stats/lm.html) uses, and reported as `NA`.
- **[`deviance()`](https://rdrr.io/r/stats/deviance.html) was
  identically `nobs`.** The score equation for the dispersion submodel
  forces `sum(dev_i/phi_i) = n` whenever it has an intercept, so the
  reported “Deviance” was the same number for every model. The default
  is now the unscaled deviance; `type = "scaled"` still gives the old
  quantity, documented.

### Statistical accuracy

- **Exact Fisher information.**
  `I(mu) = 1/(phi mu^3 (1-mu)^3) + 3/(mu(1-mu))`; the second term was
  missing, understating the information by up to 4.75x. New
  `information = c("observed", "expected")` argument on
  [`fastsimplexreg()`](https://evandeilton.github.io/fastsimplexreg/reference/fastsimplexreg.md)
  offers the closed-form expected information, which is exactly block
  diagonal in `(beta, gamma)` and positive definite by construction. The
  default stays `"observed"`, following Efron and Hinkley (1978).
- **[`psimplex()`](https://evandeilton.github.io/fastsimplexreg/reference/simplex-distribution.md)
  is a closed form.** The adaptive Gauss-Legendre quadrature is replaced
  by an exact expression, agreeing with seeded numerical integration to
  1.8e-13. `log.p` is now computed on the log scale throughout:
  `psimplex(0.15, 0.5, 0.01, log.p = TRUE)` returned `-Inf` where the
  true value is `-773.216`.
- **Randomized quantile residuals, and they are now the default.** Under
  a correct model, Pearson residuals rejected up to 100% of fits in a
  Shapiro-Wilk check – and Pearson was what
  [`summary()`](https://rdrr.io/r/base/summary.html) printed as its
  residual summary. Quantile residuals hold their nominal rate.
- **The AGHQ scaling matrix is ridged, not switched.** Replacing the
  observed curvature outright by the Fisher information made the
  objective discontinuous in `theta`, by 0.04 to 2.0 nats.

### Contracts and guards

- `na.action = na.exclude` now pads
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
  [`residuals()`](https://rdrr.io/r/stats/residuals.html) and
  [`predict()`](https://rdrr.io/r/stats/predict.html) back to
  `nrow(data)`, as `glm` does; it was accepted and ignored.
- Dispersion coefficients are named `(phi)_*` in the full parameter
  vector, as in **betareg**, so
  [`vcov()`](https://rdrr.io/r/stats/vcov.html) and
  [`confint()`](https://rdrr.io/r/stats/confint.html) can be indexed by
  name unambiguously. The per-submodel tables in
  [`summary()`](https://rdrr.io/r/base/summary.html) keep bare names.
- The packed `omega` diagonal is labelled `logchol.*` for `q >= 2`: it
  is the *conditional* standard deviation there, not the marginal one,
  and reading it as an SD understated a random slope by a factor of
  1.95.
- [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
  [`residuals()`](https://rdrr.io/r/stats/residuals.html) and
  [`predict()`](https://rdrr.io/r/stats/predict.html) carry the
  observation labels.
- `nAGQ < 5` warns; `inner_maxit < 10` is refused; a single-level
  grouping factor is refused and an all-singleton design warns.
- [`ngrps()`](https://evandeilton.github.io/fastsimplexreg/reference/ngrps.md)
  no longer breaks when **lme4** is attached afterwards.
- [`summary()`](https://rdrr.io/r/base/summary.html) now reports the
  diagnostics the fit object already stored: vcov rank, pseudo-inverse
  use, ill-conditioning, saturated observations, aliased coefficients,
  and whether inference was computed at all.

### Parallelism and the C++ backend (audit of the OpenMP layer)

- **The linear predictors are no longer computed through BLAS.**
  `X %*% beta` and `Z %*% gamma` were two `gemv` calls, serial with
  respect to the package’s own OpenMP region, and a threaded BLAS made
  them *slower*: a 1e6 x 2 matvec took 14.9 ms across 24 OpenBLAS
  threads against 5.4 ms pinned to one. They were also an Amdahl ceiling
  – 24% of an evaluation on one thread but 71% on sixteen. `eta` is now
  accumulated inside the existing parallel loop, in the same column
  order a column-major `gemv` uses, and the single-threaded result is
  **bit-identical** to the BLAS version across every `p` and link
  tested. End-to-end speed-ups on a 24-core machine:

  | problem         | 2 threads | 4    | 8    | 16   |
  |-----------------|-----------|------|------|------|
  | n = 1e5, p = 5  | 1.4x      | 2.2x | 2.1x | 2.2x |
  | n = 1e6, p = 10 | 1.9x      | 3.2x | 5.0x | 5.2x |
  | n = 5e6, p = 10 | 1.1x      | 3.3x | 5.1x | 5.1x |

  The same n = 1e5 case was previously a *loss* at every thread count
  under default BLAS settings (0.22x at 16 threads).

- **Long calls can be interrupted.** There was no interrupt check
  anywhere in the C++ backend, so `Ctrl-C` was dead for the whole of a
  `.Call` – measured at 154 s for a single-threaded fit at n = 2e6.
  `Rcpp::checkUserInterrupt()` (which throws, so destructors unwind,
  rather than `R_CheckUserInterrupt()` which longjmps) now runs at the
  top of each BFGS iteration, at the top of each finite-difference
  Hessian column, and between 65536-element chunks of the `d`/`p`/`q`
  loops.

- **Exceptions can no longer escape a parallel region.** A throw
  crossing an OpenMP structured block terminates the process rather than
  unwinding (verified: SIGABRT, the outer `catch` never runs). One
  `arma::inv()` call in the cluster loop used the throwing overload; it
  now uses the boolean form like every other decomposition there. The
  cluster loop body, which performs about twenty Armadillo allocations
  per cluster, is wrapped so that a `std::bad_alloc` under memory
  pressure becomes a clean R error instead of an abort.

- **A badly scaled design converges.** The first BFGS trial step was 1.0
  in the units of the raw gradient. With a covariate scaled by 1e5 that
  overshoots by eight orders of magnitude and the 40 available halvings
  never reach a decrease, so the fit stopped at iteration 1 with a
  log-likelihood 155.5 units below what
  [`nlminb()`](https://rdrr.io/r/stats/nlminb.html) reaches on the
  identical objective. The step is now scaled while `H` is still the
  identity – which is also 22% *cheaper* on ordinary fits.

- **`RhpcBLASctl` (Suggests) is used to pin the BLAS** for the duration
  of a mixed fit and restore it afterwards, worth a further 34% at four
  threads.

- Parallel behaviour is now tested (`tests/testthat/test-parallel.R`):
  the distribution functions are bit-identical across thread counts,
  fits agree to floating-point noise, and interruptibility has a
  regression test.

- Removed the adaptive Gauss-Legendre quadrature, dead since the CDF
  became a closed form.

### Packaging

- `Depends: R (>= 4.0.0)` – no R 4.1 feature is used.
- Suggests now declares what is actually used: `microbenchmark` (the
  shipped benchmark script) and `MASS` (a mixed-model test).
  `simplexreg` is deliberately NOT declared – it is archived on CRAN, so
  declaring it would produce a permanent NOTE;
  `inst/benchmark/run_benchmark.R` guards it with
  [`requireNamespace()`](https://rdrr.io/r/base/ns-load.html) and says
  so in its header.
- Added `inst/CITATION`; fixed the test-coverage workflow (dropped the
  archived `simplexreg` reference); dropped the redundant
  `SystemRequirements: C++17`; fixed a 404 badge in the README.

## fastsimplexreg 0.2.3

Correctness release following a full surgical audit of `R/`, `src/`,
`tests/` and `vignettes/`. All estimates, log-likelihoods, deviances,
fitted values, residuals, densities and random draws are numerically
unchanged from 0.2.2; what changed is what the package does when
something is *wrong*.

### Breaking changes

- **The distribution functions now recycle instead of erroring.**
  [`dsimplex()`](https://evandeilton.github.io/fastsimplexreg/reference/simplex-distribution.md)
  and
  [`rsimplex()`](https://evandeilton.github.io/fastsimplexreg/reference/simplex-distribution.md)
  used to require `mu`/`phi` of length 1 or `length(x)`; they now
  recycle all arguments to their common maximum length, which is the
  base-R convention (`dnorm(1:3, mean = 1:2)`).
- **`NA` and invalid parameters no longer collapse to `0`.**
  [`dsimplex()`](https://evandeilton.github.io/fastsimplexreg/reference/simplex-distribution.md)
  returned `0` (`-Inf` on the log scale) for `NA` input and for an
  out-of-domain `mu`/`phi`, turning “missing” into “impossible”. `NA`
  now propagates as `NA`, `NaN` as `NaN`, and an out-of-domain parameter
  gives `NaN` with the canonical “NaNs produced” warning.
  [`rsimplex()`](https://evandeilton.github.io/fastsimplexreg/reference/simplex-distribution.md)
  gives `NaN` with “NAs produced” instead of aborting the call.

### Inference (critical)

- **Standard errors are never reported as a confident zero again.** With
  collinear or otherwise unidentified covariates the Hessian could be
  singular or indefinite, [`solve()`](https://rdrr.io/r/base/solve.html)
  would still succeed, and `sqrt(pmax(diag(vcov), 0))` turned the
  resulting negative variances into `SE = 0` – so
  [`summary()`](https://rdrr.io/r/base/summary.html) printed `z = Inf`
  and `p = 0` for a parameter the data cannot identify, with no warning.
  Inversion now goes through an eigen-decomposition of the equilibrated
  (correlation-scale) information matrix: only the strictly positive
  spectrum is inverted (a Moore-Penrose pseudo-inverse), directions of
  negative curvature are refused rather than turned into variances,
  every parameter loading on a discarded direction gets `SE = NA`, and
  the affected parameters are named in a warning. A merely
  ill-conditioned but full-rank fit keeps its (large, honest) standard
  errors and gets a separate warning.
- New components on the fitted object: `vcov_rank`, `vcov_pseudo`,
  `vcov_condition`, `vcov_eigenvalues` and `n_saturated`.
- **Saturation of the mean link is reported.** The likelihood path
  floors the mean at `1e-12`, which zeroes those observations’
  contribution to the score. That now raises a warning naming how many
  observations are affected.

### New

- [`psimplex()`](https://evandeilton.github.io/fastsimplexreg/reference/simplex-distribution.md)
  and
  [`qsimplex()`](https://evandeilton.github.io/fastsimplexreg/reference/simplex-distribution.md)
  complete the `d`/`p`/`q`/`r` family, with `lower.tail` and `log.p`.
  The CDF uses adaptive Gauss-Legendre quadrature with panels seeded
  around the mean (so a sharply peaked density is always resolved) and
  agrees with
  [`stats::integrate()`](https://rdrr.io/r/stats/integrate.html) to
  ~1e-15;
  [`qsimplex()`](https://evandeilton.github.io/fastsimplexreg/reference/simplex-distribution.md)
  inverts it by safeguarded Newton-bisection.

### Prediction

- [`predict()`](https://rdrr.io/r/stats/predict.html) now preserves the
  length of `newdata`: rows dropped for missingness come back as `NA`
  instead of silently yielding a shorter, unaligned vector. All model
  parts are built from the same set of complete rows, so they can no
  longer end up with different row counts.
- [`predict()`](https://rdrr.io/r/stats/predict.html) on a mixed fit
  **errors** when `newdata` lacks the grouping column. It previously
  returned population-level predictions under the label of conditional
  ones, silently. Rows whose group level was not seen in the fit still
  fall back to a zero random effect, but that substitution is now
  announced.
- Population-level prediction on a mixed fit stored with `model = FALSE`
  errors instead of resolving covariates in the caller’s environment,
  where it could return predictions of the wrong length built from
  unrelated objects.
- A variable required by the model but absent from `newdata` is an
  error.

### Numerical / C++

- The reporting path
  ([`simplex_linkinv()`](https://evandeilton.github.io/fastsimplexreg/reference/simplex_linkinv.md),
  [`predict()`](https://rdrr.io/r/stats/predict.html),
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html)) no longer
  applies the likelihood path’s `1e-12` floor. It is clamped only at the
  representable boundary, so fitted means keep their full dynamic range
  (`simplex_linkinv(-40, "logit")` is `4.2e-18`, not `1e-12`) while
  remaining strictly inside the open support `(0, 1)` that the density
  and the residual formulas require.
- `simplex_mixed_ranef_cpp()` no longer reads uninitialised memory: the
  Fisher fallback for the posterior covariance consumed a buffer that
  the observed loop could leave partially filled, and the return value
  of the link map was discarded.
- The mixed model’s inner mode solver no longer commits a step its line
  search rejected, and no longer compares against the objective at an
  inadmissible point. A warm-started mode inherited from an outer trial
  point that was later rejected is now discarded when the cold start
  beats it, removing a path dependence in the objective.
- [`VarCorr()`](https://rdrr.io/pkg/nlme/man/VarCorr.html)’s `sigma`
  argument is inert for a simplex mixed model; supplying anything other
  than `1` now warns instead of being silently ignored.

### Mixed-model convergence

- **Roughly a quarter of mixed fits were reported as optimiser failures
  when they had in fact converged**, and because standard errors are
  only computed at a converged fit, those perfectly good fits silently
  lost their inference. Measured over a grid of 240
  `(J, nj, nAGQ, seed)` combinations, 26.2% ended with code 2;
  restarting the optimiser from the reported stopping point gained a
  median of 1.6e-11 in log-likelihood – it was an optimum, not a
  failure. The failure rate is now 0% over the same grid.

  The rate depends on `nAGQ`, not on the number of groups: 55% at
  `nAGQ = 3`, 47% at 5, 3% at 7 and 0% at the default 11, while being
  flat in `J` (25-29%) and in cluster size. The cause is that the
  analytic score is the exact score of the *true* marginal likelihood
  (Fisher’s identity), not of its `nAGQ`-point quadrature approximation,
  so `grad_tol` is unreachable when `nAGQ` is small and the run ends on
  a line-search failure instead.

  Two changes: when a line search fails, the inverse-Hessian
  approximation is reset and the iteration is retried once from a clean
  steepest-descent direction (it goes stale, and this recovers real
  progress – it moved 59 of 180 stopping points, every one of them to a
  *higher* log-likelihood); and the soft-convergence test, which
  required a relative change below `rel_tol` (1e-9), now uses
  `sqrt(rel_tol)`, the scale at which these runs actually flatten out. A
  genuine failure – a run that cannot take a single step – is still
  reported as code 2.

  Cost: the fixed-effects path is unaffected in results and essentially
  unaffected in speed (+3.2% function evaluations, in 3 of 24 fits, all
  `neglog`; wall clock within noise). Mixed fits that already converged
  use about 50% more function evaluations, which buys the extra accuracy
  above; fits that used to abort naturally cost more now that they run
  to completion.

- Known limitation: supplying `start` equal to the optimum itself gives
  the optimiser no iteration history from which to judge stationarity,
  so it conservatively reports code 2 rather than risk labelling a stall
  as success.

- The optimiser trace (`trace = TRUE`) now prints the objective at full
  double precision and reports the relative change per iteration. At the
  previous six significant digits, successive iterations near the
  optimum printed identically – exactly the regime the trace exists to
  diagnose.

### Testing and documentation

- New `test-inference.R` and `test-distribution-conventions.R`, plus
  regression tests pinning the convergence-reporting behaviour in both
  directions.
- The rank cut-off used by the fail-safe covariance is derived from the
  accuracy of the finite-difference Hessian (`sqrt(eps)` on the
  correlation scale) rather than from machine epsilon. With the tighter
  cut the verdict depended on which BLAS computed the Hessian: an
  exactly collinear design was flagged on Linux and passed as full rank
  on Windows.
- The analytic score is now compared directly against `numDeriv` for all
  four mean links (previously validated only indirectly, and only for
  `logit`).
- The benchmark vignette’s accuracy table is generated from the shipped
  results instead of being hard-coded, so it cannot drift from the
  figures.
- [`simulate()`](https://rdrr.io/r/stats/simulate.html) now carries the
  row names of the model frame (the previous code read
  `names(object$y)`, which is always `NULL`).

## fastsimplexreg 0.2.2

- New `benchmark` vignette comparing `fastsimplexreg` with the CRAN
  packages `simplexreg` and `betareg` on accuracy, speed (scaling to n =
  5e5), the four mean links, real data, and the mixed model, with
  professional ggplot2 figures. Estimates are numerically identical to
  `simplexreg` while fitting ~5-13x faster. The full reproducible study
  ships in `inst/benchmark/`.

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
- [`dsimplex()`](https://evandeilton.github.io/fastsimplexreg/reference/simplex-distribution.md)
  and
  [`rsimplex()`](https://evandeilton.github.io/fastsimplexreg/reference/simplex-distribution.md)
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
