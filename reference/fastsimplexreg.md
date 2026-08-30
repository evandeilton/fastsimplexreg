# Fit a Fast Simplex Regression with Variable Dispersion

Fits, by maximum likelihood, a simplex regression model with separate
submodels for the mean and the dispersion. The interface uses the
multi-part formulas of the Formula package:

`y ~ x1 + x2 | z1 + z2`

The first right-hand side component models the mean \\\mu\\; the second
component models the dispersion \\\phi\\. When the second component is
omitted, as in `y ~ x1 + x2`, the dispersion is constant (equivalent to
`| 1`).

The mean supports the `logit`, `probit`, `cloglog` and `neglog` links;
the dispersion uses a log link. The log-likelihood, the analytic score,
the link inverses and the BFGS optimiser run entirely in C++.
Matrix-vector products use Armadillo/BLAS and the per-observation loop
may use OpenMP.

## Usage

``` r
fastsimplexreg(
  formula,
  data,
  link = c("logit", "probit", "cloglog", "neglog"),
  start = NULL,
  maxit = 300L,
  rel_tol = 1e-09,
  grad_tol = 1e-06,
  n_threads = 1L,
  inference = TRUE,
  information = c("observed", "expected"),
  hessian_rel_step = 1e-05,
  trace = FALSE,
  subset = NULL,
  na.action = stats::na.omit,
  model = TRUE,
  x = FALSE,
  y = TRUE
)
```

## Arguments

- formula:

  A multi-part formula, for example `y ~ x1 + x2 | z1 + z2`. Each part
  may carry its own [`offset()`](https://rdrr.io/r/stats/offset.html)
  term: in `y ~ x1 + offset(a) | z1 + offset(b)`, `a` is added to the
  mean linear predictor and `b` to the dispersion linear predictor, each
  on its own link scale, and neither is estimated. The two are kept
  apart – unlike
  [`stats::model.offset()`](https://rdrr.io/r/stats/model.extract.html),
  which sums the offsets of every part of a multi-part formula – and are
  rebuilt from `newdata` in
  [`predict.simplex_fast()`](https://evandeilton.github.io/fastsimplexreg/reference/predict.simplex_fast.md),
  so `newdata` must supply every variable an offset uses.

- data:

  A `data.frame` containing the response and covariates.

- link:

  Character string selecting the mean link: `"logit"`, `"probit"`,
  `"cloglog"` or `"neglog"`.

- start:

  Optional numeric starting vector `c(beta, gamma)`. When `NULL`, fast
  link-specific starting values are used.

- maxit:

  Integer; the maximum number of BFGS iterations.

- rel_tol:

  Numeric; relative tolerance on the objective function.

- grad_tol:

  Numeric; tolerance on the infinity norm of the gradient.

- n_threads:

  Integer number of OpenMP threads; the loop over observations is
  parallelised. Use `0` to request all threads available to the backend.
  The linear predictors are accumulated inside that loop rather than
  through BLAS, so performance does not depend on how your BLAS is
  configured. Measured end-to-end speed-ups on a 24-core machine:

  |                 |               |       |       |        |
  |-----------------|---------------|-------|-------|--------|
  | **problem**     | **2 threads** | **4** | **8** | **16** |
  | n = 1e5, p = 5  | 1.4x          | 2.2x  | 2.1x  | 2.2x   |
  | n = 1e6, p = 10 | 1.9x          | 3.2x  | 5.0x  | 5.2x   |
  | n = 5e6, p = 10 | 1.1x          | 3.3x  | 5.1x  | 5.1x   |

  Threading pays from roughly `n = 1e5` upwards and saturates near 5x,
  at which point the loop is memory-bound. The default is `1L` so that
  the package never oversubscribes a machine it does not own.

  Results are reproducible for a fixed `n_threads`. Across DIFFERENT
  thread counts the per-thread accumulators are summed in a different
  order, which perturbs the objective at the rounding level (~1e-14
  relative); the line search then follows a slightly different path, so
  coefficients can differ by around `1e-8` relative. Compare with
  [`all.equal()`](https://rdrr.io/r/base/all.equal.html), not
  [`identical()`](https://rdrr.io/r/base/identical.html).

- inference:

  Logical; if `TRUE`, computes the information matrix, the
  variance-covariance matrix and the standard errors.

- information:

  Character; which information matrix to invert for the standard errors.
  `"observed"` (default) uses the observed information, the Hessian of
  the negative log-likelihood obtained by central differences of the
  analytic score. `"expected"` uses the exact Fisher information, which
  for the simplex is available in closed form and is block diagonal in
  \\(\beta, \gamma)\\: it needs no finite differencing, is positive
  definite by construction, and is roughly twenty times cheaper. The two
  agree asymptotically and, at \\n = 4000\\, to within 0.4\\ stays
  `"observed"` because Efron and Hinkley (1978) argue it is the better
  variance estimator for conditional inference; `"expected"` is the more
  robust choice when the observed information is ill-conditioned.

- hessian_rel_step:

  Numeric; the initial relative step for the Hessian, obtained by
  central differences of the analytic gradient.

- trace:

  Logical; if `TRUE`, prints optimiser progress.

- subset:

  Optional expression selecting a subset of observations, evaluated
  inside `data` as in [`stats::lm()`](https://rdrr.io/r/stats/lm.html) –
  for example `subset = x1 > 0`. A plain index, logical or row-name
  vector also works. An `NA` in a logical subset drops that row.

- na.action:

  A function indicating how to handle missing values.

- model:

  Logical; if `TRUE`, stores the model frame in the fitted object.

- x:

  Logical; if `TRUE`, stores the design matrices `X` and `Z`.

- y:

  Logical; if `TRUE`, stores the response in the fitted object.

## Value

An object of S3 class `"simplex_fast"`: a list whose main components are
`coefficients` (a list with `mean` and `dispersion` estimates), `par`
(the full coefficient vector), `standard_errors`, `vcov`,
`fitted.values` (fitted means), `dispersion.values` (fitted
dispersions), `linear.predictors`, `residuals` (response residuals),
`logLik`, `AIC`, `BIC`, `nobs`, `df.residual`, `convergence`, `message`,
`iterations` and the stored `terms`/`design` metadata used for
prediction. Inference diagnostics are also stored: `vcov_rank` (rank of
the observed information matrix), `vcov_pseudo` (`TRUE` when a
Moore-Penrose pseudo-inverse was required because the matrix was rank
deficient or indefinite), `vcov_eigenvalues`, and `n_saturated`
(observations whose fitted mean hit the numerical boundary of the
likelihood path). Standard errors of parameters that the data do not
identify are `NA`, never `0`.

## References

Barndorff-Nielsen, O. E. and Jorgensen, B. (1991). Some parametric
models on the simplex. *Journal of Multivariate Analysis*, **39**(1),
106–116.
[doi:10.1016/0047-259X(91)90008-P](https://doi.org/10.1016/0047-259X%2891%2990008-P)

Zhang, P., Qiu, Z. and Shi, C. (2016). simplexreg: An R Package for
Regression Analysis of Proportional Data Using the Simplex Distribution.
*Journal of Statistical Software*, **71**(11), 1–21.
[doi:10.18637/jss.v071.i11](https://doi.org/10.18637/jss.v071.i11)

Efron, B. and Hinkley, D. V. (1978). Assessing the accuracy of the
maximum likelihood estimator: observed versus expected Fisher
information. *Biometrika*, **65**(3), 457–483.
[doi:10.1093/biomet/65.3.457](https://doi.org/10.1093/biomet/65.3.457)

## See also

[`dsimplex()`](https://evandeilton.github.io/fastsimplexreg/reference/simplex-distribution.md),
[`rsimplex()`](https://evandeilton.github.io/fastsimplexreg/reference/simplex-distribution.md),
[`simplex_linkinv()`](https://evandeilton.github.io/fastsimplexreg/reference/simplex_linkinv.md),
[`predict.simplex_fast()`](https://evandeilton.github.io/fastsimplexreg/reference/predict.simplex_fast.md),
[`summary.simplex_fast()`](https://evandeilton.github.io/fastsimplexreg/reference/summary.simplex_fast.md)

## Examples

``` r
# Simulated data with variable dispersion.
set.seed(123)
n <- 500
dat <- data.frame(x1 = rnorm(n), x2 = rbinom(n, 1, 0.4), z1 = rnorm(n))
mu <- simplex_linkinv(-0.4 + 0.8 * dat$x1 - 0.5 * dat$x2, link = "logit")
phi <- exp(-1 + 0.6 * dat$z1)
dat$y <- rsimplex(n, mu, phi)

fit <- fastsimplexreg(y ~ x1 + x2 | z1, data = dat, link = "logit",
                   n_threads = 1L)
summary(fit)
#> 
#> Call:
#> fastsimplexreg(formula = y ~ x1 + x2 | z1, data = dat, link = "logit", 
#>     n_threads = 1L)
#> 
#> Quantile residuals:
#>      Min       1Q   Median       3Q      Max 
#> -2.72184 -0.67703 -0.02365  0.68464  2.62614 
#> 
#> Coefficients (mean model with logit link):
#>              Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) -0.396982   0.013501  -29.40   <2e-16 ***
#> x1           0.806678   0.008785   91.83   <2e-16 ***
#> x2          -0.525762   0.021152  -24.86   <2e-16 ***
#> 
#> Coefficients (dispersion model with log link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) -1.04037    0.06326 -16.445   <2e-16 ***
#> z1           0.63039    0.06345   9.936   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Log-likelihood: 774.8 | AIC: -1540 | BIC: -1519 
#> Deviance: 219.9 | Observations: 500 | Iterations: 20 
#> Convergence: 0 - Converged: relative objective tolerance satisfied. 
coef(fit)
#>       (Intercept)                x1                x2 (phi)_(Intercept) 
#>        -0.3969818         0.8066782        -0.5257615        -1.0403692 
#>          (phi)_z1 
#>         0.6303906 
head(predict(fit, type = "both"))
#>          mu       phi
#> 1 0.2996206 0.9318743
#> 2 0.3583205 0.3297140
#> 3 0.7027430 0.4877546
#> 4 0.2961153 0.4043428
#> 5 0.3060928 0.3142087
#> 6 0.7284007 0.3275009

# Real data: reading accuracy from the 'betareg' package.
if (requireNamespace("betareg", quietly = TRUE)) {
  data("ReadingSkills", package = "betareg")
  rs <- fastsimplexreg(accuracy ~ dyslexia + iq | dyslexia,
                       data = ReadingSkills, link = "logit")
  summary(rs)
}
#> 
#> Call:
#> fastsimplexreg(formula = accuracy ~ dyslexia + iq | dyslexia, 
#>     data = ReadingSkills, link = "logit")
#> 
#> Quantile residuals:
#>      Min       1Q   Median       3Q      Max 
#> -2.37009 -0.80217  0.14155  0.90867  1.55585 
#> 
#> Coefficients (mean model with logit link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  1.37696    0.15352   8.969  < 2e-16 ***
#> dyslexia    -0.97656    0.15485  -6.307 2.85e-10 ***
#> iq          -0.04369    0.07130  -0.613     0.54    
#> 
#> Coefficients (dispersion model with log link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)   1.4242     0.2152   6.617 3.66e-11 ***
#> dyslexia     -2.6918     0.2162 -12.450  < 2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Log-likelihood: 68.01 | AIC:  -126 | BIC: -117.1 
#> Deviance:  1538 | Observations: 44 | Iterations: 15 
#> Convergence: 0 - Converged: relative objective tolerance satisfied. 
```
