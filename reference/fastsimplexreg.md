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

  A multi-part formula, for example `y ~ x1 + x2 | z1 + z2`.

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

  Integer number of OpenMP threads. Use `0` to request all threads
  available to the backend.

- inference:

  Logical; if `TRUE`, computes the Hessian, the variance-covariance
  matrix and the standard errors.

- hessian_rel_step:

  Numeric; the initial relative step for the Hessian, obtained by
  central differences of the analytic gradient.

- trace:

  Logical; if `TRUE`, prints optimiser progress.

- subset:

  Optional vector specifying a subset of observations.

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
prediction.

## References

Barndorff-Nielsen, O. E. and Jorgensen, B. (1991). Some parametric
models on the simplex. *Journal of Multivariate Analysis*, **39**(1),
106–116.

Zhang, P., Qiu, Z. and Shi, C. (2016). simplexreg: An R Package for
Regression Analysis of Proportional Data Using the Simplex Distribution.
*Journal of Statistical Software*, **71**(11), 1–21.

## See also

[`dsimplex()`](https://evandeilton.github.io/fastsimplexreg/reference/dsimplex.md),
[`rsimplex()`](https://evandeilton.github.io/fastsimplexreg/reference/rsimplex.md),
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
#> Pearson residuals:
#>      Min       1Q   Median       3Q      Max 
#> -2.46935 -0.70189 -0.06618  0.67738  2.87579 
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
#> Deviance:   500 | Observations: 500 | Iterations: 18 
#> Convergence: 0 - Converged: relative objective tolerance satisfied. 
coef(fit)
#> (Intercept)          x1          x2 (Intercept)          z1 
#>  -0.3969817   0.8066782  -0.5257617  -1.0403689   0.6303906 
head(predict(fit, type = "both"))
#>          mu       phi
#> 1 0.2996206 0.9318746
#> 2 0.3583206 0.3297141
#> 3 0.7027430 0.4877547
#> 4 0.2961153 0.4043430
#> 5 0.3060928 0.3142088
#> 6 0.7284008 0.3275010

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
#> Pearson residuals:
#>      Min       1Q   Median       3Q      Max 
#> -2.39081 -0.62295  0.24243  0.43805  1.48447 
#> 
#> Coefficients (mean model with logit link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  1.37697    0.15352   8.969  < 2e-16 ***
#> dyslexia    -0.97657    0.15485  -6.307 2.85e-10 ***
#> iq          -0.04369    0.07130  -0.613     0.54    
#> 
#> Coefficients (dispersion model with log link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)   1.4242     0.2152   6.617 3.66e-11 ***
#> dyslexia     -2.6917     0.2162 -12.450  < 2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Log-likelihood: 68.01 | AIC:  -126 | BIC: -117.1 
#> Deviance:    44 | Observations: 44 | Iterations: 15 
#> Convergence: 0 - Converged: relative objective tolerance satisfied. 
```
