# Extractor Methods for Simplex Mixed-Model Fits

Standard extractor methods for objects of class `"simplex_fast_mixed"`
produced by
[`fastsimplexregmixed()`](https://evandeilton.github.io/fastsimplexreg/reference/fastsimplexregmixed.md).

- `coef`:

  Fixed-effect coefficients (`model = "all"`, `"mean"` or
  `"dispersion"`). Random effects are obtained with
  [`ranef()`](https://rdrr.io/pkg/nlme/man/random.effects.html).

- `vcov`:

  Covariance matrix of the estimated parameters `c(beta, gamma, omega)`.

- `logLik`:

  Maximised marginal log-likelihood, with `df = p + r + q(q+1)/2` and
  `nobs`.

- `nobs`:

  Number of observations.

- `ngrps`:

  Number of groups.

- `fitted`:

  Fitted means (conditional on the empirical-Bayes random effects) or
  fitted dispersions.

- `residuals`:

  Response, Pearson or deviance residuals, conditional on the
  empirical-Bayes random effects.

- `ranef`:

  Empirical-Bayes random-effect modes (a groups-by-`q` matrix); with
  `postVar = TRUE`, the posterior covariances are attached as the
  `"postVar"` attribute.

- `VarCorr`:

  The estimated random-effect covariance matrix \\\Sigma\\, with
  standard deviations and correlations.

- `confint`:

  Wald confidence intervals over the full parameter vector
  `c(beta, gamma, omega)`. Note that the intervals for the variance
  components are on the unconstrained log-Cholesky scale, where a Wald
  interval is defensible; on the variance scale it would not be, because
  the null lies on the boundary.

## Usage

``` r
# S3 method for class 'simplex_fast_mixed'
coef(object, model = c("all", "mean", "dispersion"), ...)

# S3 method for class 'simplex_fast_mixed'
vcov(object, ...)

# S3 method for class 'simplex_fast_mixed'
confint(object, parm, level = 0.95, ...)

# S3 method for class 'simplex_fast_mixed'
logLik(object, ...)

# S3 method for class 'simplex_fast_mixed'
nobs(object, ...)

# S3 method for class 'simplex_fast_mixed'
fitted(object, model = c("mean", "dispersion"), ...)

# S3 method for class 'simplex_fast_mixed'
residuals(object, type = c("quantile", "response", "pearson", "deviance"), ...)

# S3 method for class 'simplex_fast_mixed'
ranef(object, postVar = FALSE, ...)

# S3 method for class 'simplex_fast_mixed'
VarCorr(x, sigma = 1, ...)

# S3 method for class 'VarCorr.simplex_fast_mixed'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- object, x:

  A fitted `"simplex_fast_mixed"` object.

- model:

  For `coef`, one of `"all"`, `"mean"` or `"dispersion"`; for `fitted`,
  one of `"mean"` or `"dispersion"`.

- ...:

  Additional arguments, currently ignored.

- parm:

  For `confint`, which parameters to report: numeric positions or names,
  over the full vector `c(beta, gamma, omega)`. Defaults to all.

- level:

  For `confint`, the confidence level.

- type:

  For `residuals`, one of `"response"`, `"pearson"` or `"deviance"`.

- postVar:

  For `ranef`, logical; attach posterior covariances.

- sigma:

  For `VarCorr`, present only to match the signature of
  [`nlme::VarCorr()`](https://rdrr.io/pkg/nlme/man/VarCorr.html). A
  simplex mixed model has no residual scale parameter, so the argument
  rescales nothing; supplying anything other than `1` raises a warning
  and is ignored.

- digits:

  For the `VarCorr` print method, the number of significant digits to
  display.

## Value

`coef`, `fitted` and `residuals` return numeric vectors; `vcov` returns
a matrix; `ranef` returns a matrix; `VarCorr` returns the covariance
matrix with `stddev`/`correlation` attributes; `logLik` returns a
`"logLik"` object.

## See also

[`fastsimplexregmixed()`](https://evandeilton.github.io/fastsimplexreg/reference/fastsimplexregmixed.md)

## Examples

``` r
set.seed(1)
J <- 40; nj <- 8; n <- J * nj
dat <- data.frame(g = factor(rep(seq_len(J), each = nj)), x1 = rnorm(n))
b <- rnorm(J, 0, 0.7)[dat$g]
dat$y <- rsimplex(n, simplex_linkinv(0.3 - 0.6 * dat$x1 + b, "logit"), 1)
fit <- fastsimplexregmixed(y ~ x1, random = ~ 1 | g, data = dat,
                           nAGQ = 7, n_threads = 1)

coef(fit)
#> (Intercept)          x1 (Intercept) 
#>   0.3426502  -0.6330471   0.1200507 
VarCorr(fit)
#> Random effects covariance (group: g)
#>             Variance Std.Dev.
#> (Intercept)   0.5001   0.7072
head(ranef(fit))
#>   (Intercept)
#> 1  0.88838598
#> 2  0.89201113
#> 3  0.61477597
#> 4 -1.40691453
#> 5  0.02536362
#> 6  0.28389508
confint(fit)
#>                         2.5 %     97.5 %
#> (Intercept)        0.11826112  0.5670393
#> x1                -0.68401645 -0.5820778
#> (phi)_(Intercept) -0.04579166  0.2858931
#> logsd.(Intercept) -0.57574656 -0.1171643
ngrps(fit)
#> [1] 40
```
